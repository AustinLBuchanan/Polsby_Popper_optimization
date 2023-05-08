import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from mip_contiguity import find_fischetti_separator

# read block assignment file and store in dictionary:
#    assignment = { block_geoid : district_number }
def read_block_assignment_file(filepath, filename):
    assignment = dict()
    import pandas
    csvFile = pandas.read_csv( filepath + filename, skipinitialspace = True )
    for index, row in csvFile.iterrows():
        geoid = str( row[0] ) 
        if len(geoid) < 15: # fix issue with leading zeros
            geoid = '0' + geoid
        district_string = str( row[1] ) 
        if district_string in { 'ZZ', 'ZZZ' }:
            print("Warning: geoid =",geoid,"is unassigned.")
            assignment[geoid] = -1
            continue
        assignment[geoid] = int(district_string)
    return assignment


def get_tract_approximation(G, GB, L, U, k, block_assignment):
    
    ####################################
    # print properties of initial plan
    ####################################
    
    # retrieve districts
    districts = [ list() for j in range(k) ]
    for b in GB.nodes:
        g = GB.nodes[b]['GEOID20']
        j = block_assignment[g] - 1
        districts[j].append(b)
    
    # print their populations and connectivity status
    district_populations = [ sum( GB.nodes[i]['TOTPOP'] for i in districts[j] ) for j in range(k) ]
    connected = [ nx.is_connected( GB.subgraph(districts[j]) ) for j in range(k) ]
    print("Initial block-level plan:")
    print("# \t population \t connected?")
    for j in range(k):
        print(j,'\t',district_populations[j],'\t',connected[j])
        if not connected[j]:
            print("Component sizes:",end=" ")
            for comp in nx.connected_components( GB.subgraph(districts[j]) ):
                print(len(comp), end=', ')
            print("")
    print("min, max pop:", min(district_populations), max(district_populations) )
    print("total deviation =", max(district_populations) - min(district_populations), '\n' )
    
    ####################################
    # helpful dictionaries
    ####################################
    
    # gtn = { tract_geoid : tract_node }
    gtn = { G.nodes[i]['GEOID20'] : i for i in G.nodes }

    # btt = { block_node : its_tract_node }
    btt = { i : gtn[ GB.nodes[i]['GEOID20'][0:11] ] for i in GB.nodes }
    
    
    ####################################
    # build gurobi model
    ####################################
    
    m = gp.Model()
    x = m.addVars(G.nodes, k, vtype=GRB.BINARY)

    # assignment and population balance constraints
    m.addConstrs( gp.quicksum( x[i,j] for j in range(k) ) == 1 for i in G.nodes )
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in G.nodes ) >= L for j in range(k) )
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in G.nodes ) <= U for j in range(k) )

    # what is the "cost" to transport tract i to the block-level plan's district j?
    cost = { (i,j) : 0 for i in G.nodes for j in range(k) }        

    for j in range(k):
        dist = distance_to_vertex_set(GB, districts[j])
        for b in GB.nodes:
            i = btt[b]
            cost[i,j] += dist[b] * ( 1 + GB.nodes[b]['TOTPOP'] )

    # minimize the transportation cost (graph walk distance for block population reassignments)
    m.setObjective( gp.quicksum( cost[i,j] * x[i,j] for i in G.nodes for j in range(k)), GRB.MINIMIZE)

    # parameter settings
    m.Params.MIPGap = 0
    m.Params.IntFeasTol = 1e-7
    m.Params.FeasibilityTol = 1e-7
    m.Params.LazyConstraints = 1

    # callback info
    DG = nx.DiGraph(G)
    DG._k = k
    m._x = x
    m._DG = DG
    m._callback = cut_callback
    m._numCallbacks = 0
    m._numLazyCuts = 0

    # solve gurobi model
    m.optimize(m._callback)
    if m.solCount == 0:
        print("No solution found. Exiting.")
        return
    
    ####################################
    # print properties of output plan
    ####################################
    
    # retrieve districts
    districts = [ [ i for i in G.nodes if x[i,j].x > 0.5 ] for j in range(k) ]
    
    # print their populations and connectivity status
    district_populations = [ sum( G.nodes[i]['TOTPOP'] for i in districts[j] ) for j in range(k) ]
    connected = [ nx.is_connected( G.subgraph(districts[j]) ) for j in range(k) ]
    print("Final tract-level approximate plan:")
    print("# \t population \t connected?")
    for j in range(k):
        print(j,'\t',district_populations[j],'\t',connected[j])
        if not connected[j]:
            print("Component sizes:")
            for comp in nx.connected_components( G.subgraph(districts[j]) ):
                print(len(comp), end=', ')
            print("")
    print("min, max pop:", min(district_populations), max(district_populations) )
    print("total deviation =", max(district_populations) - min(district_populations) )

    return districts

    
# dist[i] = the distance from node i to the set of nodes S
# dist[i]=0 if i in S
# dist[i]=1 if i in N(S)
# dist[i]=2 if i in N^2(S)
# ...
# dist[i]=None if no path from i to S
#
def distance_to_vertex_set(G, S):
    dist = { i : None for i in G.nodes }
    touched = { i : False for i in G.nodes }
    child = S
    for s in S:
        touched[s] = True
    for h in range(G.number_of_nodes()):
        parent = child
        child = list()
        for p in parent:
            dist[p] = h
            for n in G.neighbors(p):
                if not touched[n]:
                    child.append(n)
                    touched[n] = True
    return dist
    

def cut_callback(m, where):
    if where == GRB.Callback.MIPSOL:
        m._numCallbacks += 1 
        DG = m._DG
        xval = m.cbGetSolution(m._x)

        for j in range(DG._k):
            
            # vertices assigned to this district (label j)
            S = [ v for v in DG.nodes if xval[v,j] > 0.5 ]
            
            # what shall we deem as the "root" of this district? call it b
            b = None
            
            # for each component that doesn't contain b, add a cut
            for component in sorted( nx.strongly_connected_components( DG.subgraph(S) ), key=len, reverse=True ):

                # what is the maximum population node in this component?
                maxp = max( DG.nodes[v]['TOTPOP'] for v in component)
                mpv = [ v for v in component if DG.nodes[v]['TOTPOP'] == maxp ][0]

                # if no root 'b' has been selected yet, pick one
                if b is None:
                    # find some vertex "b" that has largest population in this component
                    b = mpv
                    continue

                # find some vertex "a" that has largest population in this component
                a = mpv

                # get minimal a,b-separator
                C = find_fischetti_separator(DG, component, b)

                # add lazy cut
                m.cbLazy( m._x[a,j] + m._x[b,j] <= 1 + gp.quicksum( m._x[c,j] for c in C ) )
                m._numLazyCuts += 1
                
    return
    