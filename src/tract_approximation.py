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

    
    # BVAP >= 50% FOR DISTRICTS 1 AND 6
    codes = ['P0030004'] # Black or African American alone 
    codes += ['P0030011','P0030016','P0030017','P0030018','P0030019'] # Black or African American (among 2 races)
    codes += ['P0030027','P0030028','P0030029','P0030030','P0030037','P0030038','P0030039','P0030040','P0030041','P0030042'] # 3
    codes += ['P0030048','P0030049','P0030050','P0030051','P0030052','P0030053','P0030058','P0030059','P0030060','P0030061'] # 4
    codes += ['P0030064','P0030065', 'P0030066','P0030067','P0030069'] # 5
    codes += ['P0030071'] # 6

    for i in G.nodes:
        G.nodes[i]['VAP'] = G.nodes[i]['P0030001'] 
        G.nodes[i]['MVAP'] = sum( G.nodes[i][code] for code in codes )

    # Impose mvap >= 0.5 * vap for districts 2 and 7 (i.e., 1 and 6 in python)
    m.addConstr( gp.quicksum( G.nodes[i]['MVAP'] * x[i,1] for i in G.nodes ) >= 0.50 * gp.quicksum( G.nodes[i]['VAP'] * x[i,1] for i in G.nodes ) )
    m.addConstr( gp.quicksum( G.nodes[i]['MVAP'] * x[i,6] for i in G.nodes ) >= 0.50 * gp.quicksum( G.nodes[i]['VAP'] * x[i,6] for i in G.nodes ) )

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
def distance_to_vertex_set(G, S, cutoff=None):
    if cutoff==None:
        cutoff = G.number_of_nodes()
    
    dist = dict()
    touched = { i : False for i in G.nodes }
    child = S
    for s in S:
        touched[s] = True
    
    for h in range(1+cutoff):
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
    