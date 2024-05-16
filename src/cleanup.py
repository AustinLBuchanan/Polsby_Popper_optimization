import math
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from census_codes import get_census_codes
from mip_contiguity import find_fischetti_separator 
from coarsen import hop_coarsen, distance_to_vertex_set

# In the spirit of "ILP-based local search for graph partitioning"
#    by A Henzinger, A Noe, C Schulz - Journal of Experimental Algorithmics, 2020
#    https://scholar.google.com/scholar?cluster=9173497710013156715&hl=en&as_sdt=0,37

def mip_local_search(G, initial_plan, h=1, max_iterations=10, minority=None, preserve_splits=True, 
                     multidistrict_sizes=None, guess_multidistrict_sizes=True, time_limit=600, verbose=False):
    
    print(f"Applying MIP-based local search to improve the MIP warm start (hop parameter: {h})...")
    new_plan = initial_plan.copy()
    new_obj = 100 # arbitrary "big" number
    grb_time = 0
    
    # sometimes, the initial_plan contains multidistricts (e.g., with twice the requisite population).
    # rather than explicitly pass their sizes, we can usually infer them.
    if multidistrict_sizes is None and guess_multidistrict_sizes:
        multidistrict_sizes = list()
        for district in initial_plan:
            population = sum( G.nodes[i]['TOTPOP'] for i in district )
            size_lb = math.ceil( population / G._U )
            size_ub = math.floor( population / G._L )
            assert size_lb == size_ub, "Failed to guess a multidistrict_size"
            multidistrict_sizes.append( size_lb )
        if not all(el==1 for el in multidistrict_sizes):
            print("Guessed the multidistrict_sizes =",multidistrict_sizes)
    
    # LOCAL SEARCH
    print("iter \t\t obj \t\t time")
    for iteration in range(1,max_iterations+1):
        
        old_plan = new_plan.copy()
        old_obj = new_obj
        ( GH, coarsen_partition, coarsened_old_plan ) = hop_coarsen(G, old_plan, h, preserve_splits )
        GH._k = G._k
        GH._L = G._L
        GH._U = G._U
        
        (coarsened_new_plan, new_obj, runtime) = best_neighboring_plan(GH, coarsened_old_plan, h=h, minority=minority, preserve_splits=preserve_splits, multidistrict_sizes=multidistrict_sizes,time_limit=time_limit, verbose=verbose)
        
        new_plan = list()
        for coarsened_district in coarsened_new_plan:
            district = list()
            for i in coarsened_district:
                district += coarsen_partition[i]
            new_plan.append(district)
        
        grb_time += runtime
        print(iteration,'\t','{0:.8f}'.format(new_obj),'\t','{0:.2f}'.format(runtime))
        
        # because of floating point numerical issues, objectives may not be monotone, especially near the end
        if new_obj >= old_obj - 1e-8 or new_plan == old_plan:
            break
    
    return new_plan


# G = input graph (with attributes L, U, k)
# initial_plan = starting point for local search (list of lists)
# h = # hops permitted
#
def best_neighboring_plan(G, initial_plan, h=1, minority=None, preserve_splits=False, 
                          time_limit=600, verbose=False, multidistrict_sizes=None):

    k = len(initial_plan)
    if multidistrict_sizes is None:
        assert G._k == k, "Implied number of districts does not match initial_plan"
    else:
        assert G._k == sum(multidistrict_sizes), "Implied number of districts does not match multidistrict_sizes"
    
    # create a (sparse) directed "assignment graph" AG 
    # that has vertices V = V(G) + {0, 1, ..., k-1}  (which may overlap)
    # and edge (i,j) if node i can be assigned to district j in {0, 1, ..., k-1}
    nodes = list( set(G.nodes) | set(range(k)) )
    AG = nx.DiGraph()
    AG.add_nodes_from( nodes )
    
    # add edges to AG
    for j in range(k):
        district = initial_plan[j]
        dist = distance_to_vertex_set(G, district, cutoff=h)
        AG.add_edges_from( [ (i,j) for i in dist.keys() ] )
        
    # initialize MIP
    m = gp.Model()
    if not verbose:
        m.Params.OutputFlag = 0
    m.Params.TimeLimit = time_limit
        
    x = m.addVars(AG.edges, vtype=GRB.BINARY)

    # assignment constraints
    m.addConstrs( gp.quicksum( x[i,j] for j in AG.neighbors(i) ) == 1 for i in G.nodes )

    # population balance
    ms = [ 1 for j in range(k) ] if multidistrict_sizes is None else multidistrict_sizes
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in AG.predecessors(j) ) >= ms[j] * G._L for j in range(k) )
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in AG.predecessors(j) ) <= ms[j] * G._U for j in range(k) )

    # minority VAP >= 50% for Gingles districts
    if minority is not None:
        try:
            u = nx.utils.arbitrary_element(G.nodes)
            G.nodes[u]['MVAP']
            G.nodes[u]['VAP']
        except:
            codes = get_census_codes(minority)
            for i in G.nodes:
                G.nodes[i]['VAP'] = G.nodes[i]['P0030001']
                G.nodes[i]['MVAP'] = sum( G.nodes[i][code] for code in codes )

        gingles_count = 0
        if verbose:
            print("Current Gingles districts:")
        for j in range(k):
            vap = sum( G.nodes[i]['VAP'] for i in initial_plan[j] )
            mvap = sum( G.nodes[i]['MVAP'] for i in initial_plan[j] )
            if mvap >= 0.5 * vap:
                gingles_count += 1
                # multiply constraint by 2 for integer coefficients
                coef = { i : 2 * G.nodes[i]['MVAP'] - G.nodes[i]['VAP'] for i in AG.predecessors(j) }
                m.addConstr( gp.quicksum( coef[i] * x[i,j] for i in coef.keys() ) >= 0 )
                if verbose:
                    print(f"district {j} which has {minority} VAP {round(100*mvap/vap,2)}%")
                
        if verbose:
            print("Total Gingles districts:",gingles_count)

    # fix nodes based on preserve_splits status
    if preserve_splits:
        counties = { G.nodes[i]['GEOID20'][0:5] for i in G.nodes }
        support = { c : list() for c in counties }
        
        for j in range(k):
            for i in initial_plan[j]:
                c = G.nodes[i]['GEOID20'][0:5]
                if j not in support[c]:
                    support[c].append(j)
        
        nonsupport = { c : set(range(k)) - set(support[c]) for c in counties }
        
        for i in G.nodes:
            c = G.nodes[i]['GEOID20'][0:5]
            for j in nonsupport[c]:
                if j in AG.neighbors(i):
                    x[i,j].UB = 0
    
    m.update()
    
    # fix nodes that have only one choice
    num_nodes_fixed = 0
    for i in G.nodes:
        choices = [ j for j in AG.neighbors(i) if x[i,j].UB > 0.5 ]
        if len(choices) == 1:
            j = choices[0]
            x[i,j].LB = 1
            num_nodes_fixed += 1
            
    m.update()
    
    if verbose:
        print("this many nodes fixed:",num_nodes_fixed)
        print("out of a total:",G.number_of_nodes())
        print("remaining nodes:",G.number_of_nodes()-num_nodes_fixed)
    
    ###########################################
    # ADD CUT EDGE VARS/CONSTRAINTS
    ###########################################
    
    # Interpretation: y[u,v,j]=1  <=> directed edge (u,v) is cut because u->j but v!->j 
    
    # Implementation: we create y[u,v,j] variable 
    #          <=> 
    # (u,v) is directed edge and x[u,j] was created (i.e., if j in AG.neighbors(u) )
    #
    triples = list()
    DG = nx.DiGraph(G)
    for u,v in DG.edges:
        triples += [ (u,v,j) for j in AG.neighbors(u) ]
    
    if verbose:
        print("this many triples (y cut edge vars) =",len(triples))
    y = m.addVars(triples, vtype=GRB.BINARY)
    
    # add cut edge constraints
    for u,v in DG.edges:
        for j in AG.neighbors(u): # x[u,j] exists
            if j in AG.neighbors(v): # x[v,j] exists
                m.addConstr( x[u,j] - x[v,j] <= y[u,v,j] )
            else: # x[v,j] does not exist (assumed to be zero)
                m.addConstr( x[u,j] == y[u,v,j] )
                
    # fix y vars when safe
    for u,v in DG.edges:
        for j in AG.neighbors(u):
            if j in AG.neighbors(v): # x[v,j] exists
                # if x[u,j]=0 or x[v,j]=1 then y[u,v,j]=0
                if x[u,j].UB < 0.5 or x[v,j].LB > 0.5:
                    y[u,v,j].UB = 0
                # if x[u,j]=1 and x[v,j]=0 then y[u,v,j]=1
                if x[u,j].LB > 0.5 and x[v,j].UB < 0.5:
                    y[u,v,j].LB = 1
            else: # x[v,j] does not exist
                if x[u,j].UB < 0.5:
                    y[u,v,j].UB = 0
                if x[u,j].LB > 0.5:
                    y[u,v,j].LB = 1
                
    # the following undirected cut edge vars are not strictly necessary, but seem to help computationally
    is_cut = m.addVars( G.edges, vtype=GRB.BINARY )
    m.addConstrs( is_cut[u,v] == gp.quicksum( y[u,v,j] for j in AG.neighbors(u) ) for u,v in G.edges )
    m.addConstrs( is_cut[u,v] == gp.quicksum( y[v,u,j] for j in AG.neighbors(v) ) for u,v in G.edges )
                
    ###########################################
    # ADD INVERSE POLSBY-POPPER OBJECTIVE
    ###########################################
    
    # z[j] is inverse Polsby-Popper score for district j
    z = m.addVars(k, name='z')

    # objective is to minimize average of inverse Polsby-Popper scores
    m.setObjective( ( 1.0 / k ) * gp.quicksum( z[j] for j in range(k) ), GRB.MINIMIZE )
    
    # A[j] = area of district j
    A = m.addVars(k, name='A')

    # P[j] = perimeter of district j
    P = m.addVars(k, name='P')

    # add SOCP constraints relating inverse Polsby-Popper score z[j] to area and perimeter
    m.addConstrs( P[j] * P[j] <= 4 * math.pi * A[j] * z[j] for j in range(k) )

    # add constraints on areas A[j] 
    m.addConstrs( A[j] == gp.quicksum( G.nodes[i]['area'] * x[i,j] for i in AG.predecessors(j) ) for j in range(k) )

    # add constraints on perimeters P[j]
    for j in range(k):
        m.addConstr( P[j] == gp.quicksum( G.edges[u,v]['shared_perim'] * y[u,v,j] for u in AG.predecessors(j) for v in G.neighbors(u) )
                 + gp.quicksum( G.nodes[i]['boundary_perim'] * x[i,j] for i in AG.predecessors(j) if G.nodes[i]['boundary_node'] ) )
    
    # add contiguity constraints via callback
    m.Params.LazyConstraints = 1
    m._x = x
    m._k = k
    m._DG = DG
    m._AG = AG
    m._callback = AG_contiguity_callback
    
    # inject warm start (the initial plan itself)
    for j in range(k):
        for i in initial_plan[j]:
            m._x[i,j].start = 1
    
    labeling = { i : j for j in range(k) for i in initial_plan[j] }
    for u,v,j in triples:
        y[u,v,j].start = 0
        if labeling[u] == j and labeling[v] != j:
            y[u,v,j].start = 1
    
    # solve MIP 
    m.Params.IntFeasTol = 1e-7
    m.Params.FeasibilityTol = 1e-7
    m.optimize(m._callback)
    
    if m.solCount == 0:
        print("Did not find a solution. Returning initial_plan.")
        return (initial_plan, float('inf'), m.runtime) 
                            
    # double check that plan is feasible
    plan = [ [ i for i in AG.predecessors(j) if m._x[i,j].x > 0.5 ] for j in range(k) ]
    return (plan, m.objVal, m.runtime)


def AG_contiguity_callback(m, where):

    # check if LP relaxation at this branch-and-bound node has an integer solution
    if where != GRB.Callback.MIPSOL: 
        return
        
    # retrieve the LP solution
    xval = m.cbGetSolution(m._x)
    DG = m._DG
    AG = m._AG

    # check if each district is connected
    for j in range(m._k):

        # which nodes are selected?
        S = list( { i for i in AG.predecessors(j) if xval[i,j] > 0.5 } )
        if len(S)==1 or nx.is_strongly_connected( DG.subgraph(S) ):
            continue

        # for smaller connected components, find some vertex 'a' and add cut.
        max_comp = None
        
        for comp in sorted(nx.strongly_connected_components( DG.subgraph(S) ), key=len, reverse=True):
            if max_comp == None:
                max_comp = comp
                b = nx.utils.arbitrary_element( comp )
                continue

            # pick some vertex from this component
            a = nx.utils.arbitrary_element( comp )

            # add a,b-separator inequality
            C = find_fischetti_separator(DG, comp, b)
            m.cbLazy( m._x[a,j] + m._x[b,j] <= 1 + gp.quicksum( m._x[c,j] for c in C if c in AG.predecessors(j)) )
