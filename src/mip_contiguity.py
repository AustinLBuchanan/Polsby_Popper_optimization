import xprgrb as gp
from xprgrb import GRB
import networkx as nx

def lcut_callback(m, where):
    if where == GRB.Callback.MIPSOL:
        m._numCallbacks += 1 
        DG = m._DG
        xval = m.cbGetSolution(m._x)
        rval = m.cbGetSolution(m._r)

        for j in range(DG._k):
            
            # vertices assigned to this district (label j)
            S = [ v for v in DG.nodes if xval[v,j] > 0.5 ]
            
            # what shall we deem as the "root" of this district? call it b
            roots = [ i for i in DG.nodes if rval[i,j] > 0.5 ]
            b = roots[0]
            
            # for each component that doesn't contain b, add a cut
            for component in nx.strongly_connected_components(DG.subgraph(S)):
                
                if b in component: 
                    continue
                
                # find some vertex "a" that has largest population in this component
                maxp = max( DG.nodes[v]['TOTPOP'] for v in component)
                maxp_vertices = [ v for v in component if DG.nodes[v]['TOTPOP'] == maxp ]
                a = maxp_vertices[0]  
                    
                # get minimal a,b-separator
                C = find_fischetti_separator(DG, component, b)
                
                # make it a minimal *length-U* a,b-separator
                for (u,v) in DG.edges():
                    DG[u][v]['lcutweight'] = DG.nodes[u]['TOTPOP']  
                    
                # "remove" C from graph
                for c in C:
                    for node in DG.neighbors(c):
                        DG[c][node]['lcutweight'] = DG._U+1
                
                # is C\{c} a length-U a,b-separator still? If so, remove c from C
                drop_from_C = []
                for c in C:
                    
                    # temporarily add c back to graph (i.e., "remove" c from cut C)
                    for node in DG.neighbors(c):
                        DG[c][node]['lcutweight'] = DG.nodes[c]['TOTPOP']
                    
                    # what is distance from a to b in G-C ?
                    distance_from_a = nx.single_source_dijkstra_path_length(DG, a, weight='lcutweight')
                    
                    if distance_from_a[b] + DG.nodes[b]['TOTPOP'] > DG._U:
                        # c was not needed in the cut C. delete c from C
                        drop_from_C.append(c)
                    else:
                        # keep c in C. revert arc weights back to "infinity"
                        for node in DG.neighbors(c):
                            DG[c][node]['lcutweight'] = DG._U + 1
                    
                # add lazy cut
                minC = [ c for c in C if c not in drop_from_C ]
                m.cbLazy( m._x[a,j] + m._x[b,j] <= 1 + gp.quicksum( m._x[c,j] for c in minC ) )
                m._numLazyCuts += 1
                
    return


def lcut_callback_xpress_preintsol(prob, m, soltype, cutoff):

    DG = m._DG

    try:
        x = []
        prob.getlpsol(x)
    except:
        return (1, 0)

    xval = {(i,j): x[m.xmodel.getIndex(m._x[i,j])] for i in DG.nodes for j in range(DG._k)}
    rval = {(i,j): x[m.xmodel.getIndex(m._r[i,j])] for i in DG.nodes for j in range(DG._k)}

    for j in range(DG._k):

        # vertices assigned to this district (label j)
        S = [ v for v in DG.nodes if xval[v,j] > 0.5 ]

        # what shall we deem as the "root" of this district? call it b
        roots = [ i for i in DG.nodes if rval[i,j] > 0.5 ]

        # if b does not exist, skip to the next district
        if len(roots) < 1:
            continue

        b = roots[0]

        # for each component that doesn't contain b, add a cut
        for component in nx.strongly_connected_components(DG.subgraph(S)):

            if b in component:
                continue

            # find some vertex "a" that has largest population in this component
            maxp = max( DG.nodes[v]['TOTPOP'] for v in component)
            maxp_vertices = [ v for v in component if DG.nodes[v]['TOTPOP'] == maxp ]
            a = maxp_vertices[0]

            # get minimal a,b-separator
            C = find_fischetti_separator(DG, component, b)

            # make it a minimal *length-U* a,b-separator
            for (u,v) in DG.edges():
                DG[u][v]['lcutweight'] = DG.nodes[u]['TOTPOP']

            # "remove" C from graph
            for c in C:
                for node in DG.neighbors(c):
                    DG[c][node]['lcutweight'] = DG._U+1

            # is C\{c} a length-U a,b-separator still? If so, remove c from C
            drop_from_C = []
            for c in C:

                # temporarily add c back to graph (i.e., "remove" c from cut C)
                for node in DG.neighbors(c):
                    DG[c][node]['lcutweight'] = DG.nodes[c]['TOTPOP']

                # what is distance from a to b in G-C ?
                distance_from_a = nx.single_source_dijkstra_path_length(DG, a, weight='lcutweight')

                if distance_from_a[b] + DG.nodes[b]['TOTPOP'] > DG._U:
                    # c was not needed in the cut C. delete c from C
                    drop_from_C.append(c)
                else:
                    # keep c in C. revert arc weights back to "infinity"
                    for node in DG.neighbors(c):
                        DG[c][node]['lcutweight'] = DG._U + 1

            minC = [ c for c in C if c not in drop_from_C ]

            if xval[a,j] + xval[b,j] - sum(xval[c,j] for c in minC) > 1 + prob.controls.feastol:
                return (1, cutoff)

    return (0, cutoff)


def lcut_callback_xpress_optnode(prob, m, sol, lb, ub):

    DG = m._DG

    xval = {(i,j): sol[m.xmodel.getIndex(m._x[i,j])] for i in DG.nodes for j in range(DG._k)}
    rval = {(i,j): sol[m.xmodel.getIndex(m._r[i,j])] for i in DG.nodes for j in range(DG._k)}

    for j in range(DG._k):

        # vertices assigned to this district (label j)
        S = [ v for v in DG.nodes if xval[v,j] > 0.5 ]

        # what shall we deem as the "root" of this district? call it b
        roots = [ i for i in DG.nodes if rval[i,j] > 0.5 ]

        # if b does not exist, skip to the next district
        if len(roots) < 1:
            continue

        b = roots[0]

        # for each component that doesn't contain b, add a cut
        for component in nx.strongly_connected_components(DG.subgraph(S)):

            if b in component:
                continue

            # find some vertex "a" that has largest population in this component
            maxp = max( DG.nodes[v]['TOTPOP'] for v in component)
            maxp_vertices = [ v for v in component if DG.nodes[v]['TOTPOP'] == maxp ]
            a = maxp_vertices[0]

            # get minimal a,b-separator
            C = find_fischetti_separator(DG, component, b)

            # make it a minimal *length-U* a,b-separator
            for (u,v) in DG.edges():
                DG[u][v]['lcutweight'] = DG.nodes[u]['TOTPOP']

            # "remove" C from graph
            for c in C:
                for node in DG.neighbors(c):
                    DG[c][node]['lcutweight'] = DG._U+1

            # is C\{c} a length-U a,b-separator still? If so, remove c from C
            drop_from_C = []
            for c in C:

                # temporarily add c back to graph (i.e., "remove" c from cut C)
                for node in DG.neighbors(c):
                    DG[c][node]['lcutweight'] = DG.nodes[c]['TOTPOP']

                # what is distance from a to b in G-C ?
                distance_from_a = nx.single_source_dijkstra_path_length(DG, a, weight='lcutweight')

                if distance_from_a[b] + DG.nodes[b]['TOTPOP'] > DG._U:
                    # c was not needed in the cut C. delete c from C
                    drop_from_C.append(c)
                else:
                    # keep c in C. revert arc weights back to "infinity"
                    for node in DG.neighbors(c):
                        DG[c][node]['lcutweight'] = DG._U + 1

            minC = [ c for c in C if c not in drop_from_C ]

            if xval[a,j] + xval[b,j] - sum(xval[c,j] for c in minC) > 1 + prob.controls.feastol:

                # add cut

                indices = [m._x[a,j], m._x[b,j]] + [m._x[c,j] for c in minC]
                coeffs  = [1,         1]         + [-1] * len(minC)
                rhs = 1

                mcolsp, dvalp = [], []
                drhsp, status = prob.presolverow('L', indices, coeffs, rhs, prob.attributes.cols,
                                                 mcolsp, dvalp)

                if status >= 0:
                    prob.addcuts([1], ['L'], [drhsp], [0, len(mcolsp)], mcolsp, dvalp)

    return 0


def find_fischetti_separator(DG, component, b):
    neighbors_component = [ False for i in DG.nodes ]
    for i in nx.node_boundary(DG, component, None):
        neighbors_component[i] = True
    
    visited = [ False for i in DG.nodes ]
    child = [ b ]
    visited[b] = True
    
    while child:
        parent = child
        child = []
        for i in parent:
            if not neighbors_component[i]:
                for j in DG.neighbors(i):
                    if not visited[j]:
                        child.append(j)
                        visited[j] = True
    
    C = [ i for i in DG.nodes if neighbors_component[i] and visited[i] ]
    return C


def add_shir_constraints(m, DG):

    from xprgrb import solver

    # g[i,j] = amount of flow generated at node i of type j
    g = m.addVars(DG.nodes, DG._k, name='g')
    
    # f[j,u,v] = amount of flow sent across arc uv of type j
    if solver == 'gurobi':
        f = m.addVars( DG._k, DG.edges, name='f' )
    else:
        f = m.addVars([(j,u,v) for j in range(DG._k) for (u,v) in DG.edges], name='f')

    # compute big-M  
    M = most_possible_nodes_in_one_district(DG) - 1
    
    # flow can only be generated at roots
    m.addConstrs( g[i,j] <= (M+1)*m._r[i,j] for i in DG.nodes for j in range(DG._k) )
    
    # flow balance
    m.addConstrs( g[i,j] - m._x[i,j] == gp.quicksum( f[j,i,u]-f[j,u,i] for u in DG.neighbors(i)) for i in DG.nodes for j in range(DG._k) )
    
    # flow type j can enter vertex i only if (i is assigned to district j) and (i is not root of j)
    m.addConstrs( gp.quicksum( f[j,u,i] for u in DG.neighbors(i) ) <= M * (m._x[i,j]-m._r[i,j]) for i in DG.nodes for j in range(DG._k) )
    return


def add_scf_constraints(m, DG):
    
    # for big M
    most = most_possible_nodes_in_one_district(DG) 
    
    # Add flow variables: f[u,v] = amount of flow sent across arc uv 
    #  Flows are sent across arcs of DG
    m._f = m.addVars( DG.edges, name='f' )

    # if not a root, consume some flow.
    # if a root, only send out (so much) flow.
    m.addConstrs( gp.quicksum( m._f[v,u] - m._f[u,v] for v in DG.neighbors(u) ) 
                 >= 1 - ( most ) * gp.quicksum( m._r[u,j] for j in range(DG._k) ) for u in DG.nodes )

    # do not send flow across cut edges
    m.addConstrs( m._f[u,v] + m._f[v,u] <= ( most - 1 ) * ( 1 - gp.quicksum( m._y[u,v,j] for j in range(DG._k) ) ) for u,v in DG.edges )

    m.update()
    return

def most_possible_nodes_in_one_district(DG):
    cumulative_population = 0
    num_nodes = 0
    population = [ DG.nodes[i]['TOTPOP'] for i in DG.nodes ]
    for ipopulation in sorted(population):
        cumulative_population += ipopulation
        num_nodes += 1
        if cumulative_population > DG._U:
            return num_nodes - 1
    return len(population)


def connectivity_preprocess(m, DG):
    """
    Based on the connectivity of the underlying undirected graph, find
    articulation node (a.k.a. vertex cuts of cardinality 1) and, for
    any such node i, check if each connected component C created by
    removing i has a total population below the district lower
    bound. If so, C must be in the same district as i and this can be
    enforced by equalling all x[h,j] with x[i,j] for all j's and for
    all h in C.
    """

    UG = DG.to_undirected()

    if nx.is_biconnected(UG):
        return

    #import matplotlib.pyplot as plt
    #nx.draw_networkx(UG)
    #plt.show()

    VC = [i for i in nx.articulation_points(UG)]

    locked_nodes = set()

    for i in VC:
        fwstar = [j for j in UG[i]]
        UG.remove_node(i)

        #import matplotlib.pyplot as plt
        #nx.draw_networkx(UG)
        #plt.show()

        components = nx.connected_components(UG)

        for c in components:
            if len(c) >= 1 and sum(DG.nodes[h]['TOTPOP'] for h in c) < DG._L:
                locked_nodes = locked_nodes.union([h for h in c])
                m.addConstrs(m._x[i,j] == m._x[h,j] for j in range(DG._k) for h in c)
                #m.addConstrs(m._y[i,h,j] == 0       for j in range(DG._k) for h in c if (i,h) in DG.edges)

        UG.add_node(i)
        UG.add_edges_from([(i,j) for j in fwstar])

    if len(locked_nodes) > 0:
        print(f"Locked assignment for {len(locked_nodes)} nodes: {locked_nodes}")
    else:
        print("No locked assignment")
