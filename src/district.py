import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import math
from census_codes import get_census_codes
from mip_contiguity import find_fischetti_separator


def build_single_district_mip(DG, objective='polsby_popper', contiguity=None, root=None, verbose=True,
                              split_counties_limit=None, deviation_penalty=0.0,
                              minority=None, mvap_lower=0.5, mvap_upper=1.0, mvap_excess_penalty=0.0, 
                              complement_contiguity=None, complement_balance=False, complement_compactness_penalty=0.0):
    
    # sanity check the inputs
    assert objective in {'polsby_popper', 'cut_edges'}
    assert minority in { None, 'Asian', 'Black', 'Hispanic', 'Native' }
    assert contiguity in { None, 'shir', 'cut' }
    assert complement_contiguity in { None, 'cut' } # shir not supported
    assert objective=='polsby_popper' or complement_compactness_penalty==0.0, "Can use complement_compactness_penalty only for PP objective"
    if contiguity == 'shir':
        assert root is not None
        assert root in DG.nodes
        
    ##################################
    # CREATE MODEL AND MAIN VARIABLES
    ##################################
    
    m = gp.Model()
    if not verbose:
        m.Params.OutputFlag = 0
    
    # x[i] equals one when node i is selected in the district
    m._x = m.addVars(DG.nodes, name='x', vtype=GRB.BINARY)

    # y[u,v] equals one when arc (u,v) is cut because u (but not v) is selected in the district
    m._y = m.addVars(DG.edges, name='y', vtype=GRB.BINARY)

    ###########################
    # ADD MAIN CONSTRAINTS
    ###########################
    
    # add constraints saying that the district has population at least L and at most U
    m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i] for i in DG.nodes) >= DG._L )
    m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i] for i in DG.nodes) <= DG._U )
    
    # add constraints saying that the complement has population at least CL and at most CU
    if complement_balance:
        m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * (1-m._x[i]) for i in DG.nodes) >= DG._CL )
        m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * (1-m._x[i]) for i in DG.nodes) <= DG._CU )

    # add constraints saying that edge {u,v} is cut if u (but not v) is selected in the district
    m.addConstrs( m._x[u] - m._x[v] <= m._y[u,v] for u,v in DG.edges )
    
    ###########################
    # ADD OBJECTIVE 
    ###########################
    
    if objective == 'polsby_popper':
        # z is inverse Polsby-Popper score for the district
        m._z = m.addVar(name='z')

        # objective is to minimize the inverse Polsby-Popper score
        m.setObjective( m._z, GRB.MINIMIZE )

        if complement_compactness_penalty > 0:
            m._zc = m.addVar(name='zc')
            m._zc.obj = complement_compactness_penalty
        
    elif objective == 'cut_edges':
        undirected_edges = [ (i,j) for i,j in DG.edges if i < j ]
        m.addConstrs( m._y[i,j] == m._y[j,i] for i,j in undirected_edges )
        m.setObjective( gp.quicksum( m._y[i,j] for i,j in undirected_edges ), GRB.MINIMIZE )
    
    ###################################
    # ADD POLSBY-POPPER CONSTRAINTS 
    ###################################
    
    if objective == 'polsby_popper':
        
        # A = area of the district
        m._A = m.addVar(name='A')

        # P = perimeter of the district
        m._P = m.addVar(name='P')

        # add SOCP constraint relating inverse Polsby-Popper score z to area and perimeter
        m.addConstr( m._P * m._P <= 4 * math.pi * m._A * m._z )

        # add constraint on area A
        m.addConstr( m._A == gp.quicksum( DG.nodes[i]['area'] * m._x[i] for i in DG.nodes ) )

        # add constraint on perimeter P
        m.addConstr( m._P == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v] for u,v in DG.edges )
                     + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )

        if complement_compactness_penalty > 0:
            m._Ac = m.addVar(name='Ac')
            m._Pc = m.addVar(name='Pc')
            m.addConstr( m._Pc * m._Pc <= 4 * math.pi * m._Ac * m._zc )
            m.addConstr( m._Ac == gp.quicksum( DG.nodes[i]['area'] * (1-m._x[i]) for i in DG.nodes ) )
            m.addConstr( m._Pc == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v] for u,v in DG.edges )
                         + gp.quicksum( DG.nodes[i]['boundary_perim'] * (1-m._x[i]) for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )
        
        m.update()
    
    ###################################
    # ADD MINORITY CONSTRAINT(S)
    ###################################
    
    if minority is not None:
        
        codes = get_census_codes(minority)
        
        for i in DG.nodes:
        
            # voting age population (VAP)
            DG.nodes[i]['VAP'] = DG.nodes[i]['P0030001']
        
            # minority voting age population (MVAP)
            DG.nodes[i]['MVAP'] = sum( DG.nodes[i][code] for code in codes )
            
        # Idea: impose mvap >= 0.5 * vap
        m._mvap = m.addVar(name='mvap')
        m._vap = m.addVar(name='vap')

        m.addConstr( m._mvap == gp.quicksum( DG.nodes[i]['MVAP'] * m._x[i] for i in DG.nodes ) )
        m.addConstr( m._vap == gp.quicksum( DG.nodes[i]['VAP'] * m._x[i] for i in DG.nodes ) )

        excess = m.addVar(name='excess')
        excess.obj = mvap_excess_penalty
        m.addConstr( m._mvap - excess == mvap_lower * m._vap )
        
        # don't exceed, say, 80% BVAP
        if mvap_upper < 1:
            m.addConstr( m._mvap <= mvap_upper * m._vap )
        
    ###################################
    # ADD CONTIGUITY CONSTRAINTS
    ###################################
    
    m._callback = None
    m._numCallbacks = 0
    m._numLazyCuts = 0
    
    if contiguity == 'shir':
        
        m._x[root].LB = 1
        M = DG.number_of_nodes() - 1
        
        # Add flow variables: f[u,v] = amount of flow sent across arc uv 
        m._f = m.addVars( DG.edges, name='f' )
        
        # if selected but not a root, consume one unit of flow
        m.addConstrs( gp.quicksum( m._f[j,i] - m._f[i,j] for j in DG.neighbors(i) ) == m._x[i] for i in DG.nodes if i != root )

        # flow can only enter selected nodes
        m.addConstrs( gp.quicksum( m._f[j,i] for j in DG.neighbors(i) ) <= M * m._x[i] for i in DG.nodes if i != root )
        
    if contiguity == 'cut' or complement_contiguity == 'cut':
        
        m.Params.LazyConstraints = 1
        m._DG = DG
        m._root = root
        m._contiguity = contiguity
        m._complement_contiguity = complement_contiguity
        m._callback = cut_callback
        
    ##########################
    # LIMIT SPLIT COUNTIES?
    ##########################
    
    if split_counties_limit is not None:
        
        fips = list( { DG.nodes[i]['GEOID20'][0:5] for i in DG.nodes } )
        is_all = m.addVars(fips, vtype=GRB.BINARY) # is all of county selected?
        is_some = m.addVars(fips, vtype=GRB.BINARY) # is (at least) some of county selected?
        is_split = m.addVars(fips, vtype=GRB.BINARY) # is county split?

        m.addConstr( gp.quicksum( is_split ) <= split_counties_limit )
        m.addConstrs( is_split[c] == is_some[c] - is_all[c] for c in fips )

        for i in DG.nodes:
            c = DG.nodes[i]['GEOID20'][0:5]
            m.addConstr( is_all[c] <= m._x[i] )
            m.addConstr( m._x[i] <= is_some[c] )
            
    ###################################
    # IDEAL PENALTY
    ###################################
    
    if deviation_penalty > 0.0:
        
        below = m.addVar()
        below.obj = deviation_penalty
        
        above = m.addVar()
        above.obj = deviation_penalty
        
        total = sum( DG.nodes[i]['TOTPOP'] for i in DG.nodes )
        ideal = total * ( DG._L + DG._U ) / ( DG._L + DG._U + DG._CL + DG._CU )
        m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i] for i in DG.nodes) == ideal + above - below )
        
    ###################################
    # SOLVE PARAMETERS
    ###################################
    
    m.Params.MIPGap = 0.00
    m.Params.FeasibilityTol = 1e-7
    m.Params.IntFeasTol = 1e-7
    m.update()
    
    return m


def cut_callback(m, where):
    
    if where == GRB.Callback.MIPSOL:
        m._numCallbacks += 1 
        DG = m._DG
        xval = m.cbGetSolution(m._x)

        ##########################################
        # ADD CUT FOR COMPLEMENT?
        ##########################################
        
        if m._contiguity == 'cut':
            # vertices assigned to this district 
            S = [ v for v in DG.nodes if xval[v] > 0.5 ]

            # what shall we deem as the "root" of this district? call it b
            b = m._root # possibly None
            
            # for each component that doesn't contain b, add a cut
            for component in sorted( nx.strongly_connected_components( DG.subgraph(S) ), key=len, reverse=True ):

                # what is the maximum population node in this component?
                maxp = max( DG.nodes[v]['TOTPOP'] for v in component)
                mpv = [ v for v in component if DG.nodes[v]['TOTPOP'] == maxp ][0]

                # if no root 'b' has been selected yet, pick one
                if b is None:
                    # find some vertex "b" that has largest population in this component
                    b = mpv

                if b in component: 
                    continue

                # find some vertex "a" that has largest population in this component
                a = mpv

                # get minimal a,b-separator
                C = find_fischetti_separator(DG, component, b)

                # add lazy cut
                m.cbLazy( m._x[a] + m._x[b] <= 1 + gp.quicksum( m._x[c] for c in C ) )
                m._numLazyCuts += 1
            
        ##########################################
        # ADD CUT FOR COMPLEMENT?
        ##########################################
        
        if m._complement_contiguity == 'cut':
            
            # vertices assigned to the complement
            S = [ v for v in DG.nodes if xval[v] < 0.5 ]

            # what shall we deem as the "root" of the complement? call it b
            b = None
            
            # for each component that doesn't contain b, add a cut
            for component in sorted( nx.strongly_connected_components( DG.subgraph(S) ), key=len, reverse=True ):

                # what is the maximum population node in this component?
                maxp = max( DG.nodes[v]['TOTPOP'] for v in component )
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
                # replace x by 1-x in:
                #      m.cbLazy( m._x[a] + m._x[b] <= 1 + gp.quicksum( m._x[c] for c in C ) )
                # i.e., if neither a nor b is picked in district, then not all of C can be picked in district
                m.cbLazy( gp.quicksum( m._x[c] for c in C ) + 1 <= m._x[a] + m._x[b] + len(C) )
                m._numLazyCuts += 1
                
    return
