import gurobipy as gp
from gurobipy import GRB
import geopandas as gpd
import networkx as nx
import math

# The objective is to minimize the inverse PP score
def build_single_district_mip(DG, minority=None, contiguity=None, root=None, forbidden_fips=None, permitted_fips=None, max_split_counties=None):
    
    assert minority in { None, 'Asian', 'Black', 'Hispanic', 'Native' }
    assert contiguity in { None, 'shir', 'cut' }
    assert forbidden_fips is None or permitted_fips is None
    
    if contiguity == 'shir':
        assert root is not None
        assert root in DG.nodes

    ##################################
    # CREATE MODEL AND MAIN VARIABLES
    ##################################
    
    m = gp.Model()
    
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

    # add constraints saying that edge {u,v} is cut if u (but not v) is selected in the district
    m.addConstrs( m._x[u] - m._x[v] <= m._y[u,v] for u,v in DG.edges )
    
    ###########################
    # ADD OBJECTIVE 
    ###########################
    
    # z is inverse Polsby-Popper score for the district
    m._z = m.addVar(name='z')

    # objective is to minimize the inverse Polsby-Popper score
    m.setObjective( m._z, GRB.MINIMIZE )
    
    ###################################
    # ADD POLSBY-POPPER CONSTRAINTS 
    ###################################
    
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

    m.update()
    
    ###################################
    # ADD MINORITY CONSTRAINT
    ###################################
    
    if minority is not None:
        
        # Census codes available at http://starr.tamu.edu/files/2013/01/Census-Codes.pdf
        if minority == 'Asian':
            codes = ['P0030006'] # Asian alone
            codes += ['P0030013','P0030017','P0030020','P0030023','P0030024'] # Asian (among 2 races)
            codes += ['P0030028','P0030031','P0030034','P0030035','P0030037','P0030040','P0030041','P0030043','P0030044','P0030046'] # 3
            codes += ['P0030048','P0030051','P0030052','P0030054','P0030055','P0030057','P0030058','P0030059','P0030061','P0030062'] # 4
            codes += ['P0030064','P0030065','P0030067','P0030068','P0030069'] # 5
            codes += ['P0030071'] # 6
        elif minority == 'Black':
            codes = ['P0030004'] # Black or African American alone 
            codes += ['P0030011','P0030016','P0030017','P0030018','P0030019'] # Black or African American (among 2 races)
            codes += ['P0030027','P0030028','P0030029','P0030030','P0030037','P0030038','P0030039','P0030040','P0030041','P0030042'] # 3
            codes += ['P0030048','P0030049','P0030050','P0030051','P0030052','P0030053','P0030058','P0030059','P0030060','P0030061'] # 4
            codes += ['P0030064','P0030065', 'P0030066','P0030067','P0030069'] # 5
            codes += ['P0030071'] # 6
        elif minority == 'Hispanic':
            codes = ['P0040002'] # Hispanic or Latino VAP
        elif minority == 'Native':
            codes = ['P0030005'] # American Indian and Alaska Native alone
            codes += ['P0030012','P0030016','P0030020','P0030021','P0030022'] # American Indian and Alaska Native (among 2 races)
            codes += ['P0030027','P0030031','P0030032','P0030033','P0030037','P0030038','P0030039','P0030043','P0030044','P0030045'] # 3
            codes += ['P0030048','P0030049','P0030050','P0030054','P0030055','P0030056','P0030058','P0030059','P0030060','P0030062'] # 4
            codes += ['P0030064','P0030065','P0030066','P0030068','P0030069'] # 5
            codes += ['P0030071'] # 6
        
        for i in DG.nodes:
        
            # voting age population (VAP)
            DG.nodes[i]['VAP'] = DG.nodes[i]['P0030001']
        
            # minority voting age population (MVAP)
            DG.nodes[i]['MVAP'] = sum( DG.nodes[i][code] for code in codes )
            
        # Impose mvap >= 0.5 * vap
        mvap = m.addVar(name='mvap')
        vap = m.addVar(name='vap')

        m.addConstr( mvap == gp.quicksum( DG.nodes[i]['MVAP'] * m._x[i] for i in DG.nodes ) )
        m.addConstr( vap == gp.quicksum( DG.nodes[i]['VAP'] * m._x[i] for i in DG.nodes ) )

        m.addConstr( mvap >= 0.5 * vap )
    
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
        
    elif contiguity == 'cut':
        
        m.Params.LazyConstraints = 1
        m._DG = DG
        m._root = root
        m._callback = cut_callback
        
    ###################
    # FORBIDDEN FIPS
    ###################
    if forbidden_fips is not None:
        for i in DG.nodes:
            if DG.nodes[i]['GEOID20'][0:5] in forbidden_fips:
                m._x[i].UB = 0
                
    ###################
    # PERMITTED FIPS
    ###################
    if permitted_fips is not None:
        for i in DG.nodes:
            m._x[i].UB = 0
            if DG.nodes[i]['GEOID20'][0:5] in permitted_fips:
                m._x[i].UB = 1
                
    ######################
    # MAX SPLIT COUNTIES
    ######################
    if max_split_counties is not None:
        fips = list( { DG.nodes[i]['GEOID20'][0:5] for i in DG.nodes } )
        entire = m.addVars(fips, vtype=GRB.BINARY) # is all of county selected?
        some = m.addVars(fips, vtype=GRB.BINARY) # is at least some of county selected?
        split = m.addVars(fips, vtype=GRB.BINARY) # is county split?

        m.addConstr( gp.quicksum( split ) <= max_split_counties )
        m.addConstrs( split[c] == some[c] - entire[c] for c in fips )

        for i in DG.nodes:
            c = DG.nodes[i]['GEOID20'][0:5]
            m.addConstr( entire[c] <= m._x[i] )
            m.addConstr( m._x[i] <= some[c] )
    
    #######################
    # SAFE VARIABLE FIXING
    #######################
    m.update()
    for u,v in DG.edges:
        if m._x[u].UB < 0.5 or m._x[v].LB > 0.5:
            m._y[u,v].UB = 0
        if m._x[u].LB > 0.5 and m._x[v].UB < 0.5:
            m._y[u,v].LB = 1
    
    ###################################
    # SOLVE PARAMETERS
    ###################################
    m.Params.MIPGap = 0.00
    m.Params.FeasibilityTol = 1e-7
    m.Params.IntFeasTol = 1e-7
    m.update()
    
    return m


from mip_contiguity import find_fischetti_separator 

def cut_callback(m, where):
    if where == GRB.Callback.MIPSOL:
        m._numCallbacks += 1 
        DG = m._DG
        xval = m.cbGetSolution(m._x)

        # vertices assigned to this district (label j)
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
                continue
        
            if b in component: 
                continue

            # find some vertex "a" that has largest population in this component
            a = mpv

            # get minimal a,b-separator
            C = find_fischetti_separator(DG, component, b)

            # add lazy cut
            m.cbLazy( m._x[a] + m._x[b] <= 1 + gp.quicksum( m._x[c] for c in C ) )
            m._numLazyCuts += 1
                
    return


def draw_single_district( filepath, filename, G, district, zoom=False ):
    
    df = gpd.read_file( filepath + filename )
    node_with_this_geoid = { G.nodes[i]['GEOID20'] : i for i in G.nodes }
    assignment = [ -1 for i in G.nodes ]

    if zoom:
        picked = { i : None for i in G.nodes }
    else:
        picked = { i : False for i in G.nodes }
    
    for i in district:
        picked[i] = True

    for u in range(G.number_of_nodes()):
        geoid = df['GEOID20'][u]
        i = node_with_this_geoid[geoid]
        assignment[u] = picked[i]

    df['assignment'] = assignment
    my_fig = df.plot(column='assignment').get_figure()
    return 
