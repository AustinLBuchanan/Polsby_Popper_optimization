import networkx as nx
import xprgrb as gp
from xprgrb import GRB, solver
import xpress as xp
from pyproj import Proj
import random
from mip_contiguity import most_possible_nodes_in_one_district

def sq_eucl_dist(x1,y1,x2,y2):
    return (x1-x2)**2 + (y1-y2)**2

def solve_hess_model(DG):
    
    print("Solving Hess MIP to get a heuristic warm start...")
    
    # map projection
    epsg = get_epsg(DG._state)
    p = Proj("EPSG:"+epsg, preserve_units=True) 
    
    # get lat and long
    for i in DG.nodes:
        DG.nodes[i]['C_X'] = DG.nodes[i]['INTPTLON20']  # longitude of county's center
        DG.nodes[i]['C_Y'] = DG.nodes[i]['INTPTLAT20']  # latitude of county's center
        
        # projection to get (x,y) coordinates (units of KM)
        (x,y) = p( DG.nodes[i]['C_X'], DG.nodes[i]['C_Y'] )
        DG.nodes[i]['X'] = x / 1000
        DG.nodes[i]['Y'] = y / 1000  
        
    # create model 
    m = gp.Model()
    m.Params.OutputFlag = 0

    # create x[i,j] variable which equals one when county i 
    #    is assigned to (the district centered at) county j
    x = m.addVars( DG.nodes, DG.nodes, name='x', vtype=GRB.BINARY )
    
    # objective is to minimize the moment of inertia: sum (d^2 * p * x over all i and j)
    m.setObjective( gp.quicksum( round( sq_eucl_dist(DG.nodes[i]['X'],DG.nodes[i]['Y'],DG.nodes[j]['X'],DG.nodes[j]['Y']) * DG.nodes[i]['TOTPOP'] / 1000 ) * x[i,j] for i in DG.nodes for j in DG.nodes ), GRB.MINIMIZE )
    
    # add constraints saying that each county i is assigned to one district
    m.addConstrs( gp.quicksum( x[i,j] for j in DG.nodes ) == 1 for i in DG.nodes )

    # add constraint saying there should be k district centers
    m.addConstr( gp.quicksum( x[j,j] for j in DG.nodes ) == DG._k )

    # add constraints that say: if j roots a district, then its population is between L and U.
    m.addConstrs( gp.quicksum( DG.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) >= DG._L * x[j,j] for j in DG.nodes )
    m.addConstrs( gp.quicksum( DG.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) <= DG._U * x[j,j] for j in DG.nodes )

    # add coupling constraints saying that if i is assigned to j, then j is a center.
    m.addConstrs( x[i,j] <= x[j,j] for i in DG.nodes for j in DG.nodes )
    
    # add flow variables
    #    f[i,j,v] = flow across arc (i,j) that is sent from souce/root v
    if solver == 'gurobi':
        f = m.addVars( DG.edges, DG.nodes, name='f' )
    else:
        f = m.addVars([(u,v,j) for (u,v) in DG.edges for j in DG.nodes], name='f')

    # add constraints saying that if node i is assigned to node j, 
    #   then node i must consume one unit of node j's flow
    m.addConstrs( gp.quicksum( f[u,i,j] - f[i,u,j] for u in DG.neighbors(i) ) == x[i,j] for i in DG.nodes for j in DG.nodes if i != j )

    # add constraints saying that node i can receive flow of type j 
    #   only if node i is assigned to node j
    M = most_possible_nodes_in_one_district(DG) - 1
    m.addConstrs( gp.quicksum( f[u,i,j] for u in DG.neighbors(i) ) <= M * x[i,j] for i in DG.nodes for j in DG.nodes if i != j )

    # add constraints saying that node j cannot receive flow of its own type
    m.addConstrs( gp.quicksum( f[u,j,j] for u in DG.neighbors(j) ) == 0 for j in DG.nodes )

    # solve
    m.Params.TimeLimit = 60
    m.Params.IntFeasTol = 1.e-9
    m.Params.FeasibilityTol = 1.e-9
    
    # apply diagonal fixing first, to test feasibility
    for i in DG.nodes:
        for j in DG.nodes:
            if DG.nodes[i]['TOTPOP'] > DG.nodes[j]['TOTPOP']:
                gp.setUB(x[i,j], m, 0)
                
    m.optimize()
    grb_time = m.runtime
    
    # if infeasible, nothing to report
    if m.solCount <= 0:
        return (False, grb_time)
    
    # undo diagonal fixing and re-solve for compactness
    for i in DG.nodes:
        for j in DG.nodes:
            gp.setUB(x[i,j], m, 1)
        
    m.optimize()
    grb_time += m.runtime
    
    # return labeling
    roots = [ j for j in DG.nodes if gp.getsol(x[j,j], m) > 0.5 ]
    
    # maps a center *vertex* to its district *number*
    root_map = { roots[pos] : pos for pos in range(DG._k) }
    
    labeling = { i : root_map[j] for i in DG.nodes for j in DG.nodes if gp.getsol(x[i,j],m) > 0.5 }
    return (labeling, grb_time)

#     # need to sort w.r.t. ordering, because of partitioning orbitope
    
#     # for each district, find its earliest vertex in the ordering
#     # then, sort these earliest vertices to get the district labels
#     unmapped_labeling = { i : j for i in DG.nodes for j in DG.nodes if x[i,j].x > 0.5 }
#     district_map = { j : -1 for j in DG.nodes if x[j,j].x > 0.5 }
#     labeling = { i : -1 for i in DG.nodes }
#     count = 0
    
#     for i in DG._ordering:
#         j = unmapped_labeling[i]
#         if district_map[j] == -1:
#             district_map[j] = count
#             count += 1
        
#         labeling[i] = district_map[j]
    
#     return (labeling, grb_time)


def hess_heuristic(DG, impose_contiguity=True):
    
    print("Applying Hess heuristic to get a MIP warm start...")
    
    # map projection
    epsg = get_epsg(DG._state)
    p = Proj("EPSG:"+epsg, preserve_units=True) 
    
    # get lat and long
    for i in DG.nodes:
        DG.nodes[i]['C_X'] = DG.nodes[i]['INTPTLON20']  # longitude of county's center
        DG.nodes[i]['C_Y'] = DG.nodes[i]['INTPTLAT20']  # latitude of county's center
        
        # projection to get (x,y) coordinates (units of KM)
        (x,y) = p( DG.nodes[i]['C_X'], DG.nodes[i]['C_Y'] )
        DG.nodes[i]['X'] = x / 1000
        DG.nodes[i]['Y'] = y / 1000  
    
    # get initial roots
    roots = list()
    nodes = [ i for i in DG.nodes ]
    weights = [ 1 for i in DG.nodes ]
    means_x = list()
    means_y = list()

    while len(roots) < DG._k:

        r = random.choices(nodes, weights=weights, k=1)[0]
        if r in roots:
            continue

        roots.append(r)
        means_x.append( DG.nodes[r]['X'] ) 
        means_y.append( DG.nodes[r]['Y'] )

        dist = [ sq_eucl_dist(DG.nodes[i]['X'],DG.nodes[i]['Y'],DG.nodes[r]['X'],DG.nodes[r]['Y']) for i in DG.nodes ]

        if len(roots) == 1:
            # initialize weights
            weights = [ DG.nodes[i]['TOTPOP'] * dist[i] * dist[i] for i in DG.nodes ]
        else:
            # update weights
            weights = [ min( weights[i], DG.nodes[i]['TOTPOP'] * dist[i] * dist[i] ) for i in DG.nodes ]

    #print("Initial roots =",roots)  
    
    # create model 
    m = gp.Model()
    m.Params.OutputFlag = 0

    # create x[i,j] variable which equals one when county i is 
    #    assigned to district j (actually, use continuous 0-1 now)
    x = m.addVars( DG.nodes, DG._k, name='x' )
    
    # add constraints saying that each county i is assigned to one district
    m.addConstrs( gp.quicksum( x[i,j] for j in range(DG._k) ) == 1 for i in DG.nodes )

    # add constraints that say: if j roots a district, then its population is between L and U.
    m.addConstrs( gp.quicksum( DG.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) >= DG._L for j in range(DG._k) )
    m.addConstrs( gp.quicksum( DG.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) <= DG._U for j in range(DG._k) )

    total_cost = sum( weights )
    max_iterations = 10
    grb_time = 0
    #print("Iteration \t Objective \t iter_time")
    
    for iteration in range(max_iterations):

        old_total_cost = total_cost
        total_cost = 0

        # get clusters
        m.setObjective( gp.quicksum( DG.nodes[i]['TOTPOP'] * sq_eucl_dist(DG.nodes[i]['X'],DG.nodes[i]['Y'],means_x[j],means_y[j]) * x[i,j] for i in DG.nodes for j in range(DG._k) ), GRB.MINIMIZE )
        m.optimize()
        #print(iteration,"\t",'{0:.2f}'.format(m.objVal),"\t",'{0:.2f}'.format(m.runtime))

        total_cost = m.objVal
        grb_time += m.runtime
        clusters = [ [ i for i in DG.nodes if gp.getsol(x[i,j],m) > 0.5 ] for j in range(DG._k) ]

        # convergence?
        if total_cost == old_total_cost:
            break
        else:
            # get next means
            population = [ sum( DG.nodes[i]['TOTPOP'] * gp.getsol(x[i,j],m) for i in DG.nodes ) for j in range(DG._k) ]
            means_x = [ sum( DG.nodes[i]['TOTPOP'] * DG.nodes[i]['X'] * gp.getsol(x[i,j],m) for i in DG.nodes ) / population[j] for j in range(DG._k) ]
            means_y = [ sum( DG.nodes[i]['TOTPOP'] * DG.nodes[i]['Y'] * gp.getsol(x[i,j],m) for i in DG.nodes ) / population[j] for j in range(DG._k) ]

    # now solve as an IP
    global solver

    if solver == 'gurobi':
        for i in DG.nodes:
            for j in range(DG._k):
                x[i,j].vtype = GRB.BINARY
    else:
        m.xmodel.chgcoltype([x[i,j] for i in DG.nodes for j in range(DG._k)], ['B'] * DG._k * len(DG.nodes))

    # get clusters
    m.setObjective( gp.quicksum( DG.nodes[i]['TOTPOP'] * sq_eucl_dist(DG.nodes[i]['X'],DG.nodes[i]['Y'],means_x[j],means_y[j]) * x[i,j] for i in DG.nodes for j in range(DG._k)), GRB.MINIMIZE )
    m.Params.MIPGap = 0.01 # 1% gap permitted
    m.Params.TimeLimit = 60 # 60 second limit
    m.Params.IntFeasTol = 1.e-9
    m.Params.FeasibilityTol = 1.e-9
    #m.Params.OutputFlag = 1
    m.optimize()
    grb_time += m.runtime
    #print("Hess heuristic IP time =",m.runtime)

    # if infeasible, nothing to report
    if m.solCount <= 0:
        print("WARNING: Hess heuristic failed to find IP solution. Aborting...")
        return (False, grb_time)
    
    # check connectivity. If necessary, add contiguity constraints and re-solve
    connected = True
    for j in range(DG._k):
        district = [ i for i in DG.nodes if gp.getsol(x[i,j],m) > 0.5 ]
        is_conn = nx.is_strongly_connected(DG.subgraph(district)) 
        #print("In Hess heuristic, is district",j,"connected?",is_conn)
        if not is_conn:
            connected = False
            
    # add DAG constraints (which are sufficent but not necessary for contiguity);
    #    1. find 'root' of district j that is closest to its district's (x,y)-means
    #    2. find BFS tree of DG rooted at 'root'
    #    3. get associated vertex ordering ('root', ..., )
    #    4. if x[i,j]=1 then sum( x[v,j] for v in N(i) prior in ordering) >= 1.
    if impose_contiguity and not connected:
        for j in range(DG._k):
            district = [ i for i in DG.nodes if gp.getsol(x[i,j],m) > 0.5 ]
            min_dist = min( sq_eucl_dist(DG.nodes[i]['X'],DG.nodes[i]['Y'],means_x[j],means_y[j]) for i in district)
            roots = [ i for i in district if min_dist == sq_eucl_dist(DG.nodes[i]['X'],DG.nodes[i]['Y'],means_x[j],means_y[j]) ]
            root = roots[0]
            
            pos = get_bfs_position(DG, root)
            m.addConstrs( x[i,j] <= gp.quicksum( x[v,j] for v in DG.neighbors(i) if pos[v] < pos[i] ) for i in DG.nodes if i != root )
            
        m.optimize()
        grb_time += m.runtime
    
    # if infeasible, nothing to report
    if m.solCount <= 0:
        print("WARNING: Hess heuristic failed to find *contiguous* IP solution. Aborting...")
        return (False, grb_time)
    else:
        print("Success! Hess heuristic found a contiguous IP solution.")
    
    labeling = { i : j for i in DG.nodes for j in range(DG._k) if gp.getsol(x[i,j],m) > 0.5 }
    return (labeling, grb_time)
  
    
def get_bfs_position(DG, root):
    
    position = { i : -1 for i in DG.nodes }
    position[root] = 0
    children = [ root ]
    p = 1
    
    while children:
        parents = children
        children = list()
        for i in parents:
            for j in DG.neighbors(i):
                if position[j] == -1:
                    children.append(j)
                    position[j] = p
                    p += 1
                    
    return position
    
def get_epsg(state):
    for fips in states.keys():
        abbr = states[fips]['abbr']
        if abbr == state:
            epsg = states[fips]['epsg']
            return epsg
    print("ERROR: Could not find epsg.")
    return

    
states = {
    '01': {'abbr': 'AL', 'epsg': '3465', 'name': 'Alabama'},
    '02': {'abbr': 'AK', 'epsg': '3471', 'name': 'Alaska'},
    '04': {'abbr': 'AZ', 'epsg': '3478', 'name': 'Arizona'},
    '05': {'abbr': 'AR', 'epsg': '3484', 'name': 'Arkansas'},
    '06': {'abbr': 'CA', 'epsg': '3493', 'name': 'California'},
    '08': {'abbr': 'CO', 'epsg': '3501', 'name': 'Colorado'},
    '09': {'abbr': 'CT', 'epsg': '3507', 'name': 'Connecticut'},
    '10': {'abbr': 'DE', 'epsg': '3509', 'name': 'Delaware'},
    '12': {'abbr': 'FL', 'epsg': '3514', 'name': 'Florida'},
    '13': {'abbr': 'GA', 'epsg': '3518', 'name': 'Georgia'},
    '15': {'abbr': 'HI', 'epsg': '2784', 'name': 'Hawaii'},
    '16': {'abbr': 'ID', 'epsg': '3524', 'name': 'Idaho'},
    '17': {'abbr': 'IL', 'epsg': '3528', 'name': 'Illinois'},
    '18': {'abbr': 'IN', 'epsg': '3532', 'name': 'Indiana'},
    '19': {'abbr': 'IA', 'epsg': '3536', 'name': 'Iowa'},
    '20': {'abbr': 'KS', 'epsg': '3540', 'name': 'Kansas'},
    '21': {'abbr': 'KY', 'epsg': '3544', 'name': 'Kentucky'},
    '22': {'abbr': 'LA', 'epsg': '3550', 'name': 'Louisiana'},
    '23': {'abbr': 'ME', 'epsg': '3557', 'name': 'Maine'},
    '24': {'abbr': 'MD', 'epsg': '3559', 'name': 'Maryland'},
    '25': {'abbr': 'MA', 'epsg': '3585', 'name': 'Massachusetts'},
    '26': {'abbr': 'MI', 'epsg': '3587', 'name': 'Michigan'},
    '27': {'abbr': 'MN', 'epsg': '3594', 'name': 'Minnesota'},
    '28': {'abbr': 'MS', 'epsg': '3597', 'name': 'Mississippi'},
    '29': {'abbr': 'MO', 'epsg': '3602', 'name': 'Missouri'},
    '30': {'abbr': 'MT', 'epsg': '3604', 'name': 'Montana'},
    '31': {'abbr': 'NE', 'epsg': '3606', 'name': 'Nebraska'},
    '32': {'abbr': 'NV', 'epsg': '3607', 'name': 'Nevada'},
    '33': {'abbr': 'NH', 'epsg': '3613', 'name': 'NewHampshire'},
    '34': {'abbr': 'NJ', 'epsg': '3615', 'name': 'NewJersey'},
    '35': {'abbr': 'NM', 'epsg': '3617', 'name': 'NewMexico'},
    '36': {'abbr': 'NY', 'epsg': '3623', 'name': 'NewYork'},
    '37': {'abbr': 'NC', 'epsg': '3631', 'name': 'NorthCarolina'},
    '38': {'abbr': 'ND', 'epsg': '3633', 'name': 'NorthDakota'},
    '39': {'abbr': 'OH', 'epsg': '3637', 'name': 'Ohio'},
    '40': {'abbr': 'OK', 'epsg': '3639', 'name': 'Oklahoma'},
    '41': {'abbr': 'OR', 'epsg': '3645', 'name': 'Oregon'},
    '42': {'abbr': 'PA', 'epsg': '3649', 'name': 'Pennsylvania'},
    '44': {'abbr': 'RI', 'epsg': '3653', 'name': 'RhodeIsland'},
    '45': {'abbr': 'SC', 'epsg': '3655', 'name': 'SouthCarolina'},
    '46': {'abbr': 'SD', 'epsg': '3657', 'name': 'SouthDakota'},
    '47': {'abbr': 'TN', 'epsg': '3661', 'name': 'Tennessee'},
    '48': {'abbr': 'TX', 'epsg': '3669', 'name': 'Texas'},
    '49': {'abbr': 'UT', 'epsg': '3675', 'name': 'Utah'},
    '50': {'abbr': 'VT', 'epsg': '3684', 'name': 'Vermont'},
    '51': {'abbr': 'VA', 'epsg': '3685', 'name': 'Virginia'},
    '53': {'abbr': 'WA', 'epsg': '3689', 'name': 'Washington'},
    '54': {'abbr': 'WV', 'epsg': '3693', 'name': 'WestVirginia'},
    '55': {'abbr': 'WI', 'epsg': '3695', 'name': 'Wisconsin'},
    '56': {'abbr': 'WY', 'epsg': '3703', 'name': 'Wyoming'}
}
