import xprgrb as gp
from xprgrb import GRB, solver

def build_base_mip(DG):
    
    # Create model object
    m = gp.Model()
    
    # Create variables
    # x[i,j] equals one when node i is assigned to district j
    m._x = m.addVars(DG.nodes, DG._k, name='x', vtype=GRB.BINARY)

    # r[i,j] equals one when node i roots district j
    m._r = m.addVars(DG.nodes, DG._k, name='r', vtype=GRB.BINARY)

    # y[u,v,j] equals one when arc (u,v) is cut because u->j but not v->j
    if solver == 'gurobi':
        m._y = m.addVars(DG.edges, DG._k, name='y', vtype=GRB.BINARY)
    else:
        m._y = {(u,v,j): m.addVar(name=f'y_{u}_{v}_{j}', vtype=GRB.BINARY) for (u,v) in DG.edges for j in range(DG._k)}

    # add constraints saying that each node i is assigned to one district
    m.addConstrs( gp.quicksum( m._x[i,j] for j in range(DG._k)) == 1 for i in DG.nodes )
    
    # add constraints saying that if node i roots district j, then i should be in district j
    m.addConstrs( m._r[i,j] <= m._x[i,j] for i in DG.nodes for j in range(DG._k) )

    # add constraints saying that each district has population at least L and at most U
    m.addConstrs( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i,j] for i in DG.nodes) >= DG._L for j in range(DG._k) )
    m.addConstrs( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i,j] for i in DG.nodes) <= DG._U for j in range(DG._k) )

    # add constraints saying that edge {u,v} is cut if u is assigned to district j but v is not.
    m.addConstrs( m._x[u,j] - m._x[v,j] <= m._y[u,v,j] for u,v in DG.edges for j in range(DG._k) )

    m.update()
    return m
    
def add_partitioning_orbitope_constraints(m, DG):

    s = m.addVars(DG.nodes, DG._k, name='s')
    u = m.addVars(DG.nodes, DG._k, name='u')
    w = m.addVars(DG.nodes, DG._k, name='w')

    m.addConstrs(m._x[i,j] == s[i,j]-s[i,j+1] for i in DG.nodes for j in range(DG._k-1))
    m.addConstrs(m._x[i,DG._k-1] == s[i,DG._k-1] for i in DG.nodes)

    m.addConstrs(m._r[DG._ordering[0],j] == w[DG._ordering[0],j] for j in range(DG._k))
    m.addConstrs(m._r[DG._ordering[i],j] == w[DG._ordering[i],j] - w[DG._ordering[i-1],j] for i in range(1,DG.number_of_nodes()) for j in range(DG._k))

    m.addConstrs(s[i,j] <= w[i,j] for i in DG.nodes for j in range(DG._k))

    m.addConstrs(u[DG._ordering[i],j]+m._r[DG._ordering[i],j] == u[DG._ordering[i+1],j] + m._r[DG._ordering[i+1],j+1] for i in range(0,DG.number_of_nodes()-1) for j in range(DG._k-1))
    m.addConstrs(u[DG._ordering[i],DG._k-1]+m._r[DG._ordering[i],DG._k-1] == u[DG._ordering[i+1],DG._k-1] for i in range(0,DG.number_of_nodes()-1))
    m.addConstrs(u[DG._ordering[DG.number_of_nodes()-1],j]+m._r[DG._ordering[DG.number_of_nodes()-1],j] == 0 for j in range(DG._k-1))

    gp.setLB(m._r[DG._ordering[0],0], m, 1)
    m.addConstr( u[DG._ordering[DG.number_of_nodes()-1],DG._k-1] + m._r[DG._ordering[DG.number_of_nodes()-1],DG._k-1]==1 )  
    m.update()
    return
    
# We need to be careful because of partitioning orbitope constraints. 
#   For each district, find its earliest vertex in the ordering,
#   then sort these earliest vertices to get the district labels
#    
def get_orbitope_friendly_labeling(DG, unfriendly_labeling):    
    
    district_map = { j : -1 for j in range(DG._k) }
    labeling = { i : -1 for i in DG.nodes }
    count = 0
    
    for i in DG._ordering:
        j = unfriendly_labeling[i]
        
        # have we found earliest vertex from district j?
        if district_map[j] == -1:
            
            # if so, then vertex i roots district 'j' in unmapped labeling,
            #   and anything labeled 'j' should instead be relabeled 'count'
            district_map[j] = count
            count += 1
        
        labeling[i] = district_map[j]
        
    return labeling

    
def inject_warm_start(m, DG, labeling):
    global solver

    if solver == 'xpress':
        return inject_warm_start_xpress(m, DG, labeling)

    #print("In inject, labeling =",labeling)
    
    # initialize all variables to 0
    for i in DG._ordering:
        for j in range(DG._k):
            m._x[i,j].start = 0
            m._r[i,j].start = 0
            
    for u,v in DG.edges:
        for j in range(DG._k):
            m._y[u,v,j].start = 0
    
    # now inject the nonzeros of our solution
    root_found = { j : False for j in range(DG._k) }
    for i in DG._ordering:
        
        j = labeling[i]
        m._x[i,j].start = 1
        
        if not root_found[j]:
            root_found[j] = True
            m._r[i,j].start = 1
        
    for u,v in DG.edges:
        j = labeling[u]
        if labeling[v] != j:
            m._y[u,v,j].start = 1
    
    m.update()
    return


def inject_warm_start_xpress(m, DG, labeling):

    sol = {}

    # initialize all variables to 0
    for i in DG._ordering:
        for j in range(DG._k):
            sol['x', i, j] = 0
            sol['r', i, j] = 0

    for u,v in DG.edges:
        for j in range(DG._k):
            sol['y', u, v, j] = 0

    # now inject the nonzeros of our solution
    root_found = { j : False for j in range(DG._k) }
    for i in DG._ordering:

        j = labeling[i]
        sol['x', i, j] = 1

        if not root_found[j]:
            root_found[j] = True
            sol['r', i, j] = 1

    for u,v in DG.edges:
        j = labeling[u]
        if labeling[v] != j:
            sol['y', u, v, j] = 1

    solval = [sol['x', i, j]    for i in DG._ordering for j in range(DG._k)] + \
             [sol['r', i, j]    for i in DG._ordering for j in range(DG._k)] + \
             [sol['y', u, v, j] for u, v in DG.edges  for j in range(DG._k)]

    solind = [m._x[i,j]   for i in DG._ordering for j in range(DG._k)] + \
             [m._r[i,j]   for i in DG._ordering for j in range(DG._k)] + \
             [m._y[u,v,j] for u,v in DG.edges   for j in range(DG._k)]

    if m._objective == 'avepp':
        m._stored_solutions.append((0, solval, solind))
    else:
        m.xmodel.addmipsol(solval, solind, name='mysol')
