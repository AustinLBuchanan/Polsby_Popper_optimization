import gurobipy as gp
from gurobipy import GRB

def build_base_mip(DG):
    
    # Create model object
    m = gp.Model()
    
    # Create variables
    # x[i,j] equals one when node i is assigned to district j
    m._x = m.addVars(DG.nodes, DG._k, vtype=GRB.BINARY)  
    
    # r[i,j] equals one when node i roots district j
    m._r = m.addVars(DG.nodes, DG._k, vtype=GRB.BINARY)  

    # y[u,v,j] equals one when arc (u,v) is cut because u->j but not v->j
    m._y = m.addVars(DG.edges, DG._k, vtype=GRB.BINARY) 
    
    # add constraints saying that each node i is assigned to one district
    m.addConstrs( gp.quicksum( m._x[i,j] for j in range(DG._k)) == 1 for i in DG.nodes )

    # add constraints saying that each district has population at least L and at most U
    m.addConstrs( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i,j] for i in DG.nodes) >= DG._L for j in range(DG._k) )
    m.addConstrs( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i,j] for i in DG.nodes) <= DG._U for j in range(DG._k) )

    # add constraints saying that edge {u,v} is cut if u is assigned to district j but v is not.
    m.addConstrs( m._x[u,j] - m._x[v,j] <= m._y[u,v,j] for u,v in DG.edges for j in range(DG._k) )

    m.update()
    return m
    
def add_partitioning_orbitope_constraints(m, DG):
    
    s = m.addVars(DG.nodes, DG._k)
    u = m.addVars(DG.nodes, DG._k)
    w = m.addVars(DG.nodes, DG._k) 

    m.addConstrs(m._x[i,j] == s[i,j]-s[i,j+1] for i in DG.nodes for j in range(DG._k-1))
    m.addConstrs(m._x[i,DG._k-1] == s[i,DG._k-1] for i in DG.nodes)

    m.addConstrs(m._r[DG._ordering[0],j] == w[DG._ordering[0],j] for j in range(DG._k))
    m.addConstrs(m._r[DG._ordering[i],j] == w[DG._ordering[i],j] - w[DG._ordering[i-1],j] for i in range(1,DG.number_of_nodes()) for j in range(DG._k))

    m.addConstrs(m._r[i,j] <= m._x[i,j] for i in DG.nodes for j in range(DG._k))
    m.addConstrs(s[i,j] <= w[i,j] for i in DG.nodes for j in range(DG._k))

    m.addConstrs(u[DG._ordering[i],j]+m._r[DG._ordering[i],j] == u[DG._ordering[i+1],j] + m._r[DG._ordering[i+1],j+1] for i in range(0,DG.number_of_nodes()-1) for j in range(DG._k-1))
    m.addConstrs(u[DG._ordering[i],DG._k-1]+m._r[DG._ordering[i],DG._k-1] == u[DG._ordering[i+1],DG._k-1] for i in range(0,DG.number_of_nodes()-1))
    m.addConstrs(u[DG._ordering[DG.number_of_nodes()-1],j]+m._r[DG._ordering[DG.number_of_nodes()-1],j] == 0 for j in range(DG._k-1))

    m._r[DG._ordering[0],0].LB=1
    m.addConstr( u[DG._ordering[DG.number_of_nodes()-1],DG._k-1] + m._r[DG._ordering[DG.number_of_nodes()-1],DG._k-1]==1 )  
    m.update()
    return
    
def inject_warm_start(m, DG, heuristic_labeling):
    
    root = { j : -1 for j in range(DG._k) }
    for i in DG._ordering:
        for j in range(DG._k):
            m._x[i,j].start = 0
            m._r[i,j].start = 0
            
    for u,v in DG.edges:
        for j in range(DG._k):
            m._y[u,v,j].start = 0
    
    for i in DG._ordering:
        j = heuristic_labeling[i]
        m._x[i,j].start = 1
        if root[j] == -1:
            root[j] = i
            m._r[i,j].start = 1
    
    for u,v in DG.edges:
        j = heuristic_labeling[u]
        if heuristic_labeling[v] != j:
            m._y[u,v,j].start = 1
    
#     if contiguity == 'scf':
        
#         #f[u,v]=0 when u and v belong to different districts
#         for u,v in DG.edges:
#             if heuristic_labeling[u] != heuristic_labeling[v]:
#                 m._f[u,v].start = 0
        
#         # pick a spanning tree of G[district]. other edges have f[u,v]=0
#         for j in range(DG._k):
#             district = [ i for i in DG.nodes if heuristic_labeling[i] == j ]
#             tree_edges = list( nx.dfs_edges(DG.subgraph(district), source=root[j]) ) 
#             for u in district:
#                 for v in DG.neighbors(u):
#                     if v in district and (u,v) not in tree_edges:
#                         m._f[u,v].start = 0
    return