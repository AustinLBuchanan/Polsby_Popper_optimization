import networkx as nx
import gurobipy as gp

# main fixing function, that calls the others:
def do_variable_fixing(m, DG):
    do_diagonal_fixing(m, DG)
    do_L_fixing(m, DG)
    do_U_fixing(m, DG)
    do_cut_edge_fixing(m, DG)
    return

# subroutines:
def do_diagonal_fixing(m, DG):
    for p in range(DG.number_of_nodes()):
        i = DG._ordering[p]
        for j in range(p+1,DG._k):
            m._x[i,j].UB = 0
    return
            
def do_L_fixing(m, DG):
    #   First, find "back" of ordering B = {v_q, v_{q+1}, ..., v_{n-1} }
    n = DG.number_of_nodes()
    S = [ False for v in DG.nodes ]
    for p in range(n):
        v_pos = n - p - 1
        v = DG._ordering[v_pos]
        S[v] = True
        pr = reachable_population(DG, S, v)
        if pr >= DG._L:
            q = v_pos + 1
            break

    # none of the vertices at back (in B) can root a district. 
    for p in range(q,n):
        i = DG._ordering[p]
        for j in range(DG._k):
            m._r[i,j].UB = 0

    # vertex v_{q-1} cannot root districts {0, 1, ..., k-2}
    # vertex v_{q-2} cannot root districts {0, 1, ..., k-3}
    # ... 
    # vertex v_{q-t} cannot root districts {0, 1, ..., k-t-1}
    # ...
    # vertex v_{q-(k-1)} cannot root district {0}
    for t in range(1,DG._k):
        i = DG._ordering[q-t]
        for j in range(DG._k-t):
            m._r[i,j].UB = 0
    
    m.update()
    return
      
def do_U_fixing(m, DG):
    for (i,j) in DG.edges:
        DG[i][j]['ufixweight'] = DG.nodes[j]['TOTPOP'] # weight of edge (i,j) is population of its head j

    for j in range(DG._k):

        v = DG._ordering[j]
        dist = nx.shortest_path_length( DG, source=v, weight='ufixweight' )

        if j == 0:
            min_dist = DG._U+1
        else:
            min_dist = min( dist[DG._ordering[t]] + DG.nodes[v]['TOTPOP'] for t in range(j) )

        if min_dist <= DG._U:
            break

        m._r[v,j].LB = 1
        m._x[v,j].LB = 1

        for t in range(DG._k):
            if t != j:
                m._x[v,t].UB = 0

        for i in DG.nodes:
            if i != v:
                m._r[i,j].UB = 0

        for i in DG.nodes:
            if i != v and dist[i] + DG.nodes[v]['TOTPOP'] > DG._U:
                m._x[i,j].UB = 0
    
    return

# how many people are reachable from v in DG[S]? Uses BFS
def reachable_population(DG, S, v):
    pr = 0 # population reached
    if not S[v]:
        return 0
    
    visited = { i : False for i in DG.nodes }
    child = [v]
    visited[v] = True
    while child:
        parent = child
        child = list()
        for i in parent:
            pr += DG.nodes[i]['TOTPOP']
            for j in DG.neighbors(i):
                if S[j] and not visited[j]:
                    child.append(j)
                    visited[j] = True
    return pr   

def do_cut_edge_fixing(m, DG):
    for u,v in DG.edges:
        for j in range(DG._k):
            if m._x[u,j].UB < 0.5 or m._x[v,j].LB > 0.5:
                m._y[u,v,j].UB = 0
            elif m._x[u,j].LB > 0.5 and m._x[v,j].UB < 0.5:
                m._y[u,v,j].LB = 1
    return
