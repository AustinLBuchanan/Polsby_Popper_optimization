import gurobipy as gp
from gurobipy import GRB
import math

def add_cut_edges_objective(m, DG):
    # cut edges are double counted, so halve it.
    m.setObjective( 0.5 * gp.quicksum( m._y ), GRB.MINIMIZE )
    return

def add_perimeter_objective(m, DG):
    # total perimeter = internal perimeters + external perimeters
    ip = m.addVar() # internal perimeter
    ep = m.addVar() # external perimeter
    m.addConstr( ip == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v,j] for u,v in DG.edges for j in range(DG._k) ) )
    m.addConstr( ep == gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i,j] for j in range(DG._k) for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )
    m.setObjective( ip + ep, GRB.MINIMIZE )
    return

def add_inverse_Polsby_Popper_objective(m, DG):
                
    # coef * z[j] is inverse Polsby-Popper score for district j
    coef = 1 / ( 2 * math.pi )
    z = m.addVars(DG._k)

    # objective is to minimize average of inverse Polsby-Popper scores
    m.setObjective( ( 1 / DG._k ) * coef * gp.quicksum( z[j] for j in range(DG._k) ), GRB.MINIMIZE )
    
    # A[j] = area of district j
    A = m.addVars(DG._k)

    # P[j] = perimeter of district j
    P = m.addVars(DG._k)
    
    # add SOCP constraints relating inverse Polsby-Popper score z[j] to area and perimeter
    m.addConstrs( P[j] * P[j] <= 2 * A[j] * z[j] for j in range(DG._k) )

    # add constraints on areas A[j] 
    m.addConstrs( A[j] == gp.quicksum( DG.nodes[i]['area'] * m._x[i,j] for i in DG.nodes ) for j in range(DG._k) )

    # add constraints on perimeters P[j]
    for j in range(DG._k):
        m.addConstr( P[j] == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v,j] for u,v in DG.edges )
                 + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i,j] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )

    m.update()
    return

def add_average_Polsby_Popper_objective(m, DG):
    
    # z[j] / coef is inverse Polsby-Popper score for district j
    coef = 2 * math.pi 
    z = m.addVars(DG._k)

    # coef * inv_z[j] = coef / z[j] is the Polsby-Popper score for district j
    inv_z = m.addVars(DG._k) 

    # A[j] = area of district j
    A = m.addVars(DG._k)

    # P[j] = perimeter of district j
    P = m.addVars(DG._k)
    
    # objective is to maximize average Polsby-Popper score
    m.setObjective( ( 1 / DG._k ) * coef * gp.quicksum( inv_z[j] for j in range(DG._k) ), GRB.MAXIMIZE )

    # add SOCP constraints relating inverse Polsby-Popper score z[j] to area and perimeter
    m.addConstrs( P[j] * P[j] <= 2 * A[j] * z[j] for j in range(DG._k) )

    # impose inv_z = 1 / z through non-convex constraint:
    m.addConstrs( z[j] * inv_z[j] == 1 for j in range(DG._k) )

    # add constraints on areas A[j] 
    m.addConstrs( A[j] == gp.quicksum( DG.nodes[i]['area'] * m._x[i,j] for i in DG.nodes ) for j in range(DG._k) )

    # add constraints on perimeters P[j]
    for j in range(DG._k):
        m.addConstr( P[j] == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v,j] for u,v in DG.edges )
                 + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i,j] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )

    m.Params.NonConvex = 2
    m.update()
    return
