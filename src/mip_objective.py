import gurobipy as gp
from gurobipy import GRB
import math
import mip_contiguity

def add_cut_edges_objective(m, DG):
    m.setObjective( gp.quicksum( m._is_cut ), GRB.MINIMIZE )
    return


def add_perimeter_objective(m, DG):
    # minimize total perimeter = external perimeter + internal perimeter
    ep = sum( DG.nodes[i]['boundary_perim'] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) 
    m.setObjective( ep + gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v,j] for u,v in DG.edges for j in range(DG._k) ), GRB.MINIMIZE )
    return


def add_inverse_Polsby_Popper_objective(m, DG):

    # z[j] is inverse Polsby-Popper score for district j
    z = m.addVars(DG._k, name='z')

    # objective is to minimize average of inverse Polsby-Popper scores
    m.setObjective( ( 1.0 / DG._k ) * gp.quicksum( z[j] for j in range(DG._k) ), GRB.MINIMIZE )
    
    # A[j] = area of district j
    A = m.addVars(DG._k, name='A')

    # P[j] = perimeter of district j
    P = m.addVars(DG._k, name='P')

    # add SOCP constraints relating inverse Polsby-Popper score z[j] to area and perimeter
    m.addConstrs( P[j] * P[j] <= 4 * math.pi * A[j] * z[j] for j in range(DG._k) )

    # add constraints on areas A[j] 
    m.addConstrs( A[j] == gp.quicksum( DG.nodes[i]['area'] * m._x[i,j] for i in DG.nodes ) for j in range(DG._k) )

    # add constraints on perimeters P[j]
    for j in range(DG._k):
        m.addConstr( P[j] == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v,j] for u,v in DG.edges )
                 + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i,j] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )

    m.update()
    return


def add_average_Polsby_Popper_objective(m, DG):
    m._DG = DG
    
    # We use a big-M coefficient in the binary expansion,
    #   to linearize the product b[t,j] * z[j] where
    #   b[t,j] is binary and z[j] is continuous.
    # We assume that each district D_j has
    #             PP(D_j) >= 1 / big_M, i.e.,
    #   z[j] = PP^-1(D_j) <= big_M. 
    big_M = 16      
    
    # number of digits in binary expansion (controls precision)
    num_digits = 20 
    
    # z[j] is inverse Polsby-Popper score for district j
    m._z = m.addVars(DG._k, name='z') #, lb=1 )
    
    # b[t,j] is digit t in binary expansion of 1/z[j]
    b = m.addVars(num_digits, DG._k, name='b', vtype=GRB.BINARY)
    
    # w[t,j] = b[t,j] * z[j] is used for linearization purposes
    w = m.addVars(num_digits, DG._k, name='w')

    # A[j] = area of district j
    A = m.addVars(DG._k, name='A')

    # P[j] = perimeter of district j
    P = m.addVars(DG._k, name='P')

    # objective is to maximize average Polsby-Popper score
    m.setObjective( ( 1.0 / DG._k ) * gp.quicksum( (1/(2**(i+1))) * b[i,j] for j in range(DG._k) for i in range(num_digits) ), GRB.MAXIMIZE )

    # add SOCP constraints relating inverse Polsby-Popper score z[j] to area and perimeter
    m.addConstrs( P[j] * P[j] <= 4 * math.pi * A[j] * m._z[j] for j in range(DG._k) )

    # add constraints on areas A[j] 
    m.addConstrs( A[j] == gp.quicksum( DG.nodes[i]['area'] * m._x[i,j] for i in DG.nodes ) for j in range(DG._k) )

    # add constraints on perimeters P[j]
    for j in range(DG._k):
        m.addConstr( P[j] == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v,j] for u,v in DG.edges )
                 + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i,j] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )
        
    ###################################    
    # add binary expansion constraints 
    ###################################
    
    # Impose that 1 / z[j] = 0.5 * b[0,j] + 0.25 * b[1,j] + 0.125 * b[2,j] + ... 
    # Do this by multiplying both sides by z[j] and then replacing b[i,j] * z[j] by w[i,j].
    # Also, relax equality to 1 / z[j] >= ... to reduce numerical troubles.
    m.addConstrs( 1 >= gp.quicksum( (1/(2**(i+1))) * w[i,j] for i in range(num_digits) ) for j in range(DG._k) )
    
    # Impose w[i,j] <= b[i,j] * z[j]:
    m.addConstrs( w[i,j] <= big_M * b[i,j] for i in range(num_digits) for j in range(DG._k) )
    m.addConstrs( w[i,j] <= m._z[j] for i in range(num_digits) for j in range(DG._k) )
    
    # Impose w[i,j] >= b[i,j] * z[j]:
    m.addConstrs( m._z[j] + big_M * ( b[i,j] - 1 ) <= w[i,j] for i in range(num_digits) for j in range(DG._k) )
    
    m.update()
    return


def add_average_Schwartzberg_objective(m, DG):

    z = m.addVars(DG._k, name='z')  
    s = m.addVars(DG._k, name='sroot')  

    # objective is to minimize average of inverse Polsby-Popper scores
    m.setObjective( ( 1.0 / DG._k ) * gp.quicksum( s[j] for j in range(DG._k) ), GRB.MINIMIZE )

    # A[j] = area of district j
    A = m.addVars(DG._k, name='A')

    # P[j] = perimeter of district j
    P = m.addVars(DG._k, name='P')

    # add SOCP constraints relating inverse Polsby-Popper score z[j] to area and perimeter
    m.addConstrs( P[j] * P[j] <= 4 * math.pi * A[j] * z[j] for j in range(DG._k) )

    # add constraints on areas A[j]
    m.addConstrs( A[j] == gp.quicksum( DG.nodes[i]['area'] * m._x[i,j] for i in DG.nodes ) for j in range(DG._k) )

    # add constraints on perimeters P[j]
    for j in range(DG._k):
        m.addConstr( P[j] == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v,j] for u,v in DG.edges )
                 + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i,j] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )

    # Binary expansion variables
    num_digits = 20
    big_M = 17

    # b[t,j] is digit t in binary expansion of z[j]
    b = m.addVars(num_digits, DG._k, name='b', vtype=GRB.BINARY)

    # Define z as the binary expansion using b[] variables
    m.addConstrs( z[j] <= 1 + (big_M-1) * gp.quicksum(2**(-(1 + t)) * b[t,j] for t in range(num_digits)) for j in range(DG._k) )

    m.addConstrs( 1 + (big_M-1) * gp.quicksum(2**(-(1 + t)) * b[t,j] * b[t,j] for t in range(num_digits)) <= s[j]**2 for j in range(DG._k))
        
    m.update()
