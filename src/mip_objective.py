import xprgrb as gp
from xprgrb import GRB, solver
from mip_callback import xpress_branch_cb, xpress_cut_cb
import math

def add_cut_edges_objective(m, DG):
    m.setObjective( gp.quicksum( m._is_cut ), GRB.MINIMIZE )
    return

def add_perimeter_objective(m, DG):
    # minimize total perimeter = external perimeter + internal perimeter
    ep = sum( DG.nodes[i]['boundary_perim'] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) 
    m.setObjective( ep + gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v,j] for u,v in DG.edges for j in range(DG._k) ), GRB.MINIMIZE )
    return

def add_inverse_Polsby_Popper_objective(m, DG):
                
    # coef * z[j] is inverse Polsby-Popper score for district j

    coef = 1.0 / ( 2 * math.pi )
    z = m.addVars(DG._k, name='z')#, lb=2*math.pi)

    # objective is to minimize average of inverse Polsby-Popper scores
    m.setObjective( ( 1.0 / DG._k ) * coef * gp.quicksum( z[j] for j in range(DG._k) ), GRB.MINIMIZE )
    
    # A[j] = area of district j
    A = m.addVars(DG._k, name='A')

    # P[j] = perimeter of district j
    P = m.addVars(DG._k, name='P')

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
    z = m.addVars(DG._k, name='z', lb=coef )

    # coef * inv_z[j] = coef / z[j] is the Polsby-Popper score for district j
    inv_z = m.addVars(DG._k, name='invz', ub=1.0/coef)

    # A[j] = area of district j
    A = m.addVars(DG._k, name='A')

    # P[j] = perimeter of district j
    P = m.addVars(DG._k, name='P')

    # objective is to maximize average Polsby-Popper score
    m.setObjective( ( 1.0 / DG._k ) * coef * gp.quicksum( inv_z[j] for j in range(DG._k) ), GRB.MAXIMIZE )

    # add SOCP constraints relating inverse Polsby-Popper score z[j] to area and perimeter
    m.addConstrs( P[j] * P[j] <= 2 * A[j] * z[j] for j in range(DG._k) )

    if solver == 'gurobi':

        # impose inv_z = 1 / z through non-convex constraint:
        m.addConstrs( z[j] * inv_z[j] == 1 for j in range(DG._k) )
        m.Params.NonConvex = 2

    else:
        pass
        assert m.xmodel is not None
        # Do not add nonconvex constraint, but add explicit callbacks
        m.xmodel.addcbintsol(xpress_chksol_cb, DG, 1)
        m.xmodel.addcboptnode(xpress_cut_cb, DG, 1)
        m.xmodel.addcbchgbranchobject(xpress_branch_cb, DG, 1)

    # add constraints on areas A[j] 
    m.addConstrs( A[j] == gp.quicksum( DG.nodes[i]['area'] * m._x[i,j] for i in DG.nodes ) for j in range(DG._k) )

    # add constraints on perimeters P[j]
    for j in range(DG._k):
        m.addConstr( P[j] == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v,j] for u,v in DG.edges )
                 + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i,j] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )

    m.update()
    return
