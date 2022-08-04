import xprgrb as gp
from xprgrb import GRB, solver
from mip_callback import xpress_branch_cb, xpress_cut_cb, xpress_chksol_cb
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

    m._DG = DG
    
    # z[j] / coef is inverse Polsby-Popper score for district j
    coef = 2 * math.pi
    m._obj_coef = coef / DG._k
    m._z = m.addVars(DG._k, name='z', lb=coef )

    # coef * inv_z[j] = coef / z[j] is the Polsby-Popper score for district j
    m._inv_z = m.addVars(DG._k, name='invz', ub=1.0/coef)

    # A[j] = area of district j
    A = m.addVars(DG._k, name='A')

    # P[j] = perimeter of district j
    P = m.addVars(DG._k, name='P')

    # objective is to maximize average Polsby-Popper score
    m.setObjective( ( 1.0 / DG._k ) * coef * gp.quicksum( m._inv_z[j] for j in range(DG._k) ), GRB.MAXIMIZE )

    # add SOCP constraints relating inverse Polsby-Popper score z[j] to area and perimeter
    m.addConstrs( P[j] * P[j] <= 2 * A[j] * m._z[j] for j in range(DG._k) )

    if solver == 'gurobi':

        # impose inv_z = 1 / z through non-convex constraint:
        m.addConstrs( m._z[j] * m._inv_z[j] == 1 for j in range(DG._k) )
        m.Params.NonConvex = 2

    else:

        # Add initial OA cuts for z, inv_z:
        #
        # y >= g(x) ===> y >= g(x*) + g'(x*) (x - x*)
        #
        # inv_z >= 1/z* -1/(z*)^2 * (z - z*)
        # z*^2 * inv_z + z >= 2*z*
        #
        # n_init_OA = 4
        # m.addConstrs(z0**2 * m._inv_z[j] + m._z[j] >= 2*z0 for z0 in [coef*i for i in range(1, n_init_OA + 1)] for j in range(DG._k))

        # Add convex constraint y >= 1/z, which only takes care of one
        # side of the convex envelope
        m.addConstrs(m._inv_z[j] * m._z[j] >= 1 for j in range(DG._k))

        # Establish that z and inv_z should NOT be presolved away
        m.xmodel.loadsecurevecs(rowind=None, colind=[m._z[j]     for j in range(DG._k)] +
                                                    [m._inv_z[j] for j in range(DG._k)])

        assert m.xmodel is not None
        # Do not add nonconvex constraint, but add explicit callbacks
        m.xmodel.addcbpreintsol(xpress_chksol_cb, m, 1)  # Callback for checking if solution is integer
        m.xmodel.addcboptnode(xpress_cut_cb, m, 1)  # Callback for adding lcut inequalities (if specified) and secant cuts for inv_z = 1/z
        m.xmodel.addcbchgbranchobject(xpress_branch_cb, m, 1)  # Callback for branching on the z's

    # add constraints on areas A[j] 
    m.addConstrs( A[j] == gp.quicksum( DG.nodes[i]['area'] * m._x[i,j] for i in DG.nodes ) for j in range(DG._k) )

    # add constraints on perimeters P[j]
    for j in range(DG._k):
        m.addConstr( P[j] == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v,j] for u,v in DG.edges )
                 + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i,j] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )

    m.update()
    return
