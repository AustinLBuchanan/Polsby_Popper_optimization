import xprgrb as gp
from xprgrb import GRB, solver
from mip_callback import xpress_branch_cb, xpress_cut_cb, xpress_chksol_cb
import math

def add_cut_edges_objective(m, DG):

    # add strengthening vars/constraints:
    # edge {u,v} is cut <=> arc (u,v) is cut <=> arc (v,u) is cut
    undirected_edges = [ (u,v) for u,v in DG.edges if u<v ]
    m._is_cut = m.addVars( undirected_edges, name='iscut', vtype=GRB.BINARY )
    m.addConstrs( m._is_cut[min(u,v),max(u,v)] == gp.quicksum( m._y[u,v,j] for j in range(DG._k) ) for u,v in DG.edges )

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
    m._A = m.addVars(DG._k, name='A')

    # P[j] = perimeter of district j
    m._P = m.addVars(DG._k, name='P')

    # objective is to maximize average Polsby-Popper score
    m.setObjective( ( 1.0 / DG._k ) * coef * gp.quicksum( m._inv_z[j] for j in range(DG._k) ), GRB.MAXIMIZE )

    # add SOCP constraints relating inverse Polsby-Popper score z[j] to area and perimeter
    m.addConstrs( m._P[j] * m._P[j] <= 2 * m._A[j] * m._z[j] for j in range(DG._k) )

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
    m.addConstrs( m._A[j] == gp.quicksum( DG.nodes[i]['area'] * m._x[i,j] for i in DG.nodes ) for j in range(DG._k) )

    # add constraints on perimeters P[j]
    for j in range(DG._k):
        m.addConstr( m._P[j] == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v,j] for u,v in DG.edges )
                 + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i,j] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )

    m.update()
    return


def add_average_Polsby_Popper_objective_binary_expansion(m, DG):
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
    
    ## set priority on largest of binary expansion
    #for i in range(num_digits):
    #    for j in range(DG._k):
    #        b[i,j].BranchPriority = num_digits - i
    
    m.update()
    return
