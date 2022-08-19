import xprgrb as gp
from xprgrb import GRB
from mip_callback import xpress_branch_cb, xpress_cut_cb, xpress_chksol_cb
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

    # coef * z[j] is inverse Polsby-Popper score for district j

    z = m.addVars(DG._k, name='z')#, lb=1)

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


def find_bounds(DG):

    from xprgrb import solver

    m = gp.Model()

    # Create variables
    # x[i] equals one when node i is assigned to district j
    m._x = m.addVars(DG.nodes, name='x', vtype=GRB.BINARY)

    # r[i] equals one when node i roots district j
    m._r = m.addVars(DG.nodes, name='r', vtype=GRB.BINARY)

    # y[u,v] equals one when arc (u,v) is cut because u->j but not v->j
    if solver == 'gurobi':
        m._y = m.addVars(DG.edges, name='y', vtype=GRB.BINARY)
    else:
        m._y = {(u,v): m.addVar(name=f'y_{u}_{v}', vtype=GRB.BINARY) for (u,v) in DG.edges}

    # add constraints saying that if node i roots district j, then i should be in district j
    m.addConstrs( m._r[i] <= m._x[i] for i in DG.nodes)

    # district has population at least L and at most U
    m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i] for i in DG.nodes) >= DG._L)
    m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i] for i in DG.nodes) <= DG._U)

    # add constraints saying that edge {u,v} is cut if u is assigned to district j but v is not.
    m.addConstrs( m._x[u] - m._x[v] <= m._y[u,v] for u,v in DG.edges)

    # add strengthening vars/constraints:
    # edge {u,v} is cut <=> arc (u,v) is cut <=> arc (v,u) is cut
    #undirected_edges = [ (u,v) for u,v in DG.edges if u<v ]
    #m._is_cut = m.addVars( undirected_edges, name='iscut', vtype=GRB.BINARY )
    #m.addConstrs( m._is_cut[min(u,v),max(u,v)] == gp.quicksum( m._y[u,v,j] for j in range(DG._k) ) for u,v in DG.edges )

    # g[i,j] = amount of flow generated at node i of type j
    g = m.addVars(DG.nodes, name='g')

    # f[j,u,v] = amount of flow sent across arc uv of type j
    if solver == 'gurobi':
        f = m.addVars(DG.edges, name='f' )
    else:
        f = m.addVars([(u,v) for (u,v) in DG.edges], name='f')

    # compute big-M
    M = mip_contiguity.most_possible_nodes_in_one_district(DG) - 1

    # flow can only be generated at roots
    m.addConstrs( g[i] <= (M+1)*m._r[i] for i in DG.nodes)

    # flow balance
    m.addConstrs( g[i] - m._x[i] == gp.quicksum( f[i,u]-f[u,i] for u in DG.neighbors(i)) for i in DG.nodes)

    # flow type j can enter vertex i only if (i is assigned to district j) and (i is not root of j)
    m.addConstrs( gp.quicksum( f[u,i] for u in DG.neighbors(i) ) <= M * (m._x[i]-m._r[i]) for i in DG.nodes)

    m.Params.TimeLimit = 60
    m.Params.OutputFlag = 0
    # Obtain lower/upper bounds on A ########################

    m.setObjective(gp.quicksum( DG.nodes[i]['area'] * m._x[i] for i in DG.nodes), GRB.MAXIMIZE)
    m.update()
    m.optimize()
    Au = m.objVal

    m.setObjective(gp.quicksum( DG.nodes[i]['area'] * m._x[i] for i in DG.nodes), GRB.MINIMIZE)
    m.update()
    m.optimize()
    Al = m.objVal

    # Obtain lower/upper bounds on P ########################

    m.setObjective (gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v] for u,v in DG.edges )
                  + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i] for i in DG.nodes if DG.nodes[i]['boundary_node']), GRB.MAXIMIZE)
    m.update()
    m.optimize()
    Pu = m.objVal

    m.setObjective (gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v] for u,v in DG.edges )
                  + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i] for i in DG.nodes if DG.nodes[i]['boundary_node']), GRB.MINIMIZE)
    m.update()
    m.optimize()
    Pl = m.objVal

    if DG.options['zbounds'] == 'no':

        zl = Pl**2 / (4 * math.pi * Au)
        zu = Pu**2 / (4 * math.pi * Al)

    else:

        z = m.addVar(name='z')
        A = m.addVar(name='A')
        P = m.addVar(name='P')

        m.addConstr(A == gp.quicksum( DG.nodes[i]['area'] * m._x[i] for i in DG.nodes))
        m.addConstr(P == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v] for u,v in DG.edges )
                    + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i] for i in DG.nodes if DG.nodes[i]['boundary_node']))

        m.addConstr(P**2 <= 4 * math.pi * A * z)

        m.setObjective(z, GRB.MINIMIZE)
        m.update()
        m.Params.TimeLimit = 30
        m.optimize()
        zl = m.objBound

        m.addConstr(P**2 == 4 * math.pi * A * z)

        m.setObjective(z, GRB.MAXIMIZE)
        m.update()
        m.Params.NonConvex = 2
        m.Params.TimeLimit = 30
        m.optimize()
        zu = m.objBound

    return Al, Au, Pl, Pu, zl, zu


def add_average_Polsby_Popper_objective(m, DG):

    m._DG = DG

    # z[j] / coef is inverse Polsby-Popper score for district j
    # coef = 2 * math.pi

    print("Finding tight bounds on P&A")
    # Find lower/upper bounds for A and P by solving four smaller problems
    (Al, Au, Pl, Pu, zlb, zub) = find_bounds(DG)

    izlb = 1.0 / zub
    izub = 1.0 / zlb

    m._obj_coef = 1 / DG._k
    m._z = m.addVars(DG._k, name='z', lb=zlb, ub=zub)

    # coef * inv_z[j] = coef / z[j] is the Polsby-Popper score for district j
    m._inv_z = m.addVars(DG._k, name='invz', lb=izlb, ub=izub)

    # A[j] = area of district j
    m._A = m.addVars(DG._k, name='A', lb=Al, ub=Au)

    # P[j] = perimeter of district j
    m._P = m.addVars(DG._k, name='P', lb=Pl, ub=Pu)

    # objective is to maximize average Polsby-Popper score
    m.setObjective( ( 1.0 / DG._k ) * gp.quicksum( m._inv_z[j] for j in range(DG._k) ), GRB.MAXIMIZE )

    # add SOCP constraints relating inverse Polsby-Popper score z[j] to area and perimeter
    m.addConstrs( m._P[j] * m._P[j] <= 4 * math.pi * m._A[j] * m._z[j] for j in range(DG._k) )

    from xprgrb import solver

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


def add_average_Schwartzberg_objective(m, DG):

    decompose = DG._options['conedecomp'] == 'yes'  # decompose SOC using VDHL

    print("Finding tight bounds on P&A")
    # Find lower/upper bounds for A and P by solving four smaller problems
    (Al, Au, Pl, Pu, zlb, zub) = find_bounds(DG)

    z = m.addVars(DG._k, name='z')  #,     lb=zlb,            ub=zub)
    s = m.addVars(DG._k, name='sroot')  #, lb=math.sqrt(zlb), ub=math.sqrt(zub))

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

    rho = zub / (1 - 2.0**-num_digits)

    # b[t,j] is digit t in binary expansion of z[j]
    b = m.addVars(num_digits, DG._k, name='b', vtype=GRB.BINARY)

    # Define z as the binary expansion using b[] variables
    m.addConstrs(z[j] == rho * gp.quicksum(2**(-(1 + t)) * b[t,j] for t in range(num_digits)) for j in range(DG._k))

    if decompose:

        # Define the variables used in the Vielma/Dunning/Huchette/Lubin cone decomposition
        u = m.addVars(num_digits, DG._k, name='cdec')

        # Constraints for the cone decomposition
        m.addConstrs(2*rho*gp.quicksum(2**-t * u[t,j] for t in range(num_digits)) <= s[j] for j in range(DG._k))
        m.addConstrs(b[t,j]**2 <= 2*u[t,j]*s[j] for t in range(num_digits) for j in range(DG._k))

    else:
        m.addConstrs(rho * gp.quicksum(2**-(t + 1) * b[t,j]**2 for t in range(num_digits)) <= s[j]**2 for j in range(DG._k))

    m.update()
