import xprgrb as gp
from xprgrb import GRB, solver
import xpress as xp
import numpy as np
import mip_contiguity
import math


def gen_callback(m, where):


def xpress_cut_cb(prob, m):
    """
    Wrapper for the secant+OA separator and the lcut separator
    """

    # Bit 7 is set if the solution is valid
    if (prob.attributes.presolvestate & 128) == 0:
        return 0

    try:
        sol,lb,ub = get_sol_and_bounds(prob, m)
    except:
        return 0

    # Call two cut separators: for nonconvex constraint z * inv_z == 1
    # and for the lcut objective

    r1, r2 = 0, 0
    if m._objective == 'avepp':
        r1 = xpress_cut_nonconvex(prob, m, sol, lb, ub)
    if m._contiguity == 'lcut':
        r2 = mip_contiguity.lcut_callback_xpress_optnode(prob, m, sol, lb, ub)

    return r1 + r2


def get_sol_and_bounds(prob, m):
    """
    return non-presolved x, lb, ub of current node
    """

    try:
        x = []
        prob.getlpsol(x)  # postsolved solution
    except:
        return (None, None, None)

    x = np.array(x)

    lb = []
    ub = []
    prob.getlb(lb, 0, prob.attributes.cols - 1)
    prob.getub(ub, 0, prob.attributes.cols - 1)

    colmap = []
    rowmap = []
    prob.getpresolvemap(rowmap, colmap)
    colmap = np.array(colmap)

    n = len(x)

    olb = np.full(n, -1e10)
    oub = np.full(n,  1e10)
    olb[colmap] = lb
    oub[colmap] = ub

    return x, olb, oub


def get_secant_cut(lb, ub):
    """Given a function y = 1/x with x in [lb,ub], return the
    coefficients (for x, y in this order) and the rhs secant cut for
    the function, with sign <=
    """
    return ([1/lb - 1/ub, ub - lb], ub / lb - lb / ub)


def xpress_branch_cb(prob, m, branch):
    """This callback adds a branch object on z or inv_z if the current
    relaxation cannot be tightened via secant cuts (it is always
    possible to tighten one via OA cuts)"
    """

    DG = m._DG

    # Otherwise, let's see
    try:
        sol,lb,ub = get_sol_and_bounds(prob, m)
    except:
        return branch

    indx = np.array([m.xmodel.getIndex(m._x[i,j]) for i in DG.nodes for j in range(DG._k)])
    #indy = np.array([m.xmodel.getIndex(m._y[u,v,j]) for u,v in DG.edges for j in range(DG._k)])

    x = sol[indx]

    indz  = np.array([m.xmodel.getIndex(m._z[i])     for i in range(DG._k)])
    indiz = np.array([m.xmodel.getIndex(m._inv_z[i]) for i in range(DG._k)])

    # weight of the LP point vs. midpoint for branching
    alpha = .8

    zval = sol[indz]
    zlb  = lb[indz]
    zub  = ub[indz]

    midpoint = 0.5 * (zlb + zub)
    midpoint[zub > 1e5] = zval[zub > 1e5]

    brpt = alpha * zval + (1 - alpha) * midpoint

    invz = sol[indiz]

    prods = zval * invz

    # If all points are below the hyperbola branch x*y==1, no need to
    # branch on the nonconvexity
    if np.max(prods) <= 1 + prob.controls.feastol:
        return branch

    i_viol = np.argmax(prods)

    assert prods[i_viol] > 1 + prob.controls.feastol

    br = xp.branchobj(prob, isoriginal=True)

    br.addbranches(2)

    indices = [indz[i_viol], indiz[i_viol]]

    br.addbounds(0, ['U'], [indz[i_viol]],  [brpt[i_viol]])
    br.addbounds(0, ['L'], [indiz[i_viol]], [1.0 / brpt[i_viol]])

    if zlb[i_viol] > 1e-8:
        coeffs, rhs = get_secant_cut(zlb[i_viol], brpt[i_viol])
        br.addrows(0, ['L'], [rhs], [0, 2], indices, coeffs)

    br.addbounds(1, ['L'], [indz[i_viol]],  [brpt[i_viol]])
    br.addbounds(1, ['U'], [indiz[i_viol]], [1.0 / brpt[i_viol]])

    if zub[i_viol] < 1e10:
        coeffs, rhs = get_secant_cut(brpt[i_viol], zub[i_viol])
        br.addrows(1, ['L'], [rhs], [0, 2], indices, coeffs)

    return br



def xpress_cut_nonconvex(prob, m, sol, lb, ub):

    """This callback adds a secant cut and one or more OA cuts for the
    constraint

    inv_z = 1/z

    for the available bound values of z and inv_z. After tightening
    the bounds on z and inv_z when necessary, check if either an OA
    cut or a secant cut are violated, and add them. Do this for all
    districts.
    """

    DG = m._DG

    indz  = np.array([m.xmodel.getIndex(m._z[i])     for i in range(DG._k)])
    indiz = np.array([m.xmodel.getIndex(m._inv_z[i]) for i in range(DG._k)])

    zval = sol[indz]
    zlb  = lb[indz]
    zub  = ub[indz]

    invz = sol[indiz]
    izlb = lb[indiz]
    izub = ub[indiz]

    nSec = 0
    nOA = 0

    for i in range(DG._k):

        if zval[i]*invz[i] <= 1 + prob.controls.feastol or \
           zub[i] > 1.e19 or \
           (1/zlb[i] - 1/zub[i]) * zval[i] + (zub[i] - zlb[i]) * invz[i] <= \
           zub[i] / zlb[i] - zlb[i] / zub[i] + prob.controls.feastol:
            continue

        # A secant cut can be added, if violated. Renaming z as x and
        # inv_z as y,
        #
        # y - y1 <= (y2 - y1) / (x2 - x1) (x - x1)
        #
        # where (x1,x2) is (zlb,zub), (y1,y2) is (1/zlb, 1/zub)
        #
        # (x2 - x1) y - (x2 - x1) y1 <= (y2 - y1) x - (y2 - y1) x1
        # (y1 - y2) x + (x2 - x1) y  <= (y1 - y2) x1 + (x2 - x1) y1  = x2 y1 - x1 y2
        #
        # or
        #
        # (1/zlb - 1/zub) z + (zub - zlb) inv_z <= zub / zlb - zlb / zub

        indices = [indz[i], indiz[i]]
        coeffs, rhs = get_secant_cut (zlb[i], zub[i])

        mcolsp, dvalp = [], []
        drhsp, status = prob.presolverow('L', indices, coeffs, rhs, prob.attributes.cols,
                                         mcolsp, dvalp)

        if status >= 0:
            nSec += 1
            prob.addcuts([1], ['L'], [drhsp], [0, len(mcolsp)], mcolsp, dvalp)
        else:
            print("unpresolvable OA row")

    return 0


    # Future work:

    # A valid cut for inv_z = 1/z is the secant one, given we have
    # tight-enough bounds on z. We know >= 2 pi, we we must find lower
    # and upper bounds on z.

    # Because I couldn't find a way to get lb/ub of variables within a
    # callback with Gurobi (sigh), all we can use is the objective
    # function and its lower (primal) and upper (dual) bounds.

    #obj_best  = m.cbGet(GRB.Callback.MIPNODE_OBJBST)
    #obj_bound = m.cbGet(GRB.Callback.MIPNODE_OBJBND)

    # sum_i 1/z_i in [LB,UB] means what for each z_i? Obviously 1/z_i
    # <= UB, or z_i >= 1/UB, but this is a loose lower bound since 2pi
    # is certainly much higher already.

    # For z_j we have sum_{i!=j} 1/z_i + 1/z_j >= LB means 1/z_j >= LB
    # - sum_{i!=j} 1/LB(z_i) = LB - (k-1)/(2pi) =: alpha, so that z_j
    # <= 1/alpha if alpha > 0 and it should be 1/alpha > 2pi.

    # Given this bound, the secant cut for inv_z = 1/z for z in [2pi,
    # 1/alpha] = [zlb,zub] is the area through the two points
    # (zlb,1/zlb) and (zub,1/zub), i.e.
    #
    # inv_z - 1/zlb <= (1/zub - 1/zlb)/(zub - zlb) (z - zlb),
    #
    # which simplifies to
    #
    # inv_z - 1/zlb <= -1/(zlb*zub)*(z - zlb)
    #
    # zlb*zub * inv_z + z <= zlb + zub
    #
    # or, to recap,
    #
    # 2pi/alpha * inv_z + z <= 2pi + 1/alpha


    # for i,zv in enumerate(zval):
    #     izv = izval[i]

    #     alpha = obj_best - (m._DG._k - 1) / (2*math.pi)

    #     assert alpha < 1/(2*math.pi)

    #     izcoef = 2*math.pi/alpha
    #     rhs = 2*math.pi + 1/alpha

    #     if izcoef * izv + zv > rhs + eps:
    #         m.cbCut(izcoef * vars[m._inv_z[i]] + vars[m._z[i]] <= rhs)

    # return


def xpress_chksol_cb(prob, m, soltype, cutoff):
    """This callback checks if

    inv_z = 1/z

    for all districts, which is needed for feasibility of a solution
    """

    DG = m._DG

    try:
        x = []
        prob.getlpsol(x)  # Postsolved solution
    except:
        return (1,0)

    indz  = [m.xmodel.getIndex(m._z[i])     for i in range(DG._k)]
    indiz = [m.xmodel.getIndex(m._inv_z[i]) for i in range(DG._k)]

    zval = [x[ind] for ind in indz]
    invz = [x[ind] for ind in indiz]

    maxviol = max([abs(zval[i] * invz[i] - 1) for i in range(DG._k)])

    for i in range(DG._k):
        x[indiz[i]] = 1/x[indz[i]]

    objval = m._obj_coef * sum(1.0 / zval[i] for i in range(DG._k))
    if objval > m._xpress_bestobj + 1e-6:
        m._xpress_bestobj = objval
        prob.addmipsol(x)

    if maxviol > prob.controls.feastol:
        if soltype == 0:
            #print(f"node solution: {m._obj_coef * sum(1.0/zval[i] for i in range(DG._k))}, reported: {cutoff}")
            return (0, m._obj_coef * sum(1.0/zval[i] for i in range(DG._k)))
        return (1, 0)
    else:
        #print(f"computed: {m._obj_coef * sum(1.0/zval[i] for i in range(DG._k))}, reported: {cutoff}")
        return (0, cutoff) # m._obj_coef * sum(1.0/zval[i] for i in range(DG._k)))
