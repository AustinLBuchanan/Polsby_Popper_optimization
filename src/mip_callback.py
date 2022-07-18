import xprgrb as gp
from xprgrb import GRB, solver
import mip_contiguity
import math


def gen_callback(m, where):

    assert solver == 'xpress'

    raise RuntimeError("Should NOT get here, this callback is no use")

    if where != GRB.Callback.MIPNODE or \
       m.cbGet(GRB.Callback.MIPNODE_STATUS) != GRB.OPTIMAL:
        return

    # A valid cut for inv_z = 1/z is the secant one, given we have
    # tight-enough bounds on z. We know >= 2 pi, we we must find lower
    # and upper bounds on z.

    # Because I couldn't find a way to get lb/ub of variables within a
    # callback with Gurobi (sigh), all we can use is the objective
    # function and its lower (primal) and upper (dual) bounds.

    obj_best  = m.cbGet(GRB.Callback.MIPNODE_OBJBST)
    obj_bound = m.cbGet(GRB.Callback.MIPNODE_OBJBND)

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

    zval = m.cbGetNodeRel(m._z)
    izval = m.cbGetNodeRel(m._inv_z)

    eps = 1e-6

    vars = m.getVars()

    for i,zv in enumerate(zval):
        izv = izval[i]

        alpha = obj_best - (m._DG._k - 1) / (2*math.pi)

        assert alpha < 1/(2*math.pi)

        izcoef = 2*math.pi/alpha
        rhs = 2*math.pi + 1/alpha

        if izcoef * izv + zv > rhs + eps:
            m.cbCut(izcoef * vars[m._inv_z[i]] + vars[m._z[i]] <= rhs)

    return


def xpress_branch_cb(nodeprob, DG):
    """This callback adds a branch object on z or inv_z if the current
    relaxation cannot be tightened via secant cuts (it is always
    possible to tighten one via OA cuts"
    """
    pass


def xpress_cut_cb(nodeprob, DG):
    """This callback adds a secant cut and one or more OA cuts for the
    constraint

    inv_z = 1/z

    for the available bound values of z and inv_z. After tightening
    the bounds on z and inv_z when necessary, check if either an OA
    cut or a secant cut are violated, and add them. Do this for all
    districts.
    """

    pass


def xpress_chksol_cb(nodeprob, DG):
    """This callback checks if

    inv_z = 1/z

    for all districts, which is needed for feasibility of a solution
    """

    pass


# # Define callback functions: one for checking if a solution is
# # feasible (apart from linearity, all auxiliaries must be equal to
# # their respective product); one for adding new McCormick inequalities
# # for changed bounds; and finally, one for branching as we might have
# # to branch on continuous variables.

# def cbbranch(prob, aux, branch):
#     """Branch callback. Receives branch in input and, if it finds
#     continuous branches that are violated, adds them.
#     """

#     sol = []

#     if (prob.attributes.presolvestate & 128) == 0:
#         return branch

#     # Retrieve node solution
#     try:
#         prob.getlpsol(x=sol)
#     except:
#         return branch

#     lb, ub = getCBbounds(prob, len(sol))

#     assert(len(lb) == len(ub))
#     assert(len(sol) == len(lb))

#     x = prob.getVariable()  # presolved variables

#     rowmap = []
#     colmap = []

#     prob.getpresolvemap(rowmap, colmap)

#     invcolmap = [-1 for _ in lb]

#     for i, m in enumerate(colmap):
#         invcolmap[m] = i

#     # make sure all dimensions match

#     assert (len(lb) == len(ub))
#     assert (len(sol) == len(lb))
#     assert (len(invcolmap) == len(lb))

#     # Check if all auxiliaries are equal to their respective bilinear
#     # term. If so, we have a feasible solution

#     sol = np.array(sol)

#     discr = sol[Aux_ind] - sol[Aux_i] * sol[Aux_j]
#     discr[Aux_i == Aux_j] = np.maximum(0, discr[Aux_i == Aux_j])
#     maxdiscind = np.argmax(np.abs(discr))

#     if abs(discr[maxdiscind]) < eps:
#         return branch

#     i,j = Aux_i[maxdiscind], Aux_j[maxdiscind]

#     yind = prob.getIndex(aux[i, j])

#     if i == j:

#         # Test of violation is done on the original
#         # space. However, the problem variables are scrambled with invcolmap

#         if sol[i] > lb[i] + eps and \
#            sol[i] < ub[i] - eps and \
#            sol[yind] > sol[i]**2 + eps and \
#            sol[yind] - lb[i]**2 <= (ub[i] + lb[i]) * (sol[i] - lb[i]) - eps:

#             # Can't separate, must branch. Otherwise OA or secant
#             # cut separated above should be enough

#             brvarind = invcolmap[i]
#             brpoint = sol[i]
#             brvar = x[brvarind]
#             brleft = brpoint
#             brright = brpoint

#             assert(brvarind >= 0)

#             if brvar.vartype in [xp.integer, xp.binary]:
#                 brleft = math.floor(brpoint + 1e-5)
#                 brright = math.ceil(brpoint - 1e-5)

#             b = xp.branchobj(prob, isoriginal=False)

#             b.addbranches(2)

#             addrowzip(prob, b, 0, 'L', brleft,  [i], [1])
#             addrowzip(prob, b, 1, 'G', brright, [i], [1])

#             # New variable bounds are not enough, add new McCormick
#             # inequalities for y = x**2: suppose x0,y0 are the current
#             # solution values for x,y, yp = x0**2 and xu,yu = xu**2 are their
#             # upper bound, and similar for lower bound. Then these two
#             # rows must be added, one for each branch:
#             #
#             # y - yp <= (yl-yp)/(xl-x0) * (x - x0)  <===>
#             # (yl-yp)/(xl-x0) * x - y >= (yl-yp)/(xl-x0) * x0 - yp
#             #
#             # y - yp <= (yu-yp)/(xu-x0) * (x - x0)  <===>
#             # (yu-yp)/(xu-x0) * x - y >= (yu-yp)/(xu-x0) * x0 - yp
#             #
#             # Obviously do this only for finite bounds

#             ypl = brleft**2
#             ypr = brright**2

#             if lb[i] > -1e7 and sol[i] > lb[i] + eps:

#                 yl = lb[i]**2
#                 coeff = (yl - ypl) / (lb[i] - sol[i])

#                 if coeff != 0:
#                     addrowzip(prob, b, 0, 'G', coeff*sol[i] - ypl,
#                               [i, yind], [coeff, -1])

#             if ub[i] < 1e7 and sol[i] < ub[i] - eps:

#                 yu = ub[i]**2
#                 coeff = (yu - ypr) / (ub[i] - sol[i])

#                 if coeff != 0:
#                     addrowzip(prob, b, 1, 'G', coeff*sol[i] - ypr,
#                               [i, yind], [coeff, -1])

#             return b

#     else:

#         lbi0, ubi0 = lb[i], ub[i]
#         lbi1, ubi1 = lb[i], ub[i]

#         lbj0, ubj0 = lb[j], ub[j]
#         lbj1, ubj1 = lb[j], ub[j]

#         # No cut violated, must branch
#         if min(sol[i] - lb[i], ub[i] - sol[i]) / (1 + ub[i] - lb[i]) > \
#            min(sol[j] - lb[j], ub[j] - sol[j]) / (1 + ub[j] - lb[j]):
#             lbi1 = sol[i]
#             ubi0 = sol[i]
#             brvar = i
#         else:
#             lbj1 = sol[j]
#             ubj0 = sol[j]
#             brvar = j

#         alpha = 0.2

#         brvarind = invcolmap[brvar]
#         brpoint = sol[brvar]
#         brleft = brpoint
#         brright = brpoint

#         if x[brvarind].vartype in [xp.integer, xp.binary]:
#             brleft = math.floor(brpoint + 1e-5)
#             brright = math.ceil(brpoint - 1e-5)

#         b = xp.branchobj(prob, isoriginal=False)

#         b.addbranches(2)

#         addrowzip(prob, b, 0, 'L', brleft,  [brvar], [1])
#         addrowzip(prob, b, 1, 'G', brright, [brvar], [1])

#         # As for the i==j case, the variable branch is
#         # insufficient, so add updated McCormick inequalities.
#         # There are two McCormick inequalities per changed bound:
#         #
#         # y >= lb[j] * x[i] + lb[i] * x[j] - lb[j] * lb[i] ---> add to branch 1
#         # y >= ub[j] * x[i] + ub[i] * x[j] - ub[j] * ub[i] ---> add to branch 0
#         # y <= lb[j] * x[i] + ub[i] * x[j] - lb[j] * ub[i] ---> add to branch 1 if x[brvarind] == j, 0 if x[brvarind] == i
#         # y <= ub[j] * x[i] + lb[i] * x[j] - ub[j] * lb[i] ---> add to branch 1 if x[brvarind] == i, 0 if x[brvarind] == j

#         addrowzip(prob, b, 0, 'G', - ubi0 * ubj0, [yind, i, j], [1, -ubj0, -ubi0])
#         addrowzip(prob, b, 1, 'G', - lbi1 * lbj1, [yind, i, j], [1, -lbj1, -lbi1])

#         if brvarind == i:
#             addrowzip(prob, b, 0, 'L', - lbj0 * ubi0, [yind, i, j], [1, -lbj0, -ubi0])
#             addrowzip(prob, b, 1, 'L', - ubj1 * lbi1, [yind, i, j], [1, -ubj1, -lbi1])
#         else:
#             addrowzip(prob, b, 0, 'L', - ubj0 * lbi0, [yind, i, j], [1, -ubj0, -lbi0])
#             addrowzip(prob, b, 1, 'L', - lbj1 * ubi1, [yind, i, j], [1, -lbj1, -ubi1])
#         return b

#     # If no branching rule was found, return none
#     return branch
