import xprgrb as gp
from xprgrb import GRB
import mip_contiguity
import math

def gen_lcut_callback(m, where):

    if where == GRB.Callback.MIPSOL:
        return mip_contiguity.lcut_callback(m, where)
    else:
        return gen_callback(m, where)


def gen_callback(m, where):

    return

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
