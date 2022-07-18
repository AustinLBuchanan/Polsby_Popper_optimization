"""Wrapper class for Gurobi/Xpress problems. Allows for using Gurobi
notation while also creating Xpress problems if specified through
setSolver()
"""

import gurobipy
import xpress as xp

solver = 'gurobi'


def setSolver(s):
    """
    Globally set the solver between Gurobi and Xpress
    """
    global solver
    solver = s


def quicksum(arg):
    if solver == 'gurobi':
        return gurobipy.quicksum(arg)
    else:
        return xp.Sum(arg)


def getsol(var, prob):
    if solver == 'gurobi':
        return var.x
    else:
        return prob.xmodel.getSolution(var)


def getLB(var, prob):
    if solver == 'gurobi':
        return var.LB
    else:
        bds = []
        prob.xmodel.getlb(bds, var, var)
        return bds[0]


def getUB(var, prob):
    if solver == 'gurobi':
        return var.UB
    else:
        bds = []
        prob.xmodel.getub(bds, var, var)
        return bds[0]


def setLB(var, prob, value):
    if solver == 'gurobi':
        var.LB = value
    else:
        prob.xmodel.chgbounds([var], ['L'], [value])


def setUB(var, prob, value):
    if solver == 'gurobi':
        var.UB = value
    else:
        prob.xmodel.chgbounds([var], ['U'], [value])


class Params:
    def __init__(self):
        self.FeasibilityTol = None
        self.IntFeasTol = None
        self.LogToConsole = None
        self.LazyConstraints = None
        self.Method = None
        self.MIPFocus = None
        self.MIPGap = None
        self.NonConvex = None
        self.OutputFlag = None
        self.TimeLimit = None


class Model:

    def __init__(self):
        """
        Define a model independently of the solver
        """

        self.gmodel = None
        self.xmodel = None
        self.Params = Params()
        self._callback = None

        if solver == 'gurobi':
            self.gmodel = gurobipy.Model()
        elif solver == 'xpress':
            self.xmodel = xp.problem()
        else:
            raise RuntimeError('Solver must be "gurobi" or "xpress"')


    def addcallback(callback, type):
        pass


    def addConstr(self, constraint):
        """
        Add a single constraint
        """
        if solver == 'gurobi':
            return self.gmodel.addConstr(constraint)
        elif solver == 'xpress':
            self.xmodel.addConstraint(constraint)
            return constraint


    def addConstrs(self, *constraints):
        """
        Add multiple constraints
        """
        if solver == 'gurobi':
            return self.gmodel.addConstrs(*constraints)
        elif solver == 'xpress':
            self.xmodel.addConstraint(*constraints)
            return constraints


    def addVar(self, **params):
        if solver == 'gurobi':
            return self.gmodel.addVar(**params)
        elif solver == 'xpress':
            x = xp.var(**params)
            self.xmodel.addVariable(x)
            return x


    def addVars(self, *indices, **params):
        if solver == 'gurobi':
            return self.gmodel.addVars(*indices, **params)
        elif solver == 'xpress':
            if 'vtype' in params.keys():
                params['vartype'] = {'B': xp.binary, 'I': xp.integer, 'C': xp.continuous} [params['vtype']]
                del params['vtype']
            args = []
            for i,_ in enumerate(indices):
                if isinstance(indices[i], int):
                    args.append(indices[i])
                else:
                    args.append(list(indices[i]))
            x = xp.vars(*args, **params)
            self.xmodel.addVariable(x)
            return x


    def setObjective(self, *indices, **params):
        if solver == 'gurobi':
            return self.gmodel.setObjective(*indices, **params)
        elif solver == 'xpress':
            return self.xmodel.setObjective(*indices, **params)


    def optimize(self, callback=None):
        if solver == 'gurobi':

            if self.Params.FeasibilityTol  is not None: self.gmodel.Params.FeasibilityTol  = self.Params.FeasibilityTol
            if self.Params.IntFeasTol      is not None: self.gmodel.Params.IntFeasTol      = self.Params.IntFeasTol
            if self.Params.LazyConstraints is not None: self.gmodel.Params.LazyConstraints = self.Params.LazyConstraints
            if self.Params.LogToConsole    is not None: self.gmodel.Params.LogToConsole    = self.Params.LogToConsole
            if self.Params.Method          is not None: self.gmodel.Params.Method          = self.Params.Method
            if self.Params.MIPFocus        is not None: self.gmodel.Params.MIPFocus        = self.Params.MIPFocus
            if self.Params.MIPGap          is not None: self.gmodel.Params.MIPGap          = self.Params.MIPGap
            if self.Params.NonConvex       is not None: self.gmodel.Params.NonConvex       = self.Params.NonConvex
            if self.Params.OutputFlag      is not None: self.gmodel.Params.OutputFlag      = self.Params.OutputFlag
            if self.Params.TimeLimit       is not None: self.gmodel.Params.TimeLimit       = self.Params.TimeLimit

            self.gmodel.optimize(callback)

        elif solver == 'xpress':

            if self.Params.FeasibilityTol  is not None: self.xmodel.controls.feastol         = self.Params.FeasibilityTol
            if self.Params.IntFeasTol      is not None: self.xmodel.controls.miptol          = self.Params.IntFeasTol
            #if self.Params.LazyConstraints is not None: self.xmodel.controls.LazyConstraints = self.Params.LazyConstraints
            #if self.Params.LogToConsole    is not None: self.xmodel.controls.LogToConsole    = self.Params.LogToConsole
            #if self.Params.Method          is not None: self.xmodel.controls.Method          = self.Params.Method
            #if self.Params.MIPFocus        is not None: self.xmodel.controls.MIPFocus        = self.Params.MIPFocus
            if self.Params.MIPGap          is not None: self.xmodel.controls.mipgap          = self.Params.MIPGap
            #if self.Params.NonConvex       is not None: self.xmodel.controls.NonConvex       = self.Params.NonConvex
            #if self.Params.OutputFlag      is not None: self.xmodel.controls.OutputFlag      = self.Params.OutputFlag
            if self.Params.TimeLimit       is not None: self.xmodel.controls.maxtime       = self.Params.TimeLimit

            self.xmodel.solve()


    def update(self):
        if solver == 'gurobi':
            self.gmodel.update()


    def __getattr__ (self, name):

        if solver == 'gurobi':
            return getattr(self.gmodel, name)
        else:
            if name == 'status':
                return {xp.mip_infeas: gurobipy.GRB.INFEASIBLE,
                        xp.mip_optimal: gurobipy.GRB.OPTIMAL,
                        xp.mip_solution: gurobipy.GRB.TIME_LIMIT}[self.xmodel.getProbStatus()]
            elif name == 'runtime':
                return self.xmodel.attributes.time

    def __setattr__ (self, name, value):

        if name in ['gmodel', 'xmodel', 'Params', '_callback']:
            return object.__setattr__(self, name, value)
        else:
            if solver == 'gurobi':
                return object.__setattr__(self.gmodel, name, value)
            else:
                return object.__setattr__(self, name, value)


GRB = gurobipy.GRB
