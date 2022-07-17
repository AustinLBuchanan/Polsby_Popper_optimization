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


class Params:
    def __init__(self):
        self.FeasibilityTol = None
        self.IntFeasTol = None
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
            self.xmodel.addConstr(constraint) 
            return constraint


    def addConstrs(self, *constraints):
        """
        Add multiple constraints
        """
        if solver == 'gurobi':
            return self.gmodel.addConstrs(*constraints)
        elif solver == 'xpress':
            self.xmodel.addConstr(*constraints)
            return constraint


    def addVar(self, **params):
        if solver == 'gurobi':
            return self.gmodel.addVar(**params)
        elif solver == 'xpress':
            x = xp.var(**params)
            self.xmodel.addVar(x)
            return x


    def addVars(self, *indices, **params):
        if solver == 'gurobi':
            return self.gmodel.addVars(*indices, **params)
        elif solver == 'xpress':
            x = xp.vars(*indices, **params)
            self.xmodel.addVar(x)
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
            if self.Params.Method          is not None: self.gmodel.Params.Method          = self.Params.Method
            if self.Params.MIPFocus        is not None: self.gmodel.Params.MIPFocus        = self.Params.MIPFocus
            if self.Params.MIPGap          is not None: self.gmodel.Params.MIPGap          = self.Params.MIPGap
            if self.Params.NonConvex       is not None: self.gmodel.Params.NonConvex       = self.Params.NonConvex
            if self.Params.OutputFlag      is not None: self.gmodel.Params.OutputFlag      = self.Params.OutputFlag
            if self.Params.TimeLimit       is not None: self.gmodel.Params.TimeLimit       = self.Params.TimeLimit

            self.gmodel.optimize(callback)

        elif solver == 'xpress':

            if self.Params.FeasibilityTol  is not None: self.gmodel.controls.FeasibilityTol  = self.Params.FeasibilityTol
            if self.Params.IntFeasTol      is not None: self.gmodel.controls.IntFeasTol      = self.Params.IntFeasTol
            if self.Params.LazyConstraints is not None: self.gmodel.controls.LazyConstraints = self.Params.LazyConstraints
            if self.Params.Method          is not None: self.gmodel.controls.Method          = self.Params.Method
            if self.Params.MIPFocus        is not None: self.gmodel.controls.MIPFocus        = self.Params.MIPFocus
            if self.Params.MIPGap          is not None: self.gmodel.controls.MIPGap          = self.Params.MIPGap
            if self.Params.NonConvex       is not None: self.gmodel.controls.NonConvex       = self.Params.NonConvex
            if self.Params.OutputFlag      is not None: self.gmodel.controls.OutputFlag      = self.Params.OutputFlag
            if self.Params.TimeLimit       is not None: self.gmodel.controls.TimeLimit       = self.Params.TimeLimit

            self.xmodel.solve()


    def update(self):
        if solver == 'gurobi':
            self.gmodel.update()


    def __getattr__ (self, name):

        if solver == 'gurobi':
            return getattr(self.gmodel, name)
        else:
            pass


    def __setattr__ (self, name, value):

        if name in ['gmodel', 'xmodel', 'Params', '_callback']:
            return object.__setattr__(self, name, value)
        else:
            return object.__setattr__(self.gmodel, name, value)


GRB = gurobipy.GRB
