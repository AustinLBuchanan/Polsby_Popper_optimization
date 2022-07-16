import gurobipy as gp
import xpress as xp

solver = None


def setSolver(s):
    solver = s


def quicksum(arg):
    if solver == 'gurobi':
        return gp.quicksum(arg)
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

        if solver == 'gurobi':
            self.gmodel = gp.Model()
        elif solver == 'xpress':
            self.xmodel = xp.problem()
        else:
            raise RuntimeError('Solver must be "gurobi" or "xpress"')


    def addConstr(self, constraint):
        """
        Add a single constraint
        """
        if solver == 'gurobi':
            return self.gmodel.AddConstr(constraint)
        elif solver == 'xpress':
            self.xmodel.addConstr(constraint) 
            return constraint


    def addConstrs(self, *constraints)
        """
        Add multiple constraints
        """
        if solver == 'gurobi':
            return self.gmodel.AddConstrs(*constraint)
        elif solver == 'xpress':
            self.xmodel.addConstr(*constraint) 
            return constraint


    def addVar(self, **params):
        if solver == 'gurobi':
            return self.gmodel.AddVar(**params)
        elif solver == 'xpress':
            x = xp.var(**params)
            self.xmodel.addVar(x)
            return x


    def addVars(self, *indices, **params):
        if solver == 'gurobi':
            return self.gmodel.AddVars(*indices, **params)
        elif solver == 'xpress':
            x = xp.vars(*indices, **params)
            self.xmodel.addVar(x)
            return x


    def optimize(self):
        if solver == 'gurobi':
            self.gmodel.optimize()
        elif solver == 'xpress':
            self.xmodel.solve()
        

GRB = gp.GRB
