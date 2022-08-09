import xprgrb as gp
from xprgrb import GRB

def sort_by_second(val):
    return val[1]

def find_ordering(DG, B):
    V_B_with_population = [ ( i, DG.nodes[i]['TOTPOP'] ) for i in DG.nodes if i not in B ]
    V_B_with_population.sort( key=sort_by_second, reverse=True )
    return [ v for (v,p) in V_B_with_population ] + B

def solve_maxB_problem(DG):
    print("Solving the max B problem (as MIP) for use in the vertex ordering...") 
    m = gp.Model()
    m.Params.LogToConsole = 0 # keep log to a minimum
    q = DG._k
    
    # X[i,j]=1 if vertex i is assigned to bin j
    X = m.addVars(DG.nodes, q, name='X', vtype=GRB.BINARY)
    
    # B[i]=1 if vertex i is selected in set B
    B = m.addVars(DG.nodes, name='B', vtype=GRB.BINARY)
   
    # assignment constraints            
    m.addConstrs( gp.quicksum( X[i,j] for j in range(q) ) == B[i] for i in DG.nodes )
                
    # bin population should be less than L
    m.addConstrs( gp.quicksum( DG.nodes[i]['TOTPOP'] * X[i,j] for i in DG.nodes) <= DG._L-1 for j in range(q) )
    
    # bins shouldn't touch each other
    m.addConstrs( X[u,j] + B[v] <= 1 + X[v,j] for u,v in DG.edges for j in range(q) )
    
    # objective is to maximize size of set B
    m.setObjective( gp.quicksum( B ), GRB.MAXIMIZE )
    
    m.Params.MIPFocus = 1 # turn on MIPFocus
    m.Params.TimeLimit = 60 # 60-second time limit
    m.optimize()
    
    if m.status in { GRB.OPTIMAL, GRB.TIME_LIMIT }:
        B_sol = [ i for i in DG.nodes if gp.getsol(B[i],m) > 0.5 ]
    else:
        B_sol = list()
    #print("Size of B =",len(B_sol))
    return (B_sol,m.runtime)
    
 
