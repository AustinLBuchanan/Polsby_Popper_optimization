import networkx as nx
import xprgrb as gp
from xprgrb import GRB
import mip, mip_objective, mip_contiguity, mip_fixing

# In the spirit of "ILP-based local search for graph partitioning"
#    by A Henzinger, A Noe, C Schulz - Journal of Experimental Algorithmics, 2020
#    https://scholar.google.com/scholar?cluster=9173497710013156715&hl=en&as_sdt=0,37

def local_search(m, DG, labeling, radius=1):
    
    print("Applying MIP-based local search to improve the MIP warm start...")
    
    # gurobi apparently clears the callback at m.optimize(), 
    #   so we need to add it fresh each time (?)
    my_callback = m._callback 
    
    # get objective value of initial labeling
    m.Params.OutputFlag = 0
    m.Params.TimeLimit = 60 
    set_x_ub_wrt_labeling(m, DG, labeling, 0)
    mip.inject_warm_start(m, DG, labeling)
    m.optimize(my_callback)
    grb_time = m.runtime
    
    # LOCAL SEARCH
    max_iterations = 10
    print("iter \t\t obj \t\t time")
    print(0,'\t','{0:.8f}'.format(m.objVal),'\t','{0:.2f}'.format(m.runtime))
    
    for iteration in range(1,max_iterations+1):
        
        old_obj = m.objVal
        
        set_x_ub_wrt_labeling(m, DG, labeling, radius)
        m.optimize(my_callback)
        
        grb_time += m.runtime
        new_obj = m.objVal
        print(iteration,'\t','{0:.8f}'.format(new_obj),'\t','{0:.2f}'.format(m.runtime))
        
        if new_obj == old_obj:
            break
        else:
            labeling = { i : j for i in DG.nodes for j in range(DG._k) if m._x[i,j].x > 0.5 }
    
    # reset model
    m._numCallbacks = 0
    m._numLazyCuts = 0
    m._callback = my_callback
    set_x_ub(m, DG, ub=1)
    m.Params.OutputFlag = 1
    m.update()
    
    # return the final labeling
    return (labeling, new_obj, grb_time)


# set all x[i,j] upper bounds to ub 
#
def set_x_ub(m, DG, ub):
    for i in DG.nodes:
        for j in range(DG._k):
            m._x[i,j].ub = ub
    m.update()
    return


# set x[i,j].ub=1 if and only if there is a
#  vertex v with dist(i,j)<=radius and labeling[v]=j
#
def set_x_ub_wrt_labeling(m, DG, labeling, radius=1):
    set_x_ub(m, DG, ub=0)
    for i in DG.nodes:
        dist = nx.single_source_shortest_path_length(DG, source=i, cutoff=radius)
        for v in dist.keys():
            j = labeling[v]
            m._x[i,j].ub = 1
    m.update()
    return
