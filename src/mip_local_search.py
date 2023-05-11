import networkx as nx
import xprgrb as gp
from xprgrb import GRB
import mip, mip_objective, mip_contiguity, mip_fixing
from tract_approximation import distance_to_vertex_set

# In the spirit of "ILP-based local search for graph partitioning"
#    by A Henzinger, A Noe, C Schulz - Journal of Experimental Algorithmics, 2020
#    https://scholar.google.com/scholar?cluster=9173497710013156715&hl=en&as_sdt=0,37

def local_search(m, DG, labeling, radius=1, max_iterations=10, preserve_splits=False):
    
    print(f"Applying MIP-based local search to improve the MIP warm start (radius: {radius})...")
    
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

    if m.status not in [GRB.TIME_LIMIT, GRB.OPTIMAL]:
        return (None, -1e20, grb_time)

    new_obj = m.objVal
    
    # LOCAL SEARCH
    print("iter \t\t obj \t\t time")
    print(0,'\t','{0:.8f}'.format(m.objVal),'\t','{0:.2f}'.format(m.runtime))
    
    for iteration in range(1,max_iterations+1):
        
        old_obj = m.objVal
        
        set_x_ub_wrt_labeling(m, DG, labeling, radius, preserve_splits)
        m.optimize(my_callback)
        if m.status not in [GRB.TIME_LIMIT, GRB.OPTIMAL]:
            break

        grb_time += m.runtime
        new_obj = m.objVal
        print(iteration,'\t','{0:.8f}'.format(new_obj),'\t','{0:.2f}'.format(m.runtime))
        
        if new_obj == old_obj:
            break
        else:
            sol = gp.getsol(None, m, None)
            labeling = { i : j for i in DG.nodes for j in range(DG._k) if gp.getsol(m._x[i,j], m, sol) > 0.5 }
    
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
            gp.setUB(m._x[i,j], m, ub)
            #m._x[i,j].ub = ub
    m.update()
    return


# set x[i,j].ub=1 if and only if there is a
#  vertex v with dist(i,j)<=radius and labeling[v]=j 
#
def set_x_ub_wrt_labeling(m, DG, labeling, radius=1, preserve_splits=False):
    set_x_ub(m, DG, ub=0)
    
    for j in range(DG._k):
        district = [ i for i in DG.nodes if labeling[i] == j ]
        dist = distance_to_vertex_set(DG, district, cutoff=radius)
        for i in dist.keys():
            gp.setUB(m._x[i,j], m, 1)
            #m._x[i,j].ub = 1
            
    if preserve_splits:
        counties = { DG.nodes[i]['GEOID20'][0:5] for i in DG.nodes }
        support = { c : list() for c in counties }
        
        for i in labeling.keys():
            j = labeling[i]
            c = DG.nodes[i]['GEOID20'][0:5]
            if j not in support[c]:
                support[c].append(j)
        
        nonsupport = { c : set(range(DG._k)) - set(support[c]) for c in counties }
        
        for i in DG.nodes:
            c = DG.nodes[i]['GEOID20'][0:5]
            for j in nonsupport[c]:
                gp.setUB(m._x[i,j], m, 0)
                
            
    m.update()
    return
