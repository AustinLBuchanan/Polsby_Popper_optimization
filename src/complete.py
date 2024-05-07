import math
import networkx as nx
from district import build_single_district_mip

# Repeatedly bipartition the state into (nearly) equal size halves.
# Continue until each part is a district.
# In each iteration, try to split at most one county.
# For speed, use cut edges objective. This also promotes contiguity in lower-level subproblems.
# We fix the gingles_districts and apply the algorithm to the rest of the state.
# 
def recursive_bipartition_heuristic(G, gingles_districts=list()):
    
    assigned_nodes = [ i for j in range(len(gingles_districts)) for i in gingles_districts[j] ]
    unassigned_nodes = [ i for i in G.nodes if i not in assigned_nodes ]
    
    clusters = [ unassigned_nodes ]
    sizes = [ G._k - len(gingles_districts) ]
    districts = gingles_districts.copy()
    num_multidistricts = 0
    DG = nx.DiGraph(G)

    while len(sizes) > 0:

        # pick a cluster
        cluster = clusters.pop()
        size = sizes.pop()

        size1 = math.floor( size / 2 )
        size2 = math.ceil( size / 2 )

        DH = DG.subgraph(cluster)
        DH._L = size1 * G._L
        DH._U = size1 * G._U
        DH._CL = size2 * G._L
        DH._CU = size2 * G._U

        print("Using one split county, attempting to bipartition cluster into sizes:",size1,size2)

        # first, try to bipartition using 1 split county
        m = build_single_district_mip(DH, objective='cut_edges', split_counties_limit=1, deviation_penalty=0.001, verbose=False,
                                      contiguity='cut', complement_contiguity='cut', complement_balance=True)

        m.Params.TimeLimit = 3600
        m.optimize(m._callback)

        if m.solCount == 0:

            # otherwise, allow any number of split counties
            print("Without limiting splits, attempting to bipartition cluster")
            m = build_single_district_mip(DH, objective='cut_edges', deviation_penalty=0.001, verbose=False,
                                      contiguity='cut', complement_contiguity='cut', complement_balance=True)

            m.Params.TimeLimit = 3600
            m.optimize(m._callback)

            if m.solCount == 0:
                print("Unable to bipartition. Keeping as multidistrict of size =",size)
                districts.append(cluster)
                num_multidistricts += 1
                continue

        cluster1 = [ i for i in DH.nodes if m._x[i].x > 0.5 ]
        cluster2 = [ i for i in DH.nodes if m._x[i].x < 0.5 ]

        if size1 == 1:
            districts.append(cluster1)
        else:
            clusters.append(cluster1)
            sizes.append(size1)

        if size2 == 1:
            districts.append(cluster2)
        else:
            clusters.append(cluster2)
            sizes.append(size2)
            
    print(f"After recursive partitioning, we have {len(districts)} districts.")
    if num_multidistricts > 0:
        print(f"This includes {num_multidistricts} multidistricts.")
    return districts
