import networkx as nx
from district import build_single_district_mip
from coarsen import subgraph, graph_coarsen_by_county

# Repeatedly carve majority-minority districts from the graph, as able.
#
def carving_heuristic(G, minority, complement_contiguity='cut', complement_compactness=False, deviation_penalty=0.0, 
                      mvap_lower=0.5, mvap_upper=1.0, mvap_excess_penalty=0.001, polsby_popper_cutoff=0.125, verbose=False):
    
    districts = list()
    unassigned_nodes = list(G.nodes)
    
    # Coarsen all but 0, 1, or 2 counties
    for num_split in range(3):
        print(f"Seeking districts that split <= {num_split} counties.")

        for split_fips in get_ordering(G, num_split):
            print("Allowing to split:",split_fips)

            while len(districts) < G._k:
                
                # create coarsened graph instance
                S = subgraph(G, unassigned_nodes)
                H, coarsen_partition = graph_coarsen_by_county(S, granular_fips=split_fips)
                DH = nx.DiGraph(H)
                DH._L = G._L
                DH._U = G._U
                DH._CL = G._L * ( G._k - len(districts) - 1)
                DH._CU = G._U * ( G._k - len(districts) - 1)
                
                # build and solve MISOCP
                my_penalty = 1 / ( G._k - len(districts) - 1 ) if complement_compactness else 0
                m = build_single_district_mip(DH, contiguity='cut', complement_contiguity='cut',
                                              deviation_penalty=deviation_penalty, minority=minority, 
                                              mvap_lower=mvap_lower, mvap_upper=mvap_upper, 
                                              mvap_excess_penalty=mvap_excess_penalty,
                                              split_counties_limit=num_split,
                                              complement_compactness_penalty=my_penalty, verbose=verbose)

                m._z.UB = 1 / polsby_popper_cutoff
                m.Params.TimeLimit = 3600
                m.optimize(m._callback)
                if m.solCount == 0:
                    break

                print(f"District {1+len(districts)} with {minority} VAP {round(100*m._mvap.x/m._vap.x,2)}% and Polsby-Popper score {round(1/m._z.x,2)}")

                district = list()
                for i in DH.nodes:
                    if m._x[i].x > 0.5:
                        district += coarsen_partition[i]
                districts.append(district)
                unassigned_nodes = [ i for i in unassigned_nodes if i not in district ]
                assert len(unassigned_nodes)<=1 or nx.is_connected(G.subgraph(unassigned_nodes)), "Remaining subgraph is disconnected."
             
    return districts


# same idea, but without county coarsening/decomposition
#
def carving_heuristic_simple(G, minority, complement_contiguity='cut', complement_compactness=False, deviation_penalty=0.0, 
                             mvap_lower=0.5, mvap_upper=1.0, mvap_excess_penalty=0.001, polsby_popper_cutoff=0.125, verbose=False,
                             split_counties_limit=2):
    
    districts = list()
    unassigned_nodes = list(G.nodes)
    print(f"Seeking districts that split <= {split_counties_limit} counties.")

    while len(districts) < G._k:

        # create coarsened graph instance
        H = subgraph(G, unassigned_nodes)
        DH = nx.DiGraph(H)
        DH._L = G._L
        DH._U = G._U
        DH._CL = G._L * ( G._k - len(districts) - 1)
        DH._CU = G._U * ( G._k - len(districts) - 1)

        # build and solve MISOCP
        my_penalty = 1 / ( G._k - len(districts) - 1 ) if complement_compactness else 0
        m = build_single_district_mip(DH, contiguity='cut', complement_contiguity='cut',
                                      deviation_penalty=deviation_penalty, minority=minority, 
                                      mvap_lower=mvap_lower, mvap_upper=mvap_upper, 
                                      mvap_excess_penalty=mvap_excess_penalty,
                                      split_counties_limit=split_counties_limit,
                                      complement_compactness_penalty=my_penalty, verbose=verbose)

        m._z.UB = 1 / polsby_popper_cutoff
        m.Params.TimeLimit = 3600
        m.optimize(m._callback)
        if m.solCount == 0:
            break

        print(f"District {1+len(districts)} with {minority} VAP {round(100*m._mvap.x/m._vap.x,2)}% and Polsby-Popper score {round(1/m._z.x,2)}")

        district = [ i for i in DH.nodes if m._x[i].x > 0.5 ]
        districts.append(district)
        unassigned_nodes = [ i for i in unassigned_nodes if i not in district ]
        assert len(unassigned_nodes)<=1 or nx.is_connected(G.subgraph(unassigned_nodes)), "Remaining subgraph is disconnected."
             
    return districts


# two helper functions for get_ordering()
def sort_by_second(val):
    return val[1]

def fips_edge(G, i, j):
    fi = G.nodes[i]['GEOID20'][0:5]
    fj = G.nodes[j]['GEOID20'][0:5]
    return ( min(fi, fj), max(fi, fj) )

# by default, sort counties (or adjacent county pairs) from smallest pop to largest pop
def get_ordering(G, num):
    assert num <= 2, "get_ordering(G, num) is only defined for num <= 2"
    if num==0:
        return list()
    
    fips_pop = { G.nodes[i]['GEOID20'][0:5] : 0 for i in G.nodes }
    for i in G.nodes:
        f = G.nodes[i]['GEOID20'][0:5]
        fips_pop[f] += G.nodes[i]['TOTPOP']
            
    if num==1:
        fips_pop_combos = [ ( [f], fips_pop[f] ) for f in fips_pop.keys() ]
    else: # num==2
        fips_edges = { fips_edge(G,i,j) for i,j in G.edges if G.nodes[i]['GEOID20'][0:5] != G.nodes[j]['GEOID20'][0:5] }
        fips_pop_combos = [ ( [f1,f2], fips_pop[f1]+fips_pop[f2] ) for (f1,f2) in fips_edges ]
        
    fips_pop_combos.sort( key=sort_by_second )
    return [ f for (f,pop) in fips_pop_combos ]    
