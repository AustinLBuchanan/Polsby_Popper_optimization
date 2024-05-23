import json
from networkx.readwrite import json_graph

def read_graph_from_json(json_file, update_population=True, rescale_distance=True):
    with open(json_file) as f:
        data = json.load(f)
    G = json_graph.adjacency_graph(data) 

    if update_population:
        # total population
        for i in G.nodes:
            G.nodes[i]['TOTPOP'] = G.nodes[i]['P0010001'] 
    
    # change distance units from meters to 100km for better numerics in gurobi
    #  note: 100km roughly corresponds to one degree of latitude
    if rescale_distance:
        ht = 100000  
        for i in G.nodes:
            G.nodes[i]['area'] /= ( ht * ht )
            if G.nodes[i]['boundary_node']:
                G.nodes[i]['boundary_perim'] /= ht
        for i,j in G.edges:
            G.edges[i,j]['shared_perim'] /= ht
        
    return G