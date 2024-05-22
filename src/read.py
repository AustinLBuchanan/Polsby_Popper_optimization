import json
from networkx.readwrite import json_graph

def read_graph_from_json(json_file):
    with open(json_file) as f:
        data = json.load(f)
    G = json_graph.adjacency_graph(data) 

    # total population
    for i in G.nodes:
        G.nodes[i]['TOTPOP'] = G.nodes[i]['P0010001'] 
    
    # change meters to kilometers for better numerical performance in gurobi
    th = 1000
    for i in G.nodes:
        G.nodes[i]['area'] /= ( th * th )
        if G.nodes[i]['boundary_node']:
            G.nodes[i]['boundary_perim'] /= th
    for i,j in G.edges:
        G.edges[i,j]['shared_perim'] /= 1000
        
    return G