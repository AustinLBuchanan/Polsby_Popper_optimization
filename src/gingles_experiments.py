import math
from read import read_graph_from_json
from number_of_districts import number_of_districts
from gingles import gingles
from export import export_to_baf

filepath = 'C:\\districting-data-2020\\'
states = ['LA', 'MS', 'AL', 'GA']
district_types = ['SS', 'SH']

level = 'vtd'
minority = 'Black'
deviation = 0.10      # 0.10 means 10% means +/-5%

for state in states:
    for district_type in district_types:
        
        print("\n******************************")
        print(f"STARTING {state} {district_type}")
        print("******************************\n")
        
        # read graph
        filename = state + '_' + level + '.json'
        G = read_graph_from_json( filepath + filename )
        for i in G.nodes:
            G.nodes[i]['TOTPOP'] = G.nodes[i]['P0010001'] 
        print("number of nodes, edges:",G.number_of_nodes(),G.number_of_edges())
                    
        # set params
        G._k = number_of_districts[state, district_type]
        ideal_population = sum( G.nodes[i]['TOTPOP'] for i in G.nodes ) / G._k
        G._L = math.ceil(  ideal_population * (1-deviation/2) )
        G._U = math.floor( ideal_population * (1+deviation/2) )
        print("k, L, U =",G._k,G._L,G._U)
        
        # apply Gingles code
        h = 1 if district_type=='SH' else 2
        districts = gingles(G, minority=minority, h=h)
        
        # save districts 
        G._state = state
        G._level = level
        csv_filename = "gingles_" + state + "_" + district_type + ".csv"
        labeling = { i : j for j in range(len(districts)) for i in districts[j] }
        export_to_baf(G, labeling, csv_filename)
                
