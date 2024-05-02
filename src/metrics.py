import math
from census_codes import get_census_codes

def fips_support(G, districts):
    fips = { G.nodes[i]['GEOID20'][0:5] for i in G.nodes }
    fs = { f : list() for f in fips }
    for j in range(len(districts)):
        for i in districts[j]:
            f = G.nodes[i]['GEOID20'][0:5]
            if j not in fs[f]:
                fs[f].append(j)
    return fs

def number_of_counties_split(G, districts):
    fs = fips_support(G, districts)
    return sum( 1 for f in fs.keys() if len(fs[f]) > 1 )
    
def number_of_county_splits(G, districts):
    fs = fips_support(G, districts)
    return sum( ( len(fs[f]) - 1 ) for f in fs.keys() )
    
def polsby_popper(G, district, label):
    area = sum( G.nodes[i]['area'] for i in district )
    perim = sum( G.edges[u,v]['shared_perim'] for u in district for v in G.neighbors(u) if label[u]!=label[v] )
    perim += sum( G.nodes[i]['boundary_perim'] for i in district if G.nodes[i]['boundary_node'] ) 
    return 4 * math.pi * area / ( perim * perim )

def average_polsby_popper(G, districts):
    label = { i : j for j in range(len(districts)) for i in districts[j] }
    return sum( polsby_popper(G, district, label) for district in districts ) / len(districts) 

def number_of_gingles_districts(G, districts, minority):
    codes = get_census_codes(minority)
    for i in G.nodes:
        G.nodes[i]['VAP'] = G.nodes[i]['P0030001']
        G.nodes[i]['MVAP'] = sum( G.nodes[i][code] for code in codes )
    gingles_count = 0
    for district in districts:
        vap = sum( G.nodes[i]['VAP'] for i in district )
        mvap = sum( G.nodes[i]['MVAP'] for i in district )
        if mvap >= 0.5 * vap:
            gingles_count += 1
    return gingles_count

def report_metrics(G, districts, minority):
    gingles = number_of_gingles_districts(G, districts, minority)
    s1 = number_of_counties_split(G, districts)
    s2 = number_of_county_splits(G, districts)
    avepp = round( average_polsby_popper(G, districts), 4)
    print(f"  {gingles} majority-{minority} districts")
    print(f"  {s1} counties split a total of {s2} times")
    print("  average Polsby-Popper score of",avepp)
    return
    