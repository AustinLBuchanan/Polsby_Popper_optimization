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

def number_of_counties_split(G, districts, verbose=False):
    fs = fips_support(G, districts)
    if verbose:
        print("\nCounties split (by fips code):", [ f for f in fs.keys() if len(fs[f]) > 1 ])
    return sum( 1 for f in fs.keys() if len(fs[f]) > 1 )
    
def number_of_county_splits(G, districts, verbose=False):
    fs = fips_support(G, districts)
    if verbose:
        print("County splits (by fips code):", { f : len(fs[f])-1  for f in fs.keys() if len(fs[f]) > 1 })
    return sum( ( len(fs[f]) - 1 ) for f in fs.keys() )
    
def polsby_popper(G, district, label):
    area = sum( G.nodes[i]['area'] for i in district )
    perim = sum( G.edges[u,v]['shared_perim'] for u in district for v in G.neighbors(u) if label[u]!=label[v] )
    perim += sum( G.nodes[i]['boundary_perim'] for i in district if G.nodes[i]['boundary_node'] ) 
    return 4 * math.pi * area / ( perim * perim )

def average_polsby_popper(G, districts, verbose=False):
    label = { i : j for j in range(len(districts)) for i in districts[j] }
    if verbose:
        print("\nDistrict Polsby-Popper scores:")
        for p in range(len(districts)):
            print(p, round(polsby_popper(G, districts[p], label),4) )
    return sum( polsby_popper(G, district, label) for district in districts ) / len(districts) 

def number_of_gingles_districts(G, districts, minority, verbose=False):
    codes = get_census_codes(minority)
    for i in G.nodes:
        G.nodes[i]['VAP'] = G.nodes[i]['P0030001']
        G.nodes[i]['MVAP'] = sum( G.nodes[i][code] for code in codes )
    gingles_count = 0
    if verbose:
        print(f"\nDistrict {minority} percentages:")
    for p in range(len(districts)):
        district = districts[p]
        vap = sum( G.nodes[i]['VAP'] for i in district )
        mvap = sum( G.nodes[i]['MVAP'] for i in district )
        if mvap >= 0.5 * vap:
            gingles_count += 1
        if verbose:
            flag = '***' if mvap >= 0.5 * vap else ''
            print(p, round(100*mvap/vap,2), flag)
    return gingles_count

def report_metrics(G, districts, minority=None, verbose=False):
    if minority is not None:
        gingles = number_of_gingles_districts(G, districts, minority, verbose=verbose)
        print(f"-> {gingles} majority-{minority} districts")
    print(f"-> {number_of_counties_split(G, districts, verbose=verbose)} counties split a total of {number_of_county_splits(G, districts, verbose=verbose)} times")
    avepp = round( average_polsby_popper(G, districts, verbose=verbose), 4)
    print("-> average Polsby-Popper score of",avepp)
    return
    