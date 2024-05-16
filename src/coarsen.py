import networkx as nx

#########################################################
# Generic coarsen function merges each part into a node.
#########################################################

def graph_coarsen(G, partition, verbose=False):

    # check that partition is valid
    assert_valid_partition(G, partition)
    
    # create coarsened graph
    CG = nx.Graph()
    
    # add one node for each part 
    nodes = list( range( len( partition ) ) )
    CG.add_nodes_from(nodes)
    
    # add edges for adjacent parts
    assignment = { i : j for j in range(len(partition)) for i in partition[j] }
    edges = { (min(assignment[u],assignment[v]), max(assignment[u],assignment[v])) for u,v in G.edges if assignment[u] != assignment[v] } 
    CG.add_edges_from( edges )
    
    for i in G.nodes:
        if not G.nodes[i]['boundary_node']:
            G.nodes[i]['boundary_perim'] = 0
    
    additive_attr = additive_attr_default.copy()
    v = nx.utils.arbitrary_element(G.nodes)
    for attr in additive_attr_default:
        try:
            G.nodes[v][attr]
        except:
            additive_attr.remove(attr)
            
    # add node attributes
    for c in CG.nodes:
        CG.nodes[c]['boundary_node'] = False
        for attr in additive_attr:
            CG.nodes[c][attr] = 0
        
    for i in G.nodes:
        j = assignment[i]
        CG.nodes[j]['boundary_node'] |= G.nodes[i]['boundary_node'] 
        for attr in additive_attr:
            CG.nodes[j][attr] += G.nodes[i][attr]
            
    # add edge attributes
    for i,j in CG.edges:
        CG.edges[i,j]['shared_perim'] = 0
    
    for u,v in G.edges:
        i = min(assignment[u],assignment[v])
        j = max(assignment[u],assignment[v])
        if i != j:
            CG.edges[i,j]['shared_perim'] += G.edges[u,v]['shared_perim']
        
    if verbose:
        print(f"In graph_coarsen(), reduced n = {G.number_of_nodes()} to n = {CG.number_of_nodes()}")
        print(f" and m = {G.number_of_edges()} to m = {CG.number_of_edges()}")
    
    return CG


#########################################################
# Collapse each county into a single node.
# However, do not collapse those in granular_fips.
# If separate_county_components, collapse each 
#    *component* of each county into a node.
#########################################################

def graph_coarsen_by_county(G, granular_fips=list(), separate_county_components=True, verbose=False):

    # list of all fips codes in G
    fips = list( { G.nodes[i]['GEOID20'][0:5] for i in G.nodes } )
    
    # create initial partition by county
    partition_by_county = [ list() for f in fips ]
    fips_map = { fips[p] : p for p in range(len(fips)) }
    for i in G.nodes:
        f = G.nodes[i]['GEOID20'][0:5]
        partition_by_county[fips_map[f]].append(i)
    
    # create final partition
    partition = list()
    partition_geoid = list()
    for p in range(len(fips)):
        f = fips[p]
        part = partition_by_county[p]
        
        if f in granular_fips:
            for i in part:
                partition.append( [i] )
                partition_geoid.append( G.nodes[i]['GEOID20'] )
        elif separate_county_components:
            count = 0
            for comp in nx.connected_components(G.subgraph(part)):
                partition.append( list(comp) )
                partition_geoid.append( f ) # all county components get same geoid(!)
                count += 1
            if verbose and count > 1:
                print(f"Created {count} nodes for county with fips = {f}, one for each of its components.")
        else:
            partition.append(part)
            partition_geoid.append( f )

    CG = graph_coarsen(G, partition, verbose=verbose)
    for p in CG.nodes:
        CG.nodes[p]['GEOID20'] = partition_geoid[p]
    
    return (CG, partition)


# G = input graph
# plan = initial plan that we are coarsening
# h = number of hops to use for coarsening
# by_county = when coarsening, should we separate the coarsened nodes by their counties?
#
def hop_coarsen(G, plan, h, by_county, verbose=False):

    partition = list()
    district_bool = { i : False for i in G.nodes }
    coarsened_plan = [ list() for j in range(len(plan)) ]
    coarsened_nodes = 0
    
    for j in range(len(plan)):
    
        # for speed, represent the district in boolean form
        for i in plan[j]:
            district_bool[i] = True
        
        complement = [ i for i in district_bool.keys() if not district_bool[i] ]
        dist = distance_to_vertex_set(G, complement)
        
        # nodes near the district boundary will NOT be coarsened
        for i in plan[j]:
            if dist[i] <= h:
                partition.append([i]) 
                coarsened_plan[j].append(coarsened_nodes)
                coarsened_nodes += 1
        
        # nodes far from the district boundary WILL be coarsened (by connected components)
        interior = [ i for i in plan[j] if dist[i] > h ]
        if not by_county:
            for comp in nx.connected_components( G.subgraph(interior) ):
                partition.append( list(comp) )
                coarsened_plan[j].append(coarsened_nodes)
                coarsened_nodes += 1
        else:
            fips = { G.nodes[i]['GEOID20'][0:5] for i in interior }
            for f in fips:
                fips_interior = [ i for i in interior if G.nodes[i]['GEOID20'][0:5] == f ]
                for comp in nx.connected_components( G.subgraph(fips_interior) ):
                    partition.append( list(comp) )
                    coarsened_plan[j].append(coarsened_nodes)
                    coarsened_nodes += 1
        
        # reset district_bool for the next district
        for i in plan[j]:
            district_bool[i] = False
    
    CG = graph_coarsen(G, partition, verbose=verbose)
    if by_county:
        for i in CG.nodes:
            v = partition[i][0]
            CG.nodes[i]['GEOID20'] = G.nodes[v]['GEOID20'][0:5]
    
    assert CG.number_of_nodes() == coarsened_nodes
    return ( CG, partition, coarsened_plan )


# dist[i] = the distance from node i to the set of nodes S
# dist[i]=0 if i in S
# dist[i]=1 if i in N(S)
# dist[i]=2 if i in N^2(S)
# ...
# dist[i]=None if no path from i to S
#
def distance_to_vertex_set(G, S, cutoff=None):
    if cutoff==None:
        cutoff = G.number_of_nodes()
    
    dist = dict()
    touched = { i : False for i in G.nodes }
    child = S
    for s in S:
        touched[s] = True
    
    for h in range(1+cutoff):
        parent = child
        child = list()
        for p in parent:
            dist[p] = h
            for n in G.neighbors(p):
                if not touched[n]:
                    child.append(n)
                    touched[n] = True
    return dist


# Returns subgraph of DG induced by nodes.
# Importantly, also updates boundary status/length.
# Thus, it is different than DG.subgraph(nodes). 
#
def subgraph(DG, nodes):
    DH = DG.subgraph(nodes).copy()
    nodes_bool = { i : False for i in DG.nodes }
    for i in nodes:
        nodes_bool[i] = True
    
    for i in nodes:
        
        # which nodes are boundary in DH ?
        if not DH.nodes[i]['boundary_node']:
            for j in DG.neighbors(i):
                if not nodes_bool[j]:
                    DH.nodes[i]['boundary_node'] = True
                    DH.nodes[i]['boundary_perim'] = 0
                    break # break loop over neighbors j
        
        #  and what is their new exterior boundary length?
        for j in DG.neighbors(i):
            if not nodes_bool[j]:
                DH.nodes[i]['boundary_perim'] += DG.edges[i,j]['shared_perim']
    return DH


def assert_valid_partition(G, partition):
    assert set(G.nodes) == { i for j in range(len(partition)) for i in partition[j] }, "partition and graph have different node sets"
    assert G.number_of_nodes() == sum( len(part) for part in partition ), "parts are not disjoint"
    assert min( len(part) for part in partition ) > 0, "partition has an empty part"
    for part in partition:
        assert type(part) in {list, set}, "not all parts are lists/sets"
    return

additive_attr_default = ['area','boundary_perim','TOTPOP','P0010001', 'P0010002', 'P0010003', 'P0010004', 'P0010005', 'P0010006', 'P0010007', 'P0010008', 'P0010009', 'P0010010', 'P0010011', 'P0010012', 'P0010013', 'P0010014', 'P0010015', 'P0010016', 'P0010017', 'P0010018', 'P0010019', 'P0010020', 'P0010021', 'P0010022', 'P0010023', 'P0010024', 'P0010025', 'P0010026', 'P0010027', 'P0010028', 'P0010029', 'P0010030', 'P0010031', 'P0010032', 'P0010033', 'P0010034', 'P0010035', 'P0010036', 'P0010037', 'P0010038', 'P0010039', 'P0010040', 'P0010041', 'P0010042', 'P0010043', 'P0010044', 'P0010045', 'P0010046', 'P0010047', 'P0010048', 'P0010049', 'P0010050', 'P0010051', 'P0010052', 'P0010053', 'P0010054', 'P0010055', 'P0010056', 'P0010057', 'P0010058', 'P0010059', 'P0010060', 'P0010061', 'P0010062', 'P0010063', 'P0010064', 'P0010065', 'P0010066', 'P0010067', 'P0010068', 'P0010069', 'P0010070', 'P0010071', 'P0020001', 'P0020002', 'P0020003', 'P0020004', 'P0020005', 'P0020006', 'P0020007', 'P0020008', 'P0020009', 'P0020010', 'P0020011', 'P0020012', 'P0020013', 'P0020014', 'P0020015', 'P0020016', 'P0020017', 'P0020018', 'P0020019', 'P0020020', 'P0020021', 'P0020022', 'P0020023', 'P0020024', 'P0020025', 'P0020026', 'P0020027', 'P0020028', 'P0020029', 'P0020030', 'P0020031', 'P0020032', 'P0020033', 'P0020034', 'P0020035', 'P0020036', 'P0020037', 'P0020038', 'P0020039', 'P0020040', 'P0020041', 'P0020042', 'P0020043', 'P0020044', 'P0020045', 'P0020046', 'P0020047', 'P0020048', 'P0020049', 'P0020050', 'P0020051', 'P0020052', 'P0020053', 'P0020054', 'P0020055', 'P0020056', 'P0020057', 'P0020058', 'P0020059', 'P0020060', 'P0020061', 'P0020062', 'P0020063', 'P0020064', 'P0020065', 'P0020066', 'P0020067', 'P0020068', 'P0020069', 'P0020070', 'P0020071', 'P0020072', 'P0020073', 'P0030001', 'P0030002', 'P0030003', 'P0030004', 'P0030005', 'P0030006', 'P0030007', 'P0030008', 'P0030009', 'P0030010', 'P0030011', 'P0030012', 'P0030013', 'P0030014', 'P0030015', 'P0030016', 'P0030017', 'P0030018', 'P0030019', 'P0030020', 'P0030021', 'P0030022', 'P0030023', 'P0030024', 'P0030025', 'P0030026', 'P0030027', 'P0030028', 'P0030029', 'P0030030', 'P0030031', 'P0030032', 'P0030033', 'P0030034', 'P0030035', 'P0030036', 'P0030037', 'P0030038', 'P0030039', 'P0030040', 'P0030041', 'P0030042', 'P0030043', 'P0030044', 'P0030045', 'P0030046', 'P0030047', 'P0030048', 'P0030049', 'P0030050', 'P0030051', 'P0030052', 'P0030053', 'P0030054', 'P0030055', 'P0030056', 'P0030057', 'P0030058', 'P0030059', 'P0030060', 'P0030061', 'P0030062', 'P0030063', 'P0030064', 'P0030065', 'P0030066', 'P0030067', 'P0030068', 'P0030069', 'P0030070', 'P0030071', 'P0040001', 'P0040002', 'P0040003', 'P0040004', 'P0040005', 'P0040006', 'P0040007', 'P0040008', 'P0040009', 'P0040010', 'P0040011', 'P0040012', 'P0040013', 'P0040014', 'P0040015', 'P0040016', 'P0040017', 'P0040018', 'P0040019', 'P0040020', 'P0040021', 'P0040022', 'P0040023', 'P0040024', 'P0040025', 'P0040026', 'P0040027', 'P0040028', 'P0040029', 'P0040030', 'P0040031', 'P0040032', 'P0040033', 'P0040034', 'P0040035', 'P0040036', 'P0040037', 'P0040038', 'P0040039', 'P0040040', 'P0040041', 'P0040042', 'P0040043', 'P0040044', 'P0040045', 'P0040046', 'P0040047', 'P0040048', 'P0040049', 'P0040050', 'P0040051', 'P0040052', 'P0040053', 'P0040054', 'P0040055', 'P0040056', 'P0040057', 'P0040058', 'P0040059', 'P0040060', 'P0040061', 'P0040062', 'P0040063', 'P0040064', 'P0040065', 'P0040066', 'P0040067', 'P0040068', 'P0040069', 'P0040070', 'P0040071', 'P0040072', 'P0040073', 'H0010001', 'H0010002', 'H0010003', 'P0050001', 'P0050002', 'P0050003', 'P0050004', 'P0050005', 'P0050006', 'P0050007', 'P0050008', 'P0050009', 'P0050010', 'VAP', 'MVAP']
