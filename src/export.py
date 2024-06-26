import pandas
import csv
import geopandas as gpd
import matplotlib.pyplot as plt
import pathlib

datapath = pathlib.Path("C:/districting-data-2020/")

def draw_single_district( DG, district, zoom=False ):
    labeling = { i : None for i in DG.nodes } if zoom else { i : False for i in DG.nodes }
    for i in district:
        labeling[i] = True
    return export_to_png(DG, labeling, png_filename=None, resize_factor=1)


def export_to_png(DG, labeling, png_filename, resize_factor=3):
    
    print("Exporting to png...")
    
    # Read shapefile from "<state>_<level>.shp"
    filename = DG._state + '_' + DG._level + '.shp'

    # Read geopandas dataframe from file
    df = gpd.read_file( datapath / filename )
    
    # Which district is each node assigned to?
    assignment = [ -1 for i in DG.nodes ]

    # Now add the assignments to a column of the dataframe and map it
    node_with_this_geoid = { DG.nodes[i]['GEOID20'] : i for i in DG.nodes }

    # pick a position u in the dataframe
    for u in range(DG.number_of_nodes()):

        geoid = df['GEOID20'][u]

        # what node in G has this geoid?
        i = node_with_this_geoid[geoid]

        # position u in the dataframe should be given
        # the same district # that node i has in 'labeling'
        assignment[u] = labeling[i]

    # now add the assignments to a column of our dataframe and then map it
    df['assignment'] = assignment
    my_fig = df.plot(column='assignment').get_figure()
    my_fig.set_size_inches(my_fig.get_size_inches()*resize_factor)
    
    plt.axis('off')
    if png_filename is not None:
        my_fig.savefig(png_filename)
    return

def export_to_csv(DG, labeling, csv_filename):
     
    print("Exporting to csv...")
    geoid_to_label = { DG.nodes[i]['GEOID20'] : labeling[i] for i in DG.nodes }

    if DG._level == 'vtd':

        # write our csv file
        with open(csv_filename, 'w', newline="") as csvfile: 
            csvwriter = csv.writer(csvfile) 
            
            # for each row of readfile:
            for geoid, district in geoid_to_label.items():
                row = [ str(geoid), district+1 ]
                csvwriter.writerow(row) 
        return

    # else...
    # read a block assignment file to get the block GEOIDs
    readfile = pandas.read_csv(datapath / (DG._state + '_CD.baf'))

    # then write our block assignment file
    with open(csv_filename, 'w', newline="") as csvfile: 
        csvwriter = csv.writer(csvfile) 
        
        # for each row of readfile:
        for index,row in readfile.iterrows():
            
            # take the block's geoid (and scrap its assignment)
            geoid = row['BLOCKID']
            
            # make it string (and add a leading zero, if it was scrubbed)
            geoidstr = str(geoid)
            if len(geoidstr) < 15: # fixes errors with leading zero
                geoidstr = '0' + geoidstr
                
            # what are the leading digits of the geoid?
            if DG._level == 'county':
                ld = geoidstr[0:5]
            elif DG._level == 'tract':
                ld = geoidstr[0:11]
            elif DG._level == 'blockgroup':
                ld = geoidstr[0:12]
            else: # DG._level == 'block'
                ld = geoidstr
                
            row = [ geoidstr, geoid_to_label[ld]+1 ]
            csvwriter.writerow(row) 
            
    return

def append_dict_as_row(filename, dict_of_elem, field_names):
    # Open file in append mode
    with open(filename, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        dict_writer = csv.DictWriter(write_obj, fieldnames=field_names)
        # Add dictionary as word in the csv
        dict_writer.writerow(dict_of_elem)
        
    