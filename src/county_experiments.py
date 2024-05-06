import subprocess

# Alabama, Arkansas, Idaho, Iowa, Kansas, Maine, Mississippi, Montana, Nebraska, New Mexico, and West Virginia.
states = ['AR','ID','IA','KS','ME','MS','MT','NE','NM','WV']
levels = ['county']
objectives = ['invpp','avepp','schwartzb']
contiguitys = ['lcut']

for state in states:
    for level in levels:
        for objective in objectives:
            for contiguity in contiguitys:
                subprocess.call(['python','county_main.py',state,level,objective,contiguity])
                     
                
