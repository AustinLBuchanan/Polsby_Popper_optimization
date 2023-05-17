import subprocess

# Alabama, Arkansas, Idaho, Iowa, Kansas, Maine, Mississippi, Montana, Nebraska, New Mexico, and West Virginia.
states = ['AR','ID','IA','KS','ME','MS','MT','NE','NM','WV']
levels = ['county']
objectives = ['invpp','avepp','schwartzb']
contiguitys = ['shir']

for state in states:
    for level in levels:
        for objective in objectives:
            for contiguity in contiguitys:
                subprocess.call(['python','main.py',state,level,objective,contiguity])
                     
                