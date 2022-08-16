import subprocess

# Alabama, Arkansas, Idaho, Iowa, Kansas, Maine, Mississippi, Montana, Nebraska, New Mexico, and West Virginia.
states = ['AL','AR','ID','IA','KS','ME','MS','MT','NE','NM','WV']
#states = ['WV']
levels = ['county']
objectives = ['schwartzb']
contiguitys = ['lcut','scf','shir']

for state in reversed(states):
    for level in levels:
        for objective in objectives:
            for contiguity in contiguitys:
                subprocess.call(['python','main.py',state,level,objective,contiguity])
                     
                
