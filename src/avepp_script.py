from number_of_districts import congressional_districts_2020
import subprocess

# do county first

# County-feasible states only...
states = ['AL','AR','ID','IA','KS','ME','MS','MT','NE','NM','WV']
levels = ['county']
objectives = ['avepp']
contiguitys = ['lcut','scf','shir']

for state in reversed(states):
    for level in levels:
        for objective in objectives:
            for contiguity in contiguitys:
                subprocess.call(['python','main.py',state,level,objective,contiguity])
                
                
# now tract

states = [ key for key in congressional_districts_2020.keys() if congressional_districts_2020[key] >= 2 ]
levels = ['tract']
objectives = ['avepp']
contiguitys = ['lcut']

for state in reversed(states):
    for level in levels:
        for objective in objectives:
            for contiguity in contiguitys:
                subprocess.call(['python','main.py',state,level,objective,contiguity])
                     
                