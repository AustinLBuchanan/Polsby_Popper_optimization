from number_of_districts import congressional_districts_2020
import subprocess

# County-feasible states only...
states = ['AL','AR','ID','IA','KS','ME','MS','MT','NE','NM','WV']
levels = ['county']
objectives = ['aveppbe']
contiguitys = ['lcut','scf','shir']

for state in reversed(states):
    for level in levels:
        for objective in objectives:
            for contiguity in contiguitys:
                subprocess.call(['python3','main.py',state,level,objective,contiguity])

# Tract tests

states = [ key for key in congressional_districts_2020.keys() if congressional_districts_2020[key] >= 2 ]
levels = ['tract']
objectives = ['aveppbe']
contiguitys = ['lcut']

for state in reversed(states):
    for level in levels:
        for objective in objectives:
            for contiguity in contiguitys:
                subprocess.call(['python3','main.py',state,level,objective,contiguity])
