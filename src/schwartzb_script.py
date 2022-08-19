import subprocess
import sys
from number_of_districts import congressional_districts_2020

# Alabama, Arkansas, Idaho, Iowa, Kansas, Maine, Mississippi, Montana, Nebraska, New Mexico, and West Virginia.
states = ['AL','AR','ID','IA','KS','ME','MS','MT','NE','NM','WV']
levels = ['county']
objectives = ['schwartzb']
contiguitys = ['lcut','scf','shir']

for state in reversed(states):
    for level in levels:
        for objective in objectives:
            for contiguity in contiguitys:
                subprocess.call(['python','main.py',state,level,objective,contiguity] + sys.argv[1:])

# Tract tests

states = [ key for key in congressional_districts_2020.keys() if congressional_districts_2020[key] >= 2 ]
levels = ['tract']
objectives = ['schwartzb']
contiguitys = ['lcut', 'scf', 'shir']

for state in reversed(states):
    for level in levels:
        for objective in objectives:
            for contiguity in contiguitys:
                subprocess.call(['python','main.py',state,level,objective,contiguity] + sys.argv[1:])
