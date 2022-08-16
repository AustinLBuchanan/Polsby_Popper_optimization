from number_of_districts import congressional_districts_2020
import subprocess

states = [ key for key in congressional_districts_2020.keys() if congressional_districts_2020[key] >= 2 ]
levels = ['tract']
objectives = ['cut', 'perim', 'invpp', 'schwartzb']
contiguitys = ['lcut'] #,'scf','shir']

for state in reversed(states):
    for level in levels:
        for objective in objectives:
            for contiguity in contiguitys:
                subprocess.call(['python','main.py',state,level,objective,contiguity])
                     
                