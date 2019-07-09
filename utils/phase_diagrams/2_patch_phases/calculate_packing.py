import numpy as np
import pandas as pd
import argparse
import glob
import os
import re
import networkx as nx


filenpt = 'NPT_OUT.txt'
dirlist = glob.glob("mu_*")
pwd = os.getcwd()


features = ['mu', 'energy', 'topology', 'delta', 'packing_fraction'] 

df = pd.DataFrame(columns=features)

for dir in dirlist:
    print(dir)
    numbers = re.findall(r"[-+]?\d*\.\d+|\d+", dir)
    mu = numbers[0]
    energy = numbers[1]
    delta = numbers[2]
    topology = 'symm'
    p = 'Asymm'
    if re.search(p, dir):
        topology = 'Asymm'

    try:

        arr_npt = pd.read_csv('{}/{}'.format(dir, filenpt),
                                header=None,
                                delim_whitespace=True).values

        packing_fraction = arr_npt[-1,4]

        
        df = pd.concat([df,pd.DataFrame({'mu': [mu],
                                        'energy': [energy],
                                        'topology': [topology],
                                        'delta': [delta],
                                        'packing_fraction': [packing_fraction]})])

    except:
        pass

with open('packing_chains.csv', 'w') as f:
   df.to_csv(f)
