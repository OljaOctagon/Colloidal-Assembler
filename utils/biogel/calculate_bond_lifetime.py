import numpy as np
import pandas as pd
import argparse
import gel_tools as gt
import networkx as nx
import glob
import os
import matplotlib.pyplot as plt 
import seaborn as sns
import matplotlib as mpl
import matplotlib.style as style
style.use('seaborn-ticks') 
mpl.rcParams['font.family'] = "sans-serif"
#sns.set_context('poster')
plt.rcParams['axes.axisbelow'] = True

plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['figure.titlesize'] = 15 

blue_c ='#9999FF'
red_c ='#FF9999'
purple_c='#C17DCB'
green_c='#61DCB7'


from collections import defaultdict 
from functools import partial
if __name__=='__main__':

    pn_file = "patch_network.dat" 
    # network_arr format: network_arr.shape = ( frame_i, bond_rows_frame_i )
    connections = gt.read_bonds(pn_file)
    Ltime = len(connections)
    connection_soup = np.array([ ci for time_j in range(Ltime) for ci in connections[time_j] ])
    connection_soup = connection_soup[:,:2]
    sorted_connection_soup = np.sort(connection_soup)
    uniques = np.unique(sorted_connection_soup, axis=0)
    bond_sequences = defaultdict(partial(np.zeros, Ltime))

    for time in range(Ltime):
        connections[time] = np.sort(connections[time][:,:2])

    for bond in uniques:
        i = bond[0]
        j = bond[1]
        for time in range(Ltime):
            if np.equal([i,j], connections[time]).all(axis=1).any():
                bond_sequences[(i,j)][time] = 1

    from itertools import groupby
    def len_iter(items):
        return sum(1 for _ in items)

    def consecutive_one(data):
        return [len_iter(run) for val, run in groupby(data) if val]

    lifetimes=[]
    for bond in uniques:
        i = bond[0]
        j = bond[1]
        sum_of_ones = np.sum(bond_sequences[i,j])
        if (sum_of_ones == Ltime):
            print("champ:", i,j)
        lifetimes.extend(consecutive_one(bond_sequences[i,j]))

    fig,ax = plt.subplots()
    plt.xlabel("bond life-time")
    plt.ylabel("P")

    plt.hist(lifetimes, facecolor=purple_c, edgecolor='gray', bins=Ltime, lw=1, alpha=0.7, density=True, hatch='//')
    #plt.legend(loc='best')
    plt.tight_layout()

    plt.savefig("bond_lifetime.pdf")
    plt.show() 

