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
import scipy.optimize
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
import argparse

if __name__=='__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("-file", type=str, help="Bond input file")
    parser.add_argument("-system", type=str, choices=['patchy','biogel'])

    args = parser.parse_args()

    if args.system == 'patchy':
        # network_arr format: network_arr.shape = ( N_time x bond_rows_frame_i x 4 )
        connections = gt.read_bonds(pn_file)
        Ltime = len(connections)
        # all connections as a 2D array with shape  (N_timexN_connections(time) x (4)
        connection_soup = np.array([ ci for time_j in range(Ltime) for ci in connections[time_j] ])
        # throw away the information which patch is bonded
        connection_soup = connection_soup[:,:2]

        # TODO connections also only the bonding elements! otherwise error in 
        # connections[time] = np.sort(connections[time])

    if args.system == 'biogel':
        nbonds=30
        nparticles=340
        Ltime = 2500

        print("read gel")
        arr = pd.read_csv(args.file,delim_whitespace=True).values

        print("prepare data")
        # throw away links between crosslinkers on same polymer
        arr = arr[arr[:,1] != arr[:,3]]

        f_to_number=lambda tup: (tup[0]-1)*nbonds + (tup[1]-1)
        f_to_tuple=lambda c: ((c//nbonds)+1,(c%nbonds)+1)

        a=np.arange(1,nparticles+1)
        b=np.arange(1,nbonds+1)
        arr_combinations = np.array([np.meshgrid(a,b)]).T.reshape(-1,2)
        tuple_combinations = [ (i,j) for [i,j] in arr_combinations]
        map_to_number = list(map(f_to_number,tuple_combinations))

        dict_to_number = dict(zip(tuple_combinations,map_to_number))
        dict_to_tuple  = dict(zip(map_to_number, tuple_combinations))

        arr_id = np.array([
            [time, dict_to_number[(pid1,lid1)],
             dict_to_number[(pid2,lid2)]] for [time, pid1,lid1,pid2,lid2] in arr[:,:5] ])

        connections = [ arr_id[ arr_id[:,0] == time ][:,1:] for time in range(1,Ltime+1)]
        #connection_soup = arr_id[1:]


    # sort the bond ids to get unique links 
    #sorted_connection_soup = np.sort(connection_soup)
    # get the unique links 
    #uniques = np.unique(sorted_connection_soup, axis=0)

    # initialize bond sequences, keys: bond pairs, values: time series bonded states.
    # 1: is bonded 0 not bonded
    bond_sequences = defaultdict(partial(np.zeros, Ltime))
    
    #for time in range(Ltime):
    #    connections[time] = np.sort(connections[time])

    print("get bond info")
    for time in range(Ltime):
        for [i,j] in connections[time]:
            bond_sequences[(i,j)][time] = 1 

    from itertools import groupby
    def len_iter(items):
        return sum(1 for _ in items)

    def consecutive_one(data):
        return [len_iter(run) for val, run in groupby(data) if val]

    print("calculate lifetime")
    # calculate bond life times 
    lifetimes=[]
    for bond in bond_sequences.keys():
        i = bond[0]
        j = bond[1]
        sum_of_ones = np.sum(bond_sequences[i,j])
        if (sum_of_ones == Ltime):
            print("champ:", i,j)
        lifetimes.extend(consecutive_one(bond_sequences[i,j]))

    lifetimes = np.array(lifetimes)
    thresh=100
    lval = len(lifetimes[lifetimes>thresh])/len(lifetimes)

    max_lifetime = np.max(lifetimes)

    print("plot data")
    fig,ax = plt.subplots()
    plt.xlabel("bond life-time")
    plt.ylabel("P")
    plt.xlim((0,100)) 

    n,bins,patches = plt.hist(lifetimes,
            edgecolor='k', bins=(np.linspace(0,Ltime,Ltime+1)+0.1),
            lw=2, alpha=1, 
            density=True)

    fig, ax  = plt.subplots()
    plt.yscale("log")
    x = bins[1:]

    # fit exponential
    def monoExp(x,l,t,):
        return l * np.exp(-t * x) 

    p0 = (1, 1) # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(monoExp, x, n, p0)
    l,t = params

    plt.plot(x,n,c=purple_c,linestyle='dotted',label='data')
    plt.plot(x,monoExp(x,l,t), c=blue_c, linestyle='--', label='biexp. fit')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.ylim(0.00001,1)
    plt.savefig("bond_lifetime_fit.pdf")

    def cumulative_monoExp(xs,l,t):
             return 1 - l*np.exp(-t*xs)  

    percentil = 1 - cumulative_monoExp(thresh,l,t)

    data = np.array([l,t,max_lifetime,lval-percentil])
    np.savetxt("observables_bond_lifetime.dat", data, delimiter=',', newline='\n')
    np.savetxt("bond_lifetime.dat",n, newline='\n') 
