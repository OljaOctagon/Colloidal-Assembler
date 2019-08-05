import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
import re 

def basic_plot_config(xlabel, ylabel, title, ncolors,cmap):
    fig, ax = plt.subplots(figsize=(12,10))

    plt.tick_params(axis='both', labelsize=20) 
    plt.title(title, size=20)
    plt.xlabel(xlabel, size=25)
    plt.ylabel(ylabel, size=25)
    cmap=cmap
    c = cycler('color', cmap(np.linspace(0,1,ncolors)) )
    ax.set_prop_cycle(c)
    return fig,ax


def read_df(filen):
    df = pd.read_csv(filen)
    df = df.drop('Unnamed: 0', axis=1)
    df.bend = (360*df.bend/(2*np.pi))
    df.bend_parallel = (360*df.bend_parallel/(2*np.pi))
    df.bend_non_parallel = (360*df.bend_non_parallel/(2*np.pi))

    def get_sequence(slist):
        slist  = slist.split(",")
        slist = [ item.strip("[") for item in slist]
        slist = [ item.strip("]") for item in slist]
        slist = [ item.strip("'") for item in slist]
        slist = [ item.strip(" ") for item in slist]
        slist = [ item.strip("'") for item in slist]
        return slist

    def get_bond_sequence(slist):
        slist = slist.replace("[", "")
        slist = slist.replace("]", "")
        slist = slist.replace(",", "")
        slist = slist.replace("'","")
        slist = slist.replace(",", "")
        slist = slist.replace(" ", "")
        slist = slist.replace('"', '')
        slist = slist.split("\n")

        try:
            slist = [ re.findall("-?\d+\.\d+", item)[0] for item in slist]
            slist = [ float(item) for item in slist ]

        except:
            slist = None
        return slist 
    
    df['sequence'] = df.apply(lambda row: get_sequence(row['sequence']), axis=1)
    df['bond_angle'] = df.apply(lambda row: get_bond_sequence(row['bond_angle']), axis=1)
    df['bond_vec'] = df.apply(lambda row: get_bond_sequence(row['bond_vec']), axis=1)

    return df

df = read_df("chain_bond_angle.csv")
df = df[df.cluster_type == 'chain']
df = df.explode('bond_angle')

for energy in df.energy.unique():

    ncolors =  8
    cmap = plt.cm.Dark2
    xlabel="bond angle [$\degree$]"
    ylabel="P"
    fig,ax = basic_plot_config(xlabel, ylabel, "", ncolors, cmap)

    for delta in [0.2,0.3,0.5]:
        dg = df[(df.delta == delta) & (df.energy == energy)]
        arr = dg.bond_angle.values
        arr = arr[arr != None ]

        n, bins, patches = plt.hist(arr, 50, alpha=0.7, lw=2, label="$\Delta = {}$".format(1-delta), density=True)


    plt.legend(loc='best', fontsize=20) 
    plt.savefig("bond_angle_distribution_{}.pdf".format(energy), dpi=300)

        
