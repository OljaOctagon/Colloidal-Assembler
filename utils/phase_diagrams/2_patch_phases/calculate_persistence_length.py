from scipy import stats
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
    df = pd.read_pickle(filen)
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
        #slist = slist.replace("[", "")
        #slist = slist.replace("]", "")
        #slist = slist.replace(",", "")
        #slist = slist.replace("'","")
        #slist = slist.replace(",", "")
        #slist = slist.replace(" ", "")
        #slist = slist.replace('"', '')

        slist = slist.split("\n")

        slist = [ re.findall("-?\d+\.\d+", item) for item in slist]
        slist = [ np.array(item).astype(float) for item in slist ]

        return slist
    
    #df['sequence'] = df.apply(lambda row: get_sequence(row['sequence']), axis=1)
    #df['bond_vec'] = df.apply(lambda row: get_bond_sequence(row['bond_vec']), axis=1)

    return df

def get_bond_angles(bond_vec):

    #bond_vec = bond_vec.astype(float)
    
    #L = len(bond_vec) - 1
    L = (10-1)  - 1 
    
    bond_angle = np.zeros(L)

    cos_ba = np.zeros(L)
    Denom = np.zeros(L)
    for pos_i in range(0,L):
        for step_j in range(0,L-pos_i):
            norm_a = np.linalg.norm(bond_vec[pos_i])
            norm_b = np.linalg.norm(bond_vec[pos_i+step_j+1])
            cos_i = (np.dot(bond_vec[pos_i], bond_vec[pos_i+step_j+1]))/(norm_a*norm_b)
         
            cos_ba[step_j] +=  cos_i
            Denom[step_j] += 1 

    # average cosine per Length L 
    cos_ba = cos_ba/Denom
    return cos_ba


Csize=10
delta=5
df = read_df("chain_bond_angle.csv")
df = df[df.cluster_type == 'chain']
df = df[(df.cluster_size >= Csize) & (df.cluster_size<=Csize+delta)]
df = df.dropna(subset=['bond_vec'])



df['average_cos'] = df.apply(lambda row: get_bond_angles(row['bond_vec']), axis=1)
for energy in ['8.2']:
    cmap=plt.cm.rainbow
    ncolors = 7
    xlabel="L"
    ylabel="$ (\\langle cos(\\theta)\\rangle)$"
    title="Bond angle decorelation"

    fig, ax  = basic_plot_config(xlabel, ylabel, title, ncolors,cmap) 

    for delta in df.delta.unique():
        dg = df[(df.energy == energy) & (df.delta == delta)]
        arr = dg.average_cos.values
        arr = np.array([ np.array(item) for item in arr])
        arr = - np.log(arr)
        mean_arr = np.mean(arr, axis=0)
        sem = stats.sem(arr, axis=0, ddof=1, nan_policy='propagate')
        
        fit = np.polynomial.polynomial.polyfit(np.arange(1,Csize-1), mean_arr, 1, rcond=None, full=False, w=1/sem)
        print(1/fit)

        plt.errorbar(np.arange(1,Csize-1),mean_arr,
                     sem, label=delta,
                     capsize=4,
                     capthick=1,
                     lw=2,
                     ms=10)
                
    #plt.ylim([-1,0])
    plt.legend(loc='best')
    plt.savefig("cos_decorr_{}.pdf".format(energy))



