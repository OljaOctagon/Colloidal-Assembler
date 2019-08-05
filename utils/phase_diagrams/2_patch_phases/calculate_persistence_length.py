
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
    
    df['sequence'] = df.apply(lambda row: get_sequence(row['sequence']), axis=1)
    df['bond_vec'] = df.apply(lambda row: get_bond_sequence(row['bond_vec']), axis=1)

    return df

def get_bond_angles(bond_vec):

    L = len(bond_vec) - 1
    #print("L bond_vec", L)
    bond_angle = np.zeros(L)

    cos_ba = np.zeros(L)
    Denom = np.zeros(L)

    #print("L cos", len(cos_ba))
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
df = read_df("chain_bond_angle.csv")
df = df[df.cluster_type == 'chain']
df = df[df.cluster_size == Csize]

df = df[df.bond_vec != None]


df['average_cos'] = df.apply(lambda row: get_bond_angles(row['bond_vec']), axis=1)

for energy in [7.2,8.2]:

    cmap=plt.cm.rainbow
    ncolors = 7
    xlabel="L"
    ylabel="$ log(\\langle cos(\\theta)\\rangle)$"
    title="Bond angle decorelation"

    fig, ax  = basic_plot_config(xlabel, ylabel, title, ncolors,cmap) 

    for delta in df.delta.unique():
        dg = df[(df.energy == energy) & (df.delta == delta)]
        arr = dg.average_cos.values
        arr = np.array([ np.array(item) for item in arr])
        print(arr.shape)
        mean_arr = np.mean(arr, axis=0)
        std_arr = np.std(arr, axis=0)

        plt.errorbar(np.arange(1,Csize-1),np.log(mean_arr),
                     std_arr, label=delta,
                     capsize=4,
                     capthick=1,
                     lw=2,
                     ms=10)

    plt.legend(loc='best')
    plt.savefig("log_cos_theta_decorrelation_{}.pdf".format(energy))

