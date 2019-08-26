import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from cycler import cycler

df = pd.read_pickle('chains_end_to_end.csv')
df = df[df.cluster_type == 'chain']
Lx = 1
sigma_r = 0.1

def small_d(delta):
     x = np.cos(np.pi/3.)*np.fabs(0.5-delta)
     y = np.sin(np.pi/3)*np.fabs(0.5-delta)
     return np.sqrt(np.power(0.5+x,2) + np.power(y,2))*2 + 0.1

def big_d(delta):
     x = np.cos(np.pi/3.)*np.fabs(0.5-delta)
     y = np.sin(np.pi/3)*np.fabs(0.5-delta)
     return np.sqrt(np.power(0.5-x,2) + np.power(y,2))*2 + 0.1

dp_sb =  lambda delta: small_d(delta) if delta<=0.5 else big_d(delta)
dp_bs =  lambda delta: big_d(delta) if delta<=0.5 else small_d(delta)
dp_ss =  lambda delta: small_d(delta) if delta<=0.5 else small_d(delta)
dp_bb =  lambda delta: big_d(delta) if delta<=0.5 else big_d(delta)
dp_on = lambda delta: 1.1
dnp_on = lambda delta: np.sin(np.pi/3) + 0.1
dnp = lambda delta: small_d(delta)/0.5 + big_d(delta)/0.5

get_max_bond_length = {'p-off-sb': dp_sb,
                       'p-off-bs': dp_bs,
                       'p-off-s':  dp_ss,
                       'p-off-b': dp_bb,
                       'p-on': dp_on,
                       'np-on': dnp_on,
                       'np-off':dnp}

def get_contour(delta, cluster_size, sequence):
    unique, counts = np.unique(sequence, return_counts=True)
    L_contour = 0
    for i,bond_type in enumerate(unique):
        L_contour += get_max_bond_length[bond_type](delta)*counts[i]

    #add half at the beginning and end
    L_contour += (get_max_bond_length[sequence[0]](delta) - 0.1)/2.
    L_contour += (get_max_bond_length[sequence[-1]](delta) - 0.1)/2.

    return L_contour

df['L_contour'] = df.apply(
        lambda row: get_contour(float(row['delta']), row['cluster_size'], row['sequence']), axis=1)
df['scaled_e_to_e'] = df.end_to_end_distance/df['L_contour']

dg = df.groupby(['energy','delta','cluster_size']).scaled_e_to_e.agg(['mean', 'std']).reset_index() 

df[["kink_density"]] = df[["kink_density"]].apply(pd.to_numeric)
df = df[(df.kink_density>=0.2)& (df.kink_density<0.3)]

for energy in [10.2]:
    fig,ax = plt.subplots()
    cmap=plt.cm.rainbow
    c = cycler('color', cmap(np.linspace(0,1,8)) )
    ax.set_prop_cycle(c)
    for delta in np.arange(0.5):
        delta = np.round(delta,1)
        print(df.head())
        dg2 = df[(df.delta == str(delta)) & (df.energy==str(energy))]
        dg2 =  df.groupby(['cluster_size']).scaled_e_to_e.agg(['mean', 'std']).reset_index() 
        plt.errorbar(dg2.cluster_size.values,
                    dg2['mean'].values,
                    dg2['std'].values,
                    capsize=4,
                    capthick=1,
                    lw=2,
                    ms=10, 
                    label="{}".format(delta))

            
    plt.xlabel("$L$", size=20)
    plt.ylabel("end-to-end distance", size=20)
    plt.xlim([0,20])
    plt.legend(loc='best')

    plt.savefig("{}_end_to_end_test.png".format(energy))


for energy in np.arange(5.2,11.2,1):
    fig,ax = plt.subplots()
    cmap=plt.cm.rainbow
    c = cycler('color', cmap(np.linspace(0,1,8)) )
    ax.set_prop_cycle(c)

    for delta in np.arange(0.5,0.9,0.1):
        delta = np.round(delta,1)
        arr = df[(df.delta == str(delta)) & (df.energy==str(energy))].kink_density.values

        n,bins,patches = plt.hist(arr, label=str(delta))

    plt.xlabel("kink density", size=20, alpha=0.6)
    plt.ylabel("P", size=20)
    plt.legend(loc='best')

    plt.savefig("{}_kink_density.png".format(energy))
