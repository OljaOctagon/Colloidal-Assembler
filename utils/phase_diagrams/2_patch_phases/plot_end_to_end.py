import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from cycler import cycler


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
                       'np-off':dnp,
                       'np': dp_sb,
                       'p': dp_sb}

def get_contour(delta, cluster_size, sequence):
    unique, counts = np.unique(sequence, return_counts=True)
    L_contour = 0
    for i,bond_type in enumerate(unique):
        L_contour += get_max_bond_length[bond_type](delta)*counts[i]

    #add half at the beginning and end
    L_contour += (get_max_bond_length[sequence[0]](delta) - 0.1)/2.
    L_contour += (get_max_bond_length[sequence[-1]](delta) - 0.1)/2.

    return L_contour

def read_df(filen):
     df = pd.read_pickle(filen)
     df = df[df.cluster_type == 'chain']
     Lx = 1
     sigma_r = 0.1
     df['L_contour'] = df.apply(
          lambda row: get_contour(float(row['delta']), row['cluster_size'], row['sequence']), axis=1)
     df['scaled_e_to_e'] = df.end_to_end_distance/df['L_contour']
     dg = df.groupby(['energy','delta','cluster_size']).scaled_e_to_e.agg(['mean', 'std']).reset_index() 

     return dg


def plot_end_to_end(dg, ax, energy, name):
     cmap=plt.cm.jet
     c = cycler('color', cmap(np.linspace(0,1,8)) )
     ax.set_prop_cycle(c)
     ax.set_title(name, size=18)
     for delta in np.arange(0.2,0.9,0.1):
          delta = np.round(delta,1)
          dg2 = dg[(dg.delta == str(delta)) & (dg.energy==str(energy))]
          ax.errorbar(dg2.cluster_size.values,
                         dg2['mean'].values,
                         dg2['std'].values,
                         capsize=4,
                         capthick=1,
                         lw=3,
                         ms=10, 
                         label="{}".format(delta))

     ax.tick_params(axis='both', labelsize=16)
     ax.set_xticks(np.arange(2,21,2))
     ax.set_xlabel("$L$", size=20)
     ax.set_ylabel("$\\langle x_{s} - x_{e}\\rangle/l_{c}$", size=20)
     ax.set_xlim([0,20])
     ax.set_ylim([0.5,1.0])
     ax.legend(loc='best')


if __name__ == "__main__":

     f_s = 'chains_end_to_end.csv'
     f_as = '/Users/ada/Documents/Code_Development_2016/2D_patchy/chain_analysis/small_data_set/rr_chains/lege_batterie/asymm/mu_0.3/chain_bond_angle.csv'

     dg_s = read_df(f_s)
     dg_as = read_df(f_as)

     fig,(ax1,ax2) = plt.subplots(ncols=2,nrows=1, figsize=(15,7))

     energy = 8.2 
     plot_end_to_end(dg_s, ax1, energy, 'End-end distance for pl-s')
     plot_end_to_end(dg_as,ax2, energy, 'End-end distance for pl-as')


     plt.savefig("{}_pl_end_to_end.pdf".format(energy))

