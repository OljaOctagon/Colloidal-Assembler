from scipy import stats
import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
import re 
import argparse 
from scipy.optimize import curve_fit



name_dict = {
     'symm': 'pl-s',
     'asymm': 'pl-as'
}

def small_d(delta):
     x = np.cos(np.pi/3.)*np.fabs(0.5-delta)
     y = np.sin(np.pi/3)*np.fabs(0.5-delta)

     dist_x = 1 + 2*x
     dist_y = 2*y 
     dist = np.sqrt(np.power(dist_x,2) + np.power(dist_y,2))

     normed_x = dist_x/dist
     normed_y = dist_y/dist 
     return np.array([normed_x, normed_y, dist])

def big_d(delta):
     x = np.cos(np.pi/3.)*np.fabs(0.5-delta)
     y = np.sin(np.pi/3)*np.fabs(0.5-delta)

     dist_x = 1 - 2*x
     dist_y = -2*y 
     dist = np.sqrt(np.power(dist_x,2) + np.power(dist_y,2))

     normed_x = dist_x/dist
     normed_y = dist_y/dist 

     return np.array([normed_x, normed_y, dist])

dp_sb =  lambda delta: small_d(delta)[2] if delta<=0.5 else big_d(delta)[2]
dp_bs =  lambda delta: big_d(delta)[2] if delta<=0.5 else small_d(delta)[2]
dp_on = lambda delta: 1.0

def distance_pl_as(delta,a,b,c):
    return dp_sb(delta)

def distance_pl_s(delta,a,b,c):
    d =  (a*dp_sb(delta) + b*dp_bs(delta) + c*dp_on(delta))/3
    return d 

get_max_bond_length = {
    'pl-as': distance_pl_as,
    'pl-s': distance_pl_s}

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


def read_df(filen, Csize, sdelta):
    df = pd.read_pickle(filen)
    df.bend = (360*df.bend/(2*np.pi))
    df.bend_parallel = (360*df.bend_parallel/(2*np.pi))
    df.bend_non_parallel = (360*df.bend_non_parallel/(2*np.pi))
    df = df[df.cluster_type == 'chain']
    df = df[(df.cluster_size >= Csize) & (df.cluster_size<=Csize+sdelta)]
    df = df.dropna(subset=['bond_vec'])
    df['delta'] = pd.to_numeric(df['delta'])
    df['average_cos'] = df.apply(lambda row: get_bond_angles(row['bond_vec'], row['sequence'], row['delta'], Csize), axis=1)
    return df 


#Attention: 'p' and 'np' default to 'p-off-sb' because of the old versoin with which
# pl-as was run. 
t_bond = {
     'p': lambda delta: small_d(delta)[:2] if delta<=0.5 else big_d(delta)[:2],
     'np': lambda delta: small_d(delta)[:2] if delta<=0.5 else big_d(delta)[:2],
     'p-on': lambda delta: np.array([1,0]),
     'p-off-bs': lambda delta: big_d(delta)[:2] if delta<=0.5 else small_d(delta)[:2],
     'p-off-sb': lambda delta: small_d(delta)[:2] if delta<=0.5 else big_d(delta)[:2]
}


def get_bond_angles(bond_vec, seq, delta,Csize):
    L = (Csize-1)  - 1 
    bond_angle = np.zeros(L)

    cos_ba = np.zeros(L)
    Denom = np.zeros(L)
    for pos_i in range(0,L):
        for step_j in range(0,L-pos_i):


            norm_a = np.linalg.norm(bond_vec[pos_i])
            norm_b = np.linalg.norm(bond_vec[pos_i+step_j+1])
            cos_i = (np.dot(bond_vec[pos_i], bond_vec[pos_i+step_j+1]))/(norm_a*norm_b)
            arccos_i = np.arccos(cos_i)

            #print( t_bond[seq[pos_i]](delta), t_bond[seq[pos_i+step_j+1]](delta))
            arccos_compare = np.arccos(np.dot(t_bond[seq[pos_i]](delta), t_bond[seq[pos_i+step_j+1]](delta)))


            if np.isnan(arccos_compare):
                 arccos_compare = 0
            #if(delta == 0.5):
            #     print(t_bond[seq[pos_i]](delta), arccos_compare)

            cos_ii = np.cos(arccos_i - arccos_compare)
            cos_ba[step_j] +=  cos_ii  
            Denom[step_j] += 1

    # average cosine per Length L 
    cos_ba = cos_ba/Denom
    return cos_ba


#################

parser = argparse.ArgumentParser()
parser.add_argument("-size", type=int, default=10)
parser.add_argument("-srange", type=int, default=5)
args = parser.parse_args()

f_as = "/Users/ada/Documents/Code_Development_2016/2D_patchy/chain_analysis/small_data_set/rr_chains/lege_batterie/asymm/mu_0.3/chain_bond_angle.csv"
f_s = "/Users/ada/Documents/Code_Development_2016/2D_patchy/chain_analysis/small_data_set/rr_chains/lege_batterie/symm/7.2_8.2/chain_bond_angle.csv"

Csize=args.size
sdelta=args.srange
df_symm = read_df(f_s, Csize, sdelta)
df_asymm = read_df(f_as, Csize, sdelta)
df_topo = {'symm': df_symm, 'asymm': df_asymm}

a=1/3
b=1/3
c=1/3

inset = {'symm': [0.22,0.6, 0.2, 0.2],
         'asymm':[0.6,0.6, 0.2, 0.2] }


for energy in ['8.2']:

    ncolors = 3
    cmap=plt.cm.brg
    xlabel="$\Delta$"
    ylabel="$l_{p}/l_{c}$"
    fig, ax  = basic_plot_config(xlabel, ylabel, " ", ncolors,cmap)
    ax2 = []
    for topology in ['symm', 'asymm']:
         lp_lc = []
         lp_lc_err = []
         df = df_topo[topology]

         left, bottom, width, height = inset[topology]
         ax2.append(fig.add_axes([left, bottom, width, height]))

         
         ncolors = 7
         cmap=plt.cm.rainbow
         cy = cycler('color', cmap(np.linspace(0,1,ncolors)) )
         ax2[-1].set_prop_cycle(cy)
         ax2[-1].set_xlabel("$l_{d}$", size=20)
         ax2[-1].set_ylabel("$\\langle cos(\\theta)\\rangle$", size=20)
         ax2[-1].set_title("bond angle decorrelation for {}".format(name_dict[topology]))
         ax2[-1].set_ylim([0.9,1.0])


         for delta in np.round(np.arange(0.2,0.9,0.1),2):
              dg = df[(df.energy == energy) & (df.delta == delta)]
              if name_dict[topology] == 'pl-s':
                   seq_arr = df.sequence.values
                   seq_arr = np.concatenate(seq_arr).ravel()
                   unique, counts = np.unique(seq_arr, return_counts=True)
                   # for all pl-off-sb
                   N_all = np.sum(counts)
                   a= counts[unique=='p-off-sb']/(N_all)
                   b= counts[unique=='p-off-bs']/(N_all)
                   c= counts[unique=='p-on']/(N_all)

              arr = dg.average_cos.values
              arr = np.array([ np.array(item) for item in arr])
              mean_arr = np.mean(arr, axis=0)
              sem = stats.sem(arr, axis=0, ddof=1, nan_policy='propagate')

              ax2[-1].errorbar(np.arange(1,Csize-1),mean_arr,
                               sem, label=delta,
                               capsize=4,
                               capthick=1,
                               lw=2,
                               ms=10)
              mean_arr = -np.log(mean_arr)
              sem = stats.sem(arr, axis=0, ddof=1, nan_policy='propagate')
              #lp = np.polynomial.polynomial.polyfit(np.arange(1,Csize-1), mean_arr, 1, rcond=None, w=1/sem)[1]

              def line(x, a,b):
                   return x*a + b

              popt, pcov = curve_fit(line, np.arange(1,Csize-1), mean_arr, sigma=sem)
              perr = np.sqrt(np.diag(pcov))[0]

              lp = 1/popt[0]
              lp_err= lp * (perr/popt[0]) 
              sigma_contour = get_max_bond_length[name_dict[topology]](delta, a,b,c)
              sigma_contour = 1
              lp_lc.append((lp/sigma_contour)/(Csize*sigma_contour))
              lp_lc_err.append((lp_err/sigma_contour)/(Csize*sigma_contour))
              
              ax2[0].legend(bbox_to_anchor=(1.2,0.75), ncol=1, loc=2, borderaxespad=0, fontsize=10)

         ax.errorbar(np.arange(0.2,0.9,0.1), lp_lc, lp_lc_err, lw=3, marker='o', ms=5, label=name_dict[topology])
         ax.legend(bbox_to_anchor=(0.4,0.2), ncol=1, loc=2, borderaxespad=0, fontsize=20)
         ax.set_ylim([0,120])
    plt.savefig("persistence_{}.pdf".format(energy))



