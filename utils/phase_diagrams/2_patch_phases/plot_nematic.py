import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from cycler import cycler


df = pd.read_csv("nematic_op.csv")
df_as =  pd.read_csv("nematic_op_as.csv")


def plot_observable(ax, df, df_as, ylabel, parameter, name):

    arr_local = df.groupby(
            ['delta'])[parameter].agg(
                    ['mean','std']).reset_index().values


    brr_local = df_as.groupby(
            ['delta'])[parameter].agg(
                    ['mean','std']).reset_index().values

    arr_local[3] = brr_local[3]

    ax.errorbar(arr_local[:,0],
                arr_local[:,1],
                arr_local[:,2],
                label='pl-s', marker='o',
                ms=7,
                capsize=5,
                capthick=5,
                lw=1)

    ax.errorbar(brr_local[:,0],
                brr_local[:,1],
                brr_local[:,2],
                label='pl-as', marker='o',
                ms=7,
                capsize=5,
                capthick=5,
                lw=1)

    ax.errorbar(brr_local[3,0],
                brr_local[3,1],
                brr_local[3,2],
                label='pl-center', marker='o',
                ms=7,
                capsize=5,
                capthick=5,
                lw=1)
    ax.set_title(name, size=18)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_xlabel("$\Delta$", size=18)
    ax.set_ylabel(ylabel, size=18)
    ax.legend(loc='best')
    ax.set_ylim([0,1])

fig,(ax2,ax1)  = plt.subplots(nrows=1,ncols=2, figsize=(13,7))
cmap=plt.cm.Dark2
c = cycler('color', cmap(np.linspace(0,1,4)) )
ax1.set_prop_cycle(c)
ylabel = "$S_{largest}$"

plot_observable(ax1, df, df_as, ylabel, 'local_nematic_op', 'Nematic order of largest cluster')

cmap=plt.cm.jet
c = cycler('color', cmap(np.linspace(0,1,4)) )
ax2.set_prop_cycle(c)
ylabel = "$f_{largest}$"
plot_observable(ax2, df, df_as, ylabel, 'fraction_largest_nematic', 'Fraction of larget cluster')

plt.tight_layout()
plt.savefig("{}.pdf".format("pl_nematic"), dpi=300)

