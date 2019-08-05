import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from cycler import cycler


df = pd.read_csv("nematic_op.csv")
df_as =  pd.read_csv("nematic_op_as.csv")


def plot_observable(df, df_as, ylabel, parameter, name):

    arr_local = df.groupby(
            ['delta'])[parameter].agg(
                    ['mean','std']).reset_index().values


    brr_local = df_as.groupby(
            ['delta'])[parameter].agg(
                    ['mean','std']).reset_index().values

    arr_local[3] = brr_local[3]

    plt.errorbar(arr_local[:,0],
                arr_local[:,1],
                arr_local[:,2],
                label='pl-s', marker='o',
                ms=7,
                capsize=5,
                capthick=5,
                lw=1)

    plt.errorbar(brr_local[:,0],
                brr_local[:,1],
                brr_local[:,2],
                label='pl-as', marker='o',
                ms=7,
                capsize=5,
                capthick=5,
                lw=1)

    plt.errorbar(brr_local[3,0],
                brr_local[3,1],
                brr_local[3,2],
                label='pl-center', marker='o',
                ms=7,
                capsize=5,
                capthick=5,
                lw=1)

    plt.tick_params(axis='both', labelsize=20)
    plt.xlabel("$\Delta$", size=22)
    plt.ylabel(ylabel, size=22)
    plt.legend(loc='best')
    plt.ylim([0,1])
    plt.tight_layout()
    plt.savefig("{}.png".format(name), dpi=300)


fig,ax = plt.subplots()
cmap=plt.cm.Dark2
c = cycler('color', cmap(np.linspace(0,1,4)) )
ax.set_prop_cycle(c)
ylabel = "$S_{largest}$"

plot_observable(df, df_as, ylabel, 'local_nematic_op', 'local_nematic')


fig,ax = plt.subplots()
cmap=plt.cm.jet
c = cycler('color', cmap(np.linspace(0,1,4)) )
ax.set_prop_cycle(c)
ylabel = "$f_{largest}$"
plot_observable(df, df_as, ylabel, 'fraction_largest_nematic', 'fraction_largest_nematic')

