import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler


df = pd.read_csv("pl_chains.csv")
def basic_plot_config(xlabel, ylabel, title, ncolors,cmap):
    fig, ax = plt.subplots(figsize=(12,10))

    plt.tick_params(axis='both', labelsize=18) 
    plt.title(title, size=20)
    plt.xlabel(xlabel, size=25)
    plt.ylabel(ylabel, size=25)
    cmap=cmap
    c = cycler('color', cmap(np.linspace(0,1,ncolors)) )
    ax.set_prop_cycle(c)
    return fig,ax

xlabel = '$\epsilon [k_B T]$'
ylabel = '$\\langle L \\rangle$'
title = 'Average chain length $\\langle L \\rangle$ by interaction strength $\epsilon$'
ncolors = 6 
cmap = plt.cm.jet 

if __name__ == '__main__':

    fig,ax = basic_plot_config(xlabel, ylabel, title, ncolors,cmap)

    names = {'symm':'pl-s', 'Asymm':'pl-as', 'center':'pl-center'}
    z=0
    df2 = df.groupby(['delta','energy','topology']).cluster_size.agg(['mean','std']).reset_index()

    for delta in [0.2,0.5,0.8]:
        for topology in ['symm', 'Asymm']:


            if delta==0.5 and  topology=='Asymm':
                next(ax._get_lines.prop_cycler)
            

            if delta!=0.5 or topology!='Asymm':
                arr = df2[(
                    df2.topology ==topology)&(df2.delta == delta)][['energy','mean','std']].values

                if delta==0.5:
                    topology = 'center'

                plt.errorbar(arr[:,0], arr[:,1], arr[:,2],
                             marker='o', ms=10,
                             capsize=5, capthick=5, lw=2,
                             label='{}, $\Delta={}$'.format(names[topology],delta), zorder =z )
                z-=1

    plt.yticks(np.arange(0,36,3))
    plt.legend(bbox_to_anchor=(0.02,0.95), loc=2, ncol=3, fontsize=18)
    plt.xlim([5,8.5])
    plt.ylim([0,36])
    plt.xticks([5.2,6.2,7.2,8.2],[-5.2,-6.2,-7.2,-8.2])
    plt.tight_layout()


    a = plt.axes([.20, .30, .40, .40])
    a.spines["top"].set_visible(False)    
    a.spines["bottom"].set_visible(False)    
    a.spines["right"].set_visible(False)    
    a.spines["left"].set_visible(False)    
    a.get_xaxis().tick_bottom()    
    a.get_yaxis().tick_left()    

    plt.tick_params(axis='both', labelsize=20) 
    plt.title('Cluster size distribution for pl-s $\epsilon = -7.2 k_{B}T$', size=20)
    plt.xlabel('L', size=25)
    plt.ylabel('P', size=25)
    plt.xlim([1.5,15])
    cmap=plt.cm.jet
    c = cycler('color', cmap(np.linspace(0,1,6)) )
    a.set_prop_cycle(c)

    pattern = ['|','\\',"//"]
    zorder = [2,2,1,1,0,0]
    i=0

    for delta in [0.2,0.5,0.8]:
        for topology in ['symm']:
            for energy in [7.2]:
                df[(df.topology ==topology)&(
                    df.energy == energy)&(
                    df.delta == delta)].cluster_size.plot(kind='hist',
                                                     label=delta,
                                                     alpha=0.7,
                                                     density=True,
                                                     hatch=pattern[i])
                i+=1
                next(a._get_lines.prop_cycler)

    plt.tight_layout()
    plt.savefig("chain_length_rr.pdf")
    plt.show()
