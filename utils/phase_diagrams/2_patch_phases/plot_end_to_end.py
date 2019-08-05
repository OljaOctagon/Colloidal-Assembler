import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from cycler import cycler

df = pd.read_csv('chains_end_to_end.csv')


df.groupby(['energy','delta','cluster_size']).end_to_end_distance.max().reset_index()


df['max'] = df.groupby(['energy','delta','cluster_size']).end_to_end_distance.transform('max')
df['scaled_e_to_e'] = df.end_to_end_distance/df['max']

dg = df.groupby(['energy','delta','cluster_size']).scaled_e_to_e.agg(['mean', 'std']).reset_index() 


for energy in np.arange(5.2,9.2,1):
    fig,ax = plt.subplots()
    cmap=plt.cm.rainbow
    c = cycler('color', cmap(np.linspace(0,1,8)) )
    ax.set_prop_cycle(c)

    for delta in np.arange(0.2,0.9,0.1):
        delta = np.round(delta,1)
        dg2 = dg[(dg.delta == delta) & (dg.energy==energy)]
        plt.errorbar(dg2.cluster_size.values,
                     dg2['mean'].values,
                     dg2['std'].values,
                     capsize=4,
                     capthick=1,
                     lw=2,
                     ms=10, 
                     label="{}".format(delta))


    plt.xlabel("$L$", size=20)
    plt.ylabel("$\epsilon$", size=20)
    plt.xlim([0,20])
    plt.ylim([0.65,1.05])
    plt.legend(loc='best')

    plt.savefig("{}_end_to_end.pdf".format(energy))
