import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from cycler import cycler


df = pd.read_csv("yield.csv", header=None)
df.columns=['index', 'mu','energy', 'topology', 'delta', 'cluster_type', 'mean', 'std']

fig, ax = plt.subplots()

# get colormap
cmap=plt.cm.gist_rainbow_r
c = cycler('color', cmap(np.linspace(0,1,7)) )
ax.set_prop_cycle(c)

plt.xlim([0.185,0.815])
plt.ylim([0,1])
plt.tick_params(axis='both', labelsize=16)
plt.xticks([0.2,0.5,0.8])
plt.xlabel("$\Delta$", size=25)
plt.ylabel("yield [%]", size=20)
for ctype in df.cluster_type.unique():
    dg = df[ df.cluster_type == ctype ]
    arr = dg[['delta','mean','std']].values
    plt.errorbar(arr[:,0],arr[:,1], yerr=arr[:,2], label=ctype,
                 capthick=3, capsize=3,lw=3,)
plt.legend(bbox_to_anchor=(1.2,0.75), ncol=1, loc=2, borderaxespad=0, fontsize=14)
plt.tight_layout()
plt.savefig("yield_all_clusters.pdf")
plt.show()
