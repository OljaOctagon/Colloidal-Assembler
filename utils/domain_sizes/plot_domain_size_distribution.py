import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd 
sns.set(style='white')

file_id=['feq-center',
       'dmo-center',
       'dma-center',
       'dma-as2']

fig, axes = plt.subplots(2,2)

from itertools import product
Coord = list(product(np.arange(2),np.arange(2)))

for fid, coord in zip(file_id, Coord):
    ax = axes[coord[0], coord[1]]

    arr_p = pd.read_csv(fid+'-p.dat', delim_whitespace=True, header=None).values
    arr_np = pd.read_csv(fid+'-np.dat', delim_whitespace=True, header=None).values

    nbins_p = int(np.ceil(np.max(arr_p)/2.))
    nbins_np = int(np.ceil(np.max(arr_np)/2.))
    ax.set_title(fid)
    ax.hist(arr_p, bins=nbins_p, alpha=0.4, normed=True, facecolor='red', label='parallel')
    ax.hist(arr_np, bins=nbins_np, alpha=0.4, normed=True, facecolor='blue', label='non-parallel')
    ax.set_xlabel("domain sizes")
    ax.set_ylabel("P")
    handles, labels = ax.get_legend_handles_labels()

fig.legend(handles, labels, loc=(0.35,0), ncol=2)
plt.tight_layout()
plt.savefig("domain_size_histogram.png", dpi=300)
