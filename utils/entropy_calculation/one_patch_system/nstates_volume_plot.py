import numpy as np
import pandas as pd 
from matplotlib import rc 

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

arr = pd.read_csv("nstates.dat", header=None, delim_whitespace=True,
	names=['patch','N_nonzero', 'N_total', 'V_frac', 'V']).values

import matplotlib.pyplot as plt 
fig, ax = plt.subplots()
plt.plot(arr[:,0], arr[:,4], marker='o', alpha=0.6, c='b')
plt.xlabel("$\displaystyle \Delta$", size=30)
plt.ylabel("volume of states", size=18)
ax.tick_params(size=12, labelsize=14)
plt.tight_layout()
plt.show()