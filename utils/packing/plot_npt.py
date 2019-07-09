import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rc 
rc('text', usetex=True)

arr = pd.read_csv("NPT_OUT.txt", header=None, delim_whitespace=True).values

fig, ax = plt.subplots()
plt.plot(arr[:,0], arr[:,4], c='m', markeredgecolor='k', lw=2.5, marker='o')
plt.xlabel("MC time [sweeps]", size=20)
plt.ylabel(" $ \phi $", size=26)
plt.tight_layout()
plt.locator_params(axis='y', nbins=5)
ax.xaxis.get_offset_text().set_fontsize(24)
plt.tick_params(labelsize=15)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.savefig("npt.pdf")
plt.show()
