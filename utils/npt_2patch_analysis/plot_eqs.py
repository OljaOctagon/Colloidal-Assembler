import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd

import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter


sns.set(style='white')
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


#X_03 = pd.read_csv("eqs_0.2.dat", header=None, delim_whitespace=True).values
X_05 = pd.read_csv("eqs_0.3.dat", header=None, delim_whitespace=True).values
#X_07 = pd.read_csv("eqs_0.4.dat", header=None, delim_whitespace=True).values


''' 
bg_arr = np.loadtxt("domains.dat")

f_arr = np.reshape(np.array([0,0,0]), (1,3))
bg_arr = np.concatenate((f_arr,bg_arr), axis=0)
l_bg = np.reshape(bg_arr[-1], (1,3))
bg_arr = np.concatenate((bg_arr,l_bg), axis=0)

fig, ax = plt.subplots()

max_frac = np.max(bg_arr[:,1])
X=np.concatenate((f_arr, X_03), axis=0)
l_arr = np.reshape(np.array([0,0.8,0]),(1,3)) 
X = np.concatenate((X,l_arr),axis=0)

print(X.shape, bg_arr.shape)
for i in range(len(bg_arr)-1):
    ax.axvspan(X[i,1], X[i+1,1], color=plt.cm.magma_r(bg_arr[i+1,1]/max_frac), alpha=0.5)
'''
fig, ax = plt.subplots()
#plt.errorbar(X_03[:,1], X_03[:,0], xerr=X_03[:,2],
#             lw=0, ms=8, marker='o', c='#ff5e62', alpha=0.7, label="$\Delta=0.2$")

plt.errorbar(X_05[:,1], X_05[:,0], xerr=X_05[:,2],
             lw=0, ms=8, marker='s', c='#00f2a3', alpha=0.7,label="$\Delta=0.3$")

#plt.errorbar(X_07[:,1], X_07[:,0], xerr=X_07[:,2],
#             lw=0, ms=8, marker='^', c='#6300ff', alpha=0.7, label="$\Delta=0.4$")

plt.legend(loc=2,prop={'size':14})


'''
ax.axvspan(0.35,0.60, color='k', alpha=0.1)
ax.axvspan(0.60,0.70, color='k', alpha=0.3)
ax.axvspan(0.70,0.8, color='k', alpha=0.5)
ax.text(0.45, 90, "liquid", fontsize=16)
ax.text(0.61, 90, "hexatic", fontsize=16)
ax.text(0.71, 90, "solid", fontsize=16)
'''

plt.ylim([0,110])
#plt.xlim([0.35,0.7])
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlabel("$\phi$", size=30)
plt.ylabel("$P$", size=25)
plt.tight_layout()
plt.savefig("eqs_stars.png", dpi=300)

