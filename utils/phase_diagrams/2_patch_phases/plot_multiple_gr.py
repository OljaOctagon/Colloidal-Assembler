import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib as mpl
import seaborn as sns
from matplotlib import rc 
import glob
import re 

sns.set(style='white')
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

Delta=0.01
R=10
r=(np.arange(Delta,R,Delta)/2.)    

def set_basic_plot(ax):
    ax.plot(r,np.ones(len(r)), 'k--', lw=3, alpha=0.5)
    ax.set_xlabel("$ \sigma $", size=30)
    ax.tick_params(labelsize=22)
    ax.set_ylabel("$g( \sigma )$", size=30)
    ax.locator_params(nbins=6)
    ax.legend(prop={'size':18})
    ax.set_xlim([0,5])
    ax.set_ylim([0,6.5])
    return ax 

filelist = glob.glob("*_gr2d.dat")
fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1,figsize=(12,6))

ax1=set_basic_plot(ax1)
ax2=set_basic_plot(ax2)

colors={
    '0.2':'#6300ff',
    '0.3':'#00f2a3',
    '0.4':'#ff5e62'}
z=0
alpha=1.0
for ffile in filelist:
    numbers = re.findall(r"[-+]?\d*\.\d+|\d+", ffile)
    mu = numbers[0]
    energy = numbers[1]
    delta = numbers[2]
    pressure=numbers[3]
    
    print(pressure, delta)

    if pressure == '60':
        ax = ax1
    elif pressure == '100':
        ax = ax2
    else:
        ax=None

    if ax: 
        G = np.loadtxt(ffile)
        z-=1
        ax.plot(r,G,c=colors[delta], lw=2.5, alpha=alpha, label="open-boxes,  $\Delta = {}$, $P={}$".format(delta,pressure),
                zorder=z)

    plt.tight_layout()

ax1.legend(loc='upper right', fontsize=13)
ax2.legend(loc='upper right', fontsize=13)
plt.savefig("gr_open_boxes.pdf")
plt.show() 
