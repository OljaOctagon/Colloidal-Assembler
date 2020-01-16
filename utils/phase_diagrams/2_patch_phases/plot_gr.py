import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib as mpl
import seaborn as sns
from matplotlib import rc 
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("-f", default='gr.dat', help="input file")

args = parser.parse_args()


sns.set(style='white')
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

G = np.loadtxt(args.f)
Delta=0.01
R=10
r=(np.arange(Delta,R,Delta)/2.)    
plt.plot(r,np.ones(len(r)), 'k--', lw=3, alpha=0.5)
plt.plot(r,G,c="#873e84", lw=3, label="open-boxes,  $\Delta = 0.3$ , $P=100$")
plt.xlabel("$ \sigma $", size=30)
plt.tick_params(labelsize=22)
plt.ylabel("$g( \sigma )$", size=30)
plt.locator_params(nbins=6)
plt.legend(prop={'size':18})
plt.xlim([0,5])
plt.ylim([0,5])

plt.tight_layout()
plt.savefig("gr_cluster_mu_0.3Asymm_patchpos_0.3.pdf")
plt.show()      
