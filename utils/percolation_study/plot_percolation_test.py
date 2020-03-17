import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.style as style
style.use('seaborn-poster') 
mpl.rcParams['font.family'] = "sans-serif"
sns.set_context('poster')
plt.rcParams['axes.axisbelow'] = True

# dma-as1
file='percolation_test.dat'
dir='/Users/ada/Documents/Code_Development_2020/rhombi/percolations_study/dma_as1/'

df = pd.read_csv(dir+file)

pdict={'no':'r','yes':'b','test':'y'}

arr = df[['phi','T']].values
brr=df.percolation.values
state = [ pdict[item] for item in brr]
plt.title("Percolation line for dma-as1", fontsize=25)
plt.grid()
plt.xlabel("$\phi$", size=25)
plt.ylabel("$k_{b}T/\epsilon$", size=25)
plt.scatter(arr[:,0], arr[:,1], c=state, s=80)
plt.tight_layout()
plt.savefig("{}percolation_dma-as1.pdf".format(dir))
