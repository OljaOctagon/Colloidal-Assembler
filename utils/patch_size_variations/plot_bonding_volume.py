import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
import matplotlib as mpl
import matplotlib.style as style

style.use('seaborn-ticks') 
mpl.rcParams['font.family'] = "sans-serif"
#sns.set_context('poster')
plt.rcParams['axes.axisbelow'] = True

plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['figure.titlesize'] = 15 

blue_c ='#9999FF'
red_c ='#FF9999'
purple_c='#C17DCB'
green_c='#61DCB7'

sns.set_palette(sns.cubehelix_palette(16, start=0.5, rot=-0.75))

df = df.read_csv("p_bonding_vol.dat",delim_whitespace=True, names=['delta','radius', 'volume', 'max_angle', 'min_angle'])

fig,ax = plt.subplots(figsize=(10,10))
deltas=np.sort(df.delta.values)
radii=np.sort(df.radius.values)
for radius in radii:
    arr= df[df.radius == radius][['delta','volume']].values
    plt.plot(arr[:,0],arr[:,1], lw=2,marker='o', ms=7, label="radius = {}".format(radius))

plt.legend(bbox_to_anchor=(0, -0.05), ncol=3)
ax.set_xlabel("radius")
ax.set_ylabel("V_{b}")

plt.show()
