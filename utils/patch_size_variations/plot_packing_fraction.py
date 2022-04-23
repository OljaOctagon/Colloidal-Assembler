import numpy as np
import matplotlib.pyplot as plt
import glob 
delta=[0.2,0.3,0.4,0.5]
radius=[0.05,0.1,0.15,0.2]

style.use('seaborn-ticks')
mpl.rcParams['font.family'] = "sans-serif"
#sns.set_context('poster')
plt.rcParams['axes.axisbelow'] = True

plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 8

#mu0.14_d0.35_patchr0.2_x
arr_packing = []
fig,ax = plt.subplots(ncols=len(delta),nrows=len(radius),figsize=(20,22))
lines=['-','--','.-','.']*4 

for d in delta:
    for r in radius:
        packing_list = glob.glob("mu0.14_d{}_patchr{}_*/NPT_OUT.txt".format(d,r))
        print(packing_list)
        ax[d,r].set_xlabel("MC sweeps")
        ax[d,r].set_ylabel("packing fraction")
        for i, pfile in enumerate(packing_list):
            sns.set_palette(sns.cubehelix_palette(16, start=0.5, rot=-0.75))
            df = pd.read_csv(pfile, names=['time','a','b','c','phi'], delim_whitespace="True")
            arr = df[['time','phi']].values
            ax[d,r].plot(arr[:,0],arr[:,1],lw=1,ls=lines[i])


fig.savefig("packing_fraction_delta_r.pdf")
