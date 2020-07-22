import numpy as np
import numpy as np
import matplotlib.pyplot as plt 
import glob 

systems=['dma_as1','dmo_as1','dmo_s1','dmo_s2']
deltas=[0.2,0.5]

sys_dict = {"dma_as1":"double_manta_asymm_1",
            "dmo_as1":"double_mouse_asymm_1",
            "dmo_s1":"double_mouse_symm_1",
            "dmo_s2":"double_mouse_symm_2"}

fig,ax = plt.subplots(ncols=len(deltas),nrows=len(systems),figsize=(12,12))

for si, system in enumerate(systems):
    for di, delta in enumerate(deltas):
        ndir="{}/{}_phi_*_delta_{}_temp_*/psi_op.dat".format(system,sys_dict[system],delta)
        files = glob.glob(ndir)

        mean_psi = []
        std_psi = []

        mean_N = []
        std_N = []

        for file_i in files:

            phi_i = file_i.split("_")[2]
            temp_i = file_i.split("_")[6]  

            arr=np.loadtxt(file_i)
            arr = arr[-10:]

            mean_psi.append([phi_i, temp_i, np.mean(arr[:,1])])
            mean_psi.append([phi_i, temp_i, np.mean(arr[:,1])])
            mean_N.append([phi_i, temp_i, np.mean(arr[:,2])])
            std_N.append([phi_i, temp_i, np.mean(arr[:,2])])
            img = ax[si,di].scatter(mean_N[:,0], mean_N[:,1], c=mean_N[:,2], cmap='viridis')

            ax[si,di].set_xlabel("$\phi", size=25)
            ax[si,di].set_ylabel("T ")

    cbar = fig.colorbar(img, ax=ax,
                orientation='horizontal',
                fraction=0.05,
                boundaries=np.linspace(0,1500,20), ticks=[0,50,100,200,500,750,1500])
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('$\psi$', size=12)


#cmap = plt.cm.viridis
#norm = matplotlib.colors.Normalize(vmin=0, vmax=1500)
#sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#fig.colorbar(sm)
