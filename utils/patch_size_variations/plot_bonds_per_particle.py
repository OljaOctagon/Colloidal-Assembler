import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib as mpl
import os
import glob
import networkx as nx
import argparse
import seaborn as sns
import matplotlib.style as style
import pandas as pd 
from matplotlib.colors import LinearSegmentedColormap

style.use('seaborn-poster') 
mpl.rcParams['font.family'] = "sans-serif"
sns.set_context('poster')

blue_c ='#9999FF'
red_c ='#FF9999'
purple_c='#C17DCB'
green_c='#61DCB7'



def read_bonds_depricated(filen):
    
    connection_soup = pd.read_csv(filen, delim_whitespace=True,header=None).values 

    delta_id = connection_soup[:-1,0] - connection_soup[1:,0]
    cut_pos = np.where(delta_id>100)[0]  
    
    cut_pos = np.insert(cut_pos,0,-1)
    cut_pos = np.insert(cut_pos,-1,len(connection_soup)-1)
    network_list = []
    
    for i, pos_i in enumerate(cut_pos[:-1]):

        start=cut_pos[i]+1
        end = cut_pos[i+1]
        if( (end-start) > 100):
            arr=connection_soup[start:end+1,:]    
            network_list.append(arr)

    return network_list

def read_bonds(filen):
	first_line_pair = [0,0,0,0]
	cut=False
	with open(filen, 'r') as f:
		network_list = []
		for line in f:
			if "#" in line:
				network_list.append([])
				first_line_pair = [0,0,0,0]
				cut=False

			else:
				line_counter=len(network_list[-1])
				pairs = list(map(int, line.split(" ")))
				if pairs == first_line_pair or cut==True:
					cut=True
				else:
					network_list[-1].append(np.array(pairs))

				if line_counter == 0:
					first_line_pair = pairs
	network_list = [ np.array(item) for item in network_list]

	return network_list



if __name__ == '__main__':

    #bonding_volume_delta_0.4_radius_0.18
    all_dirs = glob.glob("bonding_volume_delta*")

    double_bonds =  [] 

    for filedir in all_dirs:
        network_arr = read_bonds("{}/patch_network.dat".format(filedir))

        delta=filedir.split("_")[3]
        radius=filedir.split("_")[5]

        for j, arr in enumerate(network_arr[-1:]):

            arr_id = arr[:,:2]
            arr_id_sorted = np.sort(arr_id)

            uniques, counts = np.unique(arr_id_sorted, return_counts=True, axis=0)
            av_counts = np.mean(counts)
            max_counts = np.max(counts)

            n_double = np.sum(counts>1)
            percent_n_double = n_double/len(counts)

        double_bonds.append([delta,radius,max_counts])

    double_bonds = np.array(double_bonds).astype(float)
    arr= double_bonds[double_bonds[:,2]>1]

    unique_delta_sorted = np.sort(np.unique(arr[:,0]))

    db_curve=[]
    for delta in unique_delta_sorted:
        min_r=np.min(arr[arr[:,0]==delta][:,1])
        db_curve.append([delta,min_r])

    db_curve = np.array(db_curve)
    plt.plot(db_curve[:,0],db_curve[:,1],marker='o',markersize=12, linestyle="--",lw=3,c='#23bd92')
    # want a tuple of postions (delta1,radius1, delta2, radius2, ...) 

    #plot this curve
    plt.xlabel("$\Delta$")
    plt.ylabel("r")
    plt.savefig("curve_double_bond_per_particles.pdf")

    import pandas as pd

    
    #cmap = plt.cm.GnBu_d  # define the colormap
    # extract all colors from the .jet map

    cmap = sns.cubehelix_palette(8, start=.5, rot=-0.75, as_cmap=True)

    bounds = np.linspace(0.5,3.5,4)
    print(bounds)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    fig,ax=plt.subplots(figsize=(13,10))
    df_bonds = pd.DataFrame(data=double_bonds, columns=["delta","r","nbonds"])
    df_bonds = df_bonds.pivot("delta", "r", "nbonds")
    ax = sns.heatmap(df_bonds,annot=True, cmap=cmap,
                     norm=norm, linewidth=0.5, linecolor='k',cbar=True, 
                     cbar_kws={"drawedges":True, "ticks":[1,2,3]})


    plt.xlabel("$\Delta$")
    plt.ylabel("r")

    plt.tight_layout()
    plt.savefig("poly_bond_per_particles.pdf")
