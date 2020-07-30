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


style.use('seaborn-poster') 
mpl.rcParams['font.family'] = "sans-serif"
sns.set_context('poster')

def read_bonds(filen):
    
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

if __name__ == '__main__':

    # get all check point values and sort them
    checkpoints= glob.glob("Box*.bin")
    check_point_values = np.sort(
    [ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])
    # import the bonds file:
    # file format:
    # ------------------------
    # #new time
    # particle_id1 particle_id2  patch_id1 patch_id2]
    #  ....
    # --------------------------

    # network_arr format: network_arr.shape = ( frame_i, bond_rows_frame_i )
    network_arr = read_bonds("patch_network.dat")
    # patch position calculation

    print("test:", len(network_arr), len(check_point_values))

    for j, arr in enumerate(network_arr[-1:]):
        
        arr_id = arr[:,:2]
        
        print("arr", arr_id.shape)
        arr_id_sorted = np.sort(arr_id)

        uniques, counts = np.unique(arr_id_sorted, return_counts=True, axis=0)
        
        av_counts = np.mean(counts)
        max_counts = np.max(counts)

        n_double = np.sum(counts>1)
        percent_n_double = n_double/len(counts)
        
        print("av counts", av_counts)
        print("max_counts", max_counts)
        print("n_double", n_double)
        print("percent_n_double", percent_n_double)

