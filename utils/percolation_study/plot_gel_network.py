import numpy as np
import pandas as pd
import argparse
import gel_tools as gt
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.style as style

blue_c ='#9999FF'
red_c ='#FF9999'
purple_c='#C17DCB'

def get_pdist(vpos):
    l = len(vpos)
    pdist = np.zeros((l,l,2))

    for i in range(2):
        p = np.reshape(vpos[:,i], (l,1))
        pdist = p - p.transpose()

    N_pdist = np.sqrt(
        np.power( pdist[:,:,0], 2)
        + np.power( pdist[:,:,1], 2))

    return pdist, N_pdist

def is_distance_smaller(pos_i, id1,id2,cutoff):
    dist = pos_i[id1] - pos_i[id2]
    Ndist = np.linalg.norm(dist)
    if Ndist < cutoff:
        return True
    else:
        return False

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, default='./')
    parser.add_argument('-f', type=str, default='patch_network.dat')
    args = parser.parse_args()

    N_particles = 1500 
    valence=4
    T=0.15
    phi=0.5
    run_id="dma-as1_r001"

    # 1: paralelel
    # 0: non parallel
    bond_dict = {(0,0): 1, (0,1): 0, (1,1):1,
                (2,2): 1, (2,3): 0, (3,3):1,
                (0,2): 0, (0,3): 1, (1,2):1,
                (1,3):0}


    # read in connections, format: particle id1, particle id2, patch id1, patch id2
    dir_name = args.d
    f_name = args.f

    # get all check point values and sort them
    checkpoints= glob.glob("{}Box*.bin".format(dir_name))
    check_point_values = np.sort(
    [ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])

    pn_file = dir_name+f_name

    # network_arr format: network_arr.shape = ( frame_i, bond_rows_frame_i )
    connections = gt.read_bonds(pn_file)

    results_dir = 'graph_plots_colored'
    # make frame directory if it doesn't exist 
    if not os.path.isdir(dir_name+results_dir):
        os.mkdir(dir_name+results_dir)

    for j,val in enumerate(check_point_values[1:]):

         print("Checkpoint ....",val)
         # Initialize new results

         connections_j = connections[j]
         pos_j = np.fromfile("{}positions_{}.bin".format(dir_name, val))
         pos_j = np.reshape(pos_j, (-1,3))
         pos_j = pos_j[:,:2]

         fig,ax = plt.subplots()
         ax.set_title("Frame {}".format(j))
         plt.xlim((-1,50))
         plt.ylim((-1,60))
         plt.axis("equal")
         plt.axis('off')

         for connection in connections_j:
             id1=connection[0]
             id2=connection[1]
             if is_distance_smaller(pos_j, id1,id2,5):
                 X = point = np.array([pos_j[id1],pos_j[id2]])

                 plt.plot(X[:,0],X[:,1], lw=1, c=purple_c, alpha=1)

                 id_list =[1311,1355]
                 if id1 in id_list and id2 in id_list:
                     plt.plot(X[:,0],X[:,1], lw=1, c=red_c, alpha=1)

         ax.scatter(pos_j[:,0], pos_j[:,1], s=0.5, c='k', alpha=0.7)

         plt.savefig("{}{}/frame_{}.png".format(dir_name,results_dir,val), dpi=500)
         plt.close()


# Calculate distances

# plot gel network: connect center to center points, plot rhombi in the background maybe
# color according to p/np bond
# color according to degree (branch like, chain like, crystal like )
# color all 3p loops
# simplify graph 3p loops are one element
# -- calculate which p are members of 3 p loops
# -- replace these partilces ids with virutal partilces N>1600
# --- degree of virtual loop particles 
