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
style.use('seaborn-poster') 
mpl.rcParams['font.family'] = "sans-serif"
sns.set_context('poster')

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

    for j,val in enumerate(check_point_values[-10:]):

        network_arr[j]
        uniques, counts = np.unique(network_arr[j], return_counts=True, axis=1)

        av_counts = np.mean(counts)
        max_counts = np.max(counts)

        n_double = np.sum(counts>2)
        percent_n_double = n_double/len(counts)

        print(av_counts, max_counts, n_double, percent_n_double)




