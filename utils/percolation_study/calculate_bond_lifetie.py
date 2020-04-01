
import numpy as np
import pandas as pd
import argparse
import gel_tools as gt
import networkx as nx
import glob
import os

from collections import defaultdict 

if __name__=='__main__':

    N_particles = 1500 

    # get all check point values and sort them
    checkpoints= glob.glob("{}Box*.bin".format(dir_name))
    check_point_values = np.sort(
    [ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])

    pn_file = "patch_network.dat" 
    # network_arr format: network_arr.shape = ( frame_i, bond_rows_frame_i )
    connections = gt.read_bonds(pn_file)
    Ltime = len(connections)
    connection_soup = np.array([ ci for ci in connections[time_j] for time_j in range(Ltime)])
    connection_soup = connection_soup[:2]
    sorted_connection_soup = np.sort(connection_soup)
    uniques = np.unique(sorted_connection_soup, axis=0)

    bond_sequences = defaultdict(np.zeros(Ltime))
    for bond in uniques:
        i = bond[0]
        j = bond[1]

        for time in range(Ltime):
            if [i,j] in connections[time,:2]
            bond_sequences[(i,j)][time] = 1

    print(bond_sequences)
