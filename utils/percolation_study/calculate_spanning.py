import numpy as np
import pandas as pd
import argparse
import gel_tools as gt
import networkx as nx
import glob
import os

def calculate_pbc_images(pos_i,box_l,particles_max_domain):

    N_largest = len(particles_max_domain)
    # all periodic images of largest cluster
    N_images = 9
    N_virtual = N_largest*N_images
    virtual_pos = np.zeros((N_virtual,2))
    sign_array = np.array([[0,0],
                            [-1,0],
                            [-1,-1],
                            [0,-1],
                            [1,-1],
                            [1,0],
                            [1,1],
                            [0,1],
                            [-1,1]])

    for k,pid in enumerate(particles_max_domain):
        for image_j in range(N_images):
            virtual_pos[k+image_j*N_largest] = pos_i[pid] + sign_array[image_j]*box_l

    return virtual_pos

def get_pdist(vpos):
    l = len(vpos)
    pdist = np.zeros((l,l,2))
    print(l)
    for i in range(2):
        p = np.reshape(vpos[:,i], (l,1))
        p = p - p.transpose()

    N_pdist = np.sqrt(
        np.power( pdist[:,:,0], 2)
        + np.power( pdist[:,:,:1], 2))

    return pdist, N_pdist


if __name__ == '__main__':
    # get all check point values and sort them
    checkpoints= glob.glob("Box*.bin")
    check_point_values = np.sort(
    [ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])

    pn_file = "patch_network.dat"
    connections = gt.read_bonds("patch_network.dat")

    for j,val in enumerate(check_point_values[1:]):
        pos_i = np.fromfile("positions_{}.bin".format(val))
        pos_i = np.reshape(pos_i, (-1,3))
        pos_i = pos_i[:,:2]

        N_particles=len(pos_i)

        box_all = np.fromfile("Box_{}.bin".format(val))
        box_l = box_all[3:5]
        print("box", box_l)

        # Make a graph for time j: 
        connections_j = connections[j]
        G=nx.Graph()
        G.add_edges_from(connections_j[:,:2])
        particles_max_domain = gt.get_particles_in_largest_cluster(G)

        virtual_pos = calculate_pbc_images(pos_i,box_l,particles_max_domain)
        pdist, Ndist = get_pdist(virtual_pos)
        print("Ndist", Ndist)
        max_dist = np.max(Ndist)
        box_diag = np.linalg.norm(3*box_l)
        print(box_diag)
        frac_max_dist = max_dist/box_diag

        print("spanning length",max_dist)
        print("fraction of spanning",frac_max_dist)
