import numpy as np
import pandas as pd
import argparse
import gel_tools as gt
import networkx as nx
import glob
import os

# calculate if spanning:
# use the virtual position
# get pdist and get patch dist
# calculate all connections again
# get largest cluster


# NOTE: probalbly deprecated
def calculate_pbc_images(pos_i,box_l,patches_i, particles_max_domain):

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

# NOTE: probalbly deprecated
def get_pdist(vpos):
    l = len(vpos)
    pdist = np.zeros((l,l,2))
    for i in range(2):
        p = np.reshape(vpos[:,i], (l,1))
        pdist = p - p.transpose()

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


    frac_largest = []
    virtual_frac_largest = []

    for j,val in enumerate(check_point_values[1:]):
        pos_i = np.fromfile("positions_{}.bin".format(val))
        pos_i = np.reshape(pos_i, (-1,3))
        pos_i = pos_i[:,:2]
        orient_i = np.fromfile("orientations_{}.bin".format(val))
        orient_i = np.reshape(orient_i, (-1,5))[:,4]
        N_particles=len(pos_i)
        box_all = np.fromfile("Box_{}.bin".format(val))
        box_l = box_all[3:5]

        # Make a graph for time j: 
        connections_j = connections[j]
        G=nx.Graph()
        G.add_edges_from(connections_j[:,:2])
        particles_max_domain = gt.get_particles_in_largest_cluster(G)
        N_largest = len(particles_max_domain)

        frac_largest_i = N_largest/N_particles
        frac_largest.append(frac_largest_i)
        virtual_frac_largest_i = 0 
        if frac_largest_i > 0.5:

            patches_i = np.zeros((N,4,2))
            sin60 = np.sin(np.pi/3.)
            cos60 = np.cos(np.pi/3.)
            cutoff = np.sqrt(np.power((1 + cos60),2) + np.power(sin60,2)) + 0.1
            patch_cutoff = 0.1

            ax0 = np.array([[1,cos60],[0,sin60]])
            edges = np.zeros((4,2))
            ax_n = np.zeros((2,2))

            for pid in particles_max_domain:
                rotmat_i = rotation_matrix(orient_i[pid])
                ax_n = get_orient(ax0, rotmat_i)

                edges[0] = get_edge_points(pos_i[pid],ax_n,np.array([-1,-1]))
                edges[1] = get_edge_points(pos_i[pid],ax_n,np.array([+1,-1]))
                edges[2] = get_edge_points(pos_i[pid],ax_n,np.array([+1,+1]))
                edges[3] = get_edge_points(pos_i[pid],ax_n,np.array([-1,+1]))

                # dma as1 type

                # patch type 1 
                patches_i[pid,0,:] = edges[0] + 0.2*(edges[3]-edges[0])
                patches_i[pid,1,:] = edges[2] + 0.8*(edges[3]-edges[2])

                # patch type 2 
                patches_i[pid,2,:] = edges[0] + 0.2*(edges[1]-edges[0])
                patches_i[pid,3,:] = edges[1] + 0.2*(edges[2]-edges[1])


           for pid1 in particles_max_domain:
               for pid2 in particles_max_domain:


            









''' 
            # calculate virtual connections of image particles 
            for connection in connections_j:
                id1_0 = connection[0]
                id2_0 = connection[1]
                if id1_0 in particles_max_domain and id2_0 in particles_max_domain:
                    for image_j in range(N_images):
                        id1 = id1_0 + image_j*N_largest
                        id2 = id2_0 + image_j*N_largest
                        virtual_connections.append([id1,id2])

            # make a new graph of virtual particles, calculate the largest
            G_virtual  = nx.Graph()
            G_virtual.add_edges_from(virtual_connections)
            virtual_max_domain = gt.get_particles_in_largest_cluster(G_virtual)
            virtual_frac_largest_i = len(virtual_max_domain)/(N_particles*N_images)

        virtual_frac_largest.append(virtual_frac_largest_i)

    with open("spanning.dat", 'w') as f:
        f.write("time,fraction_largest, fraction_largest_virtual\n")
        for j,val in enumerate(check_point_values[1:]):
            f.write("{},{},{}\n".format(val,frac_largest[j], virtual_frac_largest[j]))
'''
