import numpy as np
import pandas as pd
import argparse
import gel_tools as gt
import networkx as nx
import glob
import os
from itertools import islice
from itertools import cycle 

def push(i,j,rot,cell_lists):
    print(i,j)
    network_list=[]
    for klist in cell_lists:
        
        nklist=klist
        if rot==-1:
            nklist=klist[::-1]

        for ki,k in enumerate(nklist):
            im_j = j + k*N_particles
            print("list", nklist)
            starts_at_k = islice(cycle(nklist),ki+1,None)
            l = next(starts_at_k) 
            print("k, l", k,l)

            im_i = i + l*N_particles
            network_list.append([im_i,im_j])

    return network_list 

if __name__ == '__main__':

    # get all check point values and sort them
    checkpoints= glob.glob("Box*.bin")
    check_point_values = np.sort(
    [ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])

    pn_file = "patch_network.dat"
    connections = gt.read_bonds("patch_network.dat")


    frac_largest = []
    virtual_frac_largest = []


    cell_dict = {'top-bottom': [[2,1,8],[3,0,7],[4,5,6]],
                    'left-right':[[2,3,4],[1,0,5],[8,7,6]],
                    'diag-top-bottom':[[2,0,6],[3,5,8],[4,1,7]],
                    'diag-bottom-top': [[8,0,4],[1,3,6],[2,7,5]]}

    for j,val in enumerate(check_point_values[-1:]):
        print("checkpoint time ", val)
        pos = np.fromfile("positions_{}.bin".format(val))
        pos = np.reshape(pos, (-1,3))
        pos = pos[:,:2]
        N_particles=len(pos)
        box_all = np.fromfile("Box_{}.bin".format(val))
        box_l = box_all[3:5]

        # Make a graph for time j: 
        connections_j = connections[-1]
        connections_j = connections_j[:,:2]
        G=nx.Graph()
        G.add_edges_from(connections_j)
        particles_max_domain = gt.get_particles_in_largest_cluster(G)
        N_largest = len(particles_max_domain)

        frac_largest_i = N_largest/N_particles
        print(frac_largest_i)
        frac_largest.append(frac_largest_i)

        virtual_frac_largest_i = 0 
        virtual_patch_network = []
        N_images = 9

        for conn_j in connections_j:
            i,j = conn_j

            # Calculate distances between next neighbours
            # (bonded particles according to patch_network) 
            pdist = pos[j]-pos[i]

            x = pdist[0]
            y=  pdist[1]

            # define some shorthands 
            x_abs = np.fabs(x)
            x_sign = np.sign(x)
            y_abs = np.fabs(y)
            y_sign= np.sign(y)

            box_x=box_l[0]/2
            box_y=box_l[1]/2

            # if pbc not hit:
            # particles remain neighbours and are filled into new network edge array
            if (x_abs < box_x) and (y_abs < box_y):
                for k in range(N_images):
                    virtual_patch_network.append([i+k*N_particles, j+k*N_particles])

            # pbc are hit and particles are neighbours with image according to cell lists

            else: 
                # top-bottom
                if (x_abs < box_x) and (y_abs > box_y):
                    
                    print("top-bottom") 
                    if y_sign > 0:
                        rot = 1
                        print("rot=1")
                    if y_sign < 0:
                        print("rot=-1")
                        rot = -1

                    network_list_images = push(i,j,rot,cell_dict['top-bottom'])
                    virtual_patch_network.extend(network_list_images)

                # left-right
                if (x_abs > box_x) and (y_abs < box_y):
                    
                    print("left-right") 
                    # to the left
                    if x_sign > 0:
                        rot=-1
                    # to the right
                    if x_sign < 0:
                        rot=1

                    network_list_images = push(i,j,rot,cell_dict['left-right'])
                    virtual_patch_network.extend(network_list_images)

                # diagognal  
                if (x_abs > box_x) and (y_abs > box_y):
                    # left top to right bottom 
                    if x_sign != y_sign:
                        print("diagonal left top to right bottom") 
                        if x_sign > 0 and y_sign < 0:
                            rot = -1
                        if x_sign < 0 and y_sign > 0:
                            rot = 1
                        network_list_images = push(i,j,rot,cell_dict['diag-top-bottom'])
                        virtual_patch_network.extend(network_list_images)

                    # left bottom to right top 
                    if x_sign == y_sign:
                        print("left bottom to right top") 
                        if x_sign > 0:
                            rot = 1
                        if x_sign < 0:
                            rot = -1
                        network_list_images = push(i,j,rot,cell_dict['diag-bottom-top'])
                        virtual_patch_network.extend(network_list_images)

        print("Calculate bonded neighbour graph")
        # get virtual Graph
        G_virtual  = nx.Graph()
        G_virtual.add_edges_from(virtual_patch_network)

        print("Calculate largest cluster")
        virtual_max_domain = gt.get_particles_in_largest_cluster(G_virtual)
        virtual_frac_largest_i = len(virtual_max_domain)/(N_particles*N_images)

    virtual_frac_largest.append(virtual_frac_largest_i)

    print(frac_largest, virtual_frac_largest, check_point_values)
    with open("spanning.dat", 'w') as f:
        f.write("time,fraction_largest, fraction_largest_virtual\n")
        for j,val in enumerate(check_point_values[-1:]):
            f.write("{},{},{}\n".format(val,frac_largest[j], virtual_frac_largest[j]))
