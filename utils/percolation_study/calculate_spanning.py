import numpy as np
import pandas as pd
import argparse
import gel_tools as gt
import networkx as nx
import glob
import os

def get_edge_points(pos_i,ax_n,sign_p):
    edge_n = np.zeros(2)
    edge_n = pos_i + sign_p[0]*ax_n[:,0]/2. + sign_p[1]*ax_n[:,1]/2.

    return edge_n

def get_orient(v, rot_mat):
    return rot_mat.dot(v)


def rotation_matrix(theta):
    rot_mat = np.zeros((2,2))

    rot_mat[0,0] = np.cos(theta) 
    rot_mat[0,1] = -np.sin(theta)
    rot_mat[1,0] = np.sin(theta)
    rot_mat[1,1] = np.cos(theta)

    return rot_mat


def get_patches(pid,orient_pid,pos_pid):

    patches_pid = np.zeros((4,2))
    sin60 = np.sin(np.pi/3.)
    cos60 = np.cos(np.pi/3.)

    ax0 = np.array([[1,cos60],[0,sin60]])
    edges = np.zeros((4,2))
    ax_n = np.zeros((2,2))

    rotmat_pid = rotation_matrix(orient_pid)
    ax_n = get_orient(ax0, rotmat_pid)

    edges[0] = get_edge_points(pos_pid,ax_n,np.array([-1,-1]))
    edges[1] = get_edge_points(pos_pid,ax_n,np.array([+1,-1]))
    edges[2] = get_edge_points(pos_pid,ax_n,np.array([+1,+1]))
    edges[3] = get_edge_points(pos_pid,ax_n,np.array([-1,+1]))

    # dma as1 type

    # patch type 1 
    patches_pid[0,:] = edges[0] + 0.2*(edges[3]-edges[0])
    patches_pid[1,:] = edges[2] + 0.8*(edges[3]-edges[2])

    # patch type 2 
    patches_pid[2,:] = edges[0] + 0.2*(edges[1]-edges[0])
    patches_pid[3,:] = edges[1] + 0.2*(edges[2]-edges[1])

    return patches_pid



def calculate_pbc_images(pos_i,orient_i,box_l):

    # Periodic images of all particles 
    N_images = 9
    N_patches = 4 
    N_virtual = N_particles*N_images
    virtual_pos = np.zeros((N_virtual,2))
    virtual_orient = np.zeros((N_virtual))

    sign_array = np.array([[0,0],
                            [-1,0],
                            [-1,-1],
                            [0,-1],
                            [1,-1],
                            [1,0],
                            [1,1],
                            [0,1],
                            [-1,1]])

    for pid in (range(N_particles)):
        for image_j in range(N_images):
            new_id = pid+image_j*N_particles
            virtual_pos[new_id] = pos_i[pid] + sign_array[image_j]*box_l
            virtual_orient[new_id] = orient_i[pid]

    return virtual_pos, virtual_orient


if __name__ == '__main__':

    # get all check point values and sort them
    checkpoints= glob.glob("Box*.bin")
    check_point_values = np.sort(
    [ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])

    pn_file = "patch_network.dat"
    connections = gt.read_bonds("patch_network.dat")


    frac_largest = []
    virtual_frac_largest = []

    for j,val in enumerate(check_point_values[-1:]):
        print("checkpoint time ", val)
        pos_i = np.fromfile("positions_{}.bin".format(val))
        pos_i = np.reshape(pos_i, (-1,3))
        pos_i = pos_i[:,:2]
        orient_i = np.fromfile("orientations_{}.bin".format(val))
        orient_i = np.reshape(orient_i, (-1,5))[:,4]
        N_particles=len(pos_i)
        box_all = np.fromfile("Box_{}.bin".format(val))
        box_l = box_all[3:5]

        # Make a graph for time j: 
        connections_j = connections[-1]
        G=nx.Graph()
        G.add_edges_from(connections_j[:,:2])
        particles_max_domain = gt.get_particles_in_largest_cluster(G)
        N_largest = len(particles_max_domain)

        frac_largest_i = N_largest/N_particles
        print(frac_largest_i)
        frac_largest.append(frac_largest_i)
        virtual_frac_largest_i = 0 
        if frac_largest_i > 0.5:
            sin60 = np.sin(np.pi/3.)
            cos60 = np.cos(np.pi/3.)
            N_images = 9 

            print("Calculate virtual positions")
            cutoff = np.sqrt(np.power((1 + cos60),2) + np.power(sin60,2)) + 0.1
            patch_cutoff = 0.1
            virtual_pos, virtual_orient = calculate_pbc_images(pos_i,orient_i,box_l)
            virtual_box = box_l*3
            N_patches = 4
            N_virtual = N_particles*N_images

            virtual_patch_network = []

            print("Calculate next neighbours and bonds")
            # Next get all particles that are next neighbours, and calc if patches meet
            for pid in range(N_virtual):
                patches_pid = get_patches(pid,virtual_orient[pid],virtual_pos[pid])
                for nid in range(N_virtual):
                    if pid<nid:
                        pdist = virtual_pos[pid] - virtual_pos[nid]
                        pdist = pdist - virtual_box*np.rint(pdist/virtual_box)
                        Ndist = np.linalg.norm(pdist)
                        if Ndist<cutoff: 
                            patches_nid = get_patches(nid,virtual_orient[nid],virtual_pos[nid])
                            for patch_pid in range(N_patches):
                                for patch_nid in range(N_patches):
                                    p = patch_pid - patch_nid
                                    patch_dist = p - virtual_box*np.rint(p/virtual_box)
                                    N_patch_dist = np.linalg.norm(patch_dist)
                                    if N_patch_dist < patch_cutoff:
                                        virtual_patch_network.append([pid,nid])


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
