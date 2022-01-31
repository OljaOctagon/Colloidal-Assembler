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


def rotation_matrix(theta):
    rot_mat = np.zeros((2,2))

    rot_mat[0,0] = np.cos(theta) 
    rot_mat[0,1] = -np.sin(theta)
    rot_mat[1,0] = np.sin(theta)
    rot_mat[1,1] = np.cos(theta)

    return rot_mat


def get_orient(v, rot_mat):
    return rot_mat.dot(v)

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


def read_config(fdir, val):
    pos_i = np.fromfile("{}/positions_{}.bin".format(fdir, val))
    pos_i = np.reshape(pos_i, (-1,3))
    pos_i = pos_i[:,:2]

    orient_i = np.fromfile("{}/orientations_{}.bin".format(fdir, val))
    orient_i = np.reshape(orient_i, (-1,5))[:,4]


    return pos_i, orient_i


def get_rhombi_vertices(pos_i,orient_i,i):

    def get_edge_points(pos_i,ax_n,sign_p):
        vertex_n = np.zeros(2)
        vertex_n = pos_i + sign_p[0]*ax_n[:,0]/2. + sign_p[1]*ax_n[:,1]/2.
        return vertex_n

    sin60 = np.sin(np.pi/3.)
    cos60 = np.cos(np.pi/3.)

    ax0 = np.array([[1,cos60],[0,sin60]])
    vertices = np.zeros((4,2))
    ax_n = np.zeros((2,2))

    rotmat_i = rotation_matrix(orient_i[i])
    ax_n = get_orient(ax0, rotmat_i)

    vertices[0] = get_edge_points(pos_i[i],ax_n,np.array([-1,-1]))
    vertices[1] = get_edge_points(pos_i[i],ax_n,np.array([+1,-1]))
    vertices[2] = get_edge_points(pos_i[i],ax_n,np.array([+1,+1]))
    vertices[3] = get_edge_points(pos_i[i],ax_n,np.array([-1,+1]))

    return vertices




def get_intersect(pos_i,all_vertices,id_i,sx,sy,r_sphere,a_rhombi,rhombi_long_diagonal):

    intersect=False 
    dist  = pos_i[id_i] - sp 
    dist = dist - np.array([box_lx,box_ly])*np.fint(
        dist/np.array([box_lx,box_ly]))

    # test outer sphere of rhombi as coarse overlap test
    if np.abs(dist)>rhombi_long_diagonal:
        intersect=False 

    else:

        nearest_point = sp + (dist/np.linalg.norm(dist))*r_sphere 

        e01 = (all_vertices[id_i][1] - all_vertices[id_i][0])/a_rhombi
        e02 = (all_vertices[id_i][3] - all_vertices[id_i][0])/a_rhombi

        nearest_point_rhs = np.array([np.dot(nearest_point,e01),
            np.dot(nearest_point,e01)])

        length_nearest_point_rhs = np.linalg.norm(nearest_point_rhs)

        if length_nearest_point_rhs < 1:
            intersect=True
    
    return intersect 


if __name__ == '__main__':
    
    # get all check point values and sort them
    filedir = "/Users/ada/Documents/Code_Development_2020/rhombi/percolations_study/vsc3/percolation_runs/copy_dir"
    #filedir_TU = "/home/carina/Documents/2D_patchy/percolation_study/vsc3/runs/runs"

    phi=0.125
    delta=0.2
    T=0.01 
    #ptypes= ['double_manta_asymm_1', 'double_mouse_asymm_1', 'double_mouse_symm_1', 'double_mouse_symm_2']
    ptype = 'double_mouse_asymm_1'
    filei= dirn="{}/{}_phi_{}_delta_{}_temp_{}".format(filedir,ptype,phi,delta,T)

    val=17600000
    box_all = np.fromfile("{}/Box_{}.bin".format(filei, val))
    box_x_center = box_all[0]  
    box_y_center = box_all[1]    

    box_lx = box_all[3]
    box_ly = box_all[4]



    pos_i, orient_i = read_config(filei, val)
    N_particles = len(pos_i)

    N_vertices = []

    # Make cell lists 
    r_sphere=0.5
    a_rhombi = 1.0
    alpha = np.pi/3
    rhombi_long_diagonal = a_rhombi * np.sqrt(np.power(1+np.cos(alpha)+np.sin(alpha),2))

    L_cell_min = r_sphere + rhombi_long_diagonal/2 

    Nx_cell = int(np.floor(box_lx/L_cell_min))
    Ny_cell = int(np.floor(box_ly/L_cell_min))

    L_cell_x = box_lx/Nx_cell
    L_cell_y = box_ly/Ny_cell 


    cell_ID_of_particle = np.array([int.np.floor(pos_i[:,0]/L_cell_x)
        , int(np.floor(pos_i[:,1]/L_cell_y))

    from collections import defaultdict 

    cell_list=defaultdict(list)

    # fill cell list 
    for i,entry in enumerate(cell_ID_of_particle):
        cell_list[(entry[0],entry[1])].append(i)

    all_vertices = []
    for i in range(N_particles):
        all_vertices.append(get_rhombi_vertices(pos_i,orient_i,i))


    print("box x-length: {}, box y_lenght: {}, box_center: {},{}".format(box_lx,box_ly, box_x_center, box_y_center))

    print("Number of cells in x-direction: {}, side-length of cells in x-direction: {}".format(Nx_cell,L_Cell_x))

    print("Number of cells in y-direction: {}, side-length of cells in y-direction: {}".format(Ny_cell,L_Cell_y))


    Nsteps = 100
    pore_cloud_centers = []


    print("Box start and end coordinates:")
    box_xstart = box_x_center - box_lx/2.
    box_xend = box_x_center + box_lx/2

    box_ystart = box_y_center - box_ly/2.
    box_yend = box_y_center + box_ly/2. 

    print("Box xstart, Box xend: {},{}".format(box_xstart,box_xend))
    print("Box ystart, Box yend: {},{}".format(box_ystart,box_yend))


    print("number of generated spheres: {}".format(Nsteps))
    for istep in Nsteps:

        print(istep)
        # pick center of sphere 
        sx = np.random.uniform(box_xstart,box_xend)
        sy = np.random.uniform(box_ystart,box_yend)
        sp=np.array([sx,sy])

        cell_of_sphere = (int.np.floor(sx/L_cell_x),
            int(np.floor(sy/L_cell_y)))

        overlap_candidates = []
        # add particles of original cell 
        overlap_candidates.append(cell_list[cell_of_sphere])

        # get particles from neighbour cells 
        neigbour_cells = [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,-1),(-1,1),(1,1)]

        for ncell in neighbour_cells: 
            ncell = cell_list + neighbour_cells 
            ncell = (left_neighbour[0]%Nx_cell, left_neighbour[1]%Ny_cell)
            overlap_candidates.append(cell_list[ncell])

        global_intersect=False 
        ic=0
        while (global_intersect == False and ic < len(overlap_candidates)):
            
            id_i = overlap_candidates[ic]       

            global_intersect=get_intersect(pos_i,all_vertices,id_i,sp, r_sphere,a_rhombi,rhombi_long_diagonal)

            ic=ic+1 

        if global_intersect==False:
            pore_cloud_centers.append(sp)


    print("Fraction of non_overlapping spheres: {}".format(len(pore_cloud_centers)/Nsteps))

    print("plotting")
    fig,ax = plt.subplots(figsize=(20,20))
    ax.set_aspect('equal', 'box')
    polygon_color = 'cyan'
    #polygon_color = "#04CC80"

    for i in range(N_particles):
        vertices = get_rhombi_vertices(pos_i, orient_i, i)
        polygon = patches.Polygon(vertices,
            linewidth=0.1,
            edgecolor='k',
            facecolor=polygon_color, alpha=0.7)

        ax.add_patch(polygon)   
    

    for sphere_i in pore_cloud_centers:
        patches.Circle((sphere_i[0],sphere_i[1]), radius=r_sphere,color='red',edgecolor='k', alpha=0.5)
    

    plt.axis("equal")
    #plt.axis('off')

    plt.savefig("{}/rhombi_pore_test.png".format(filedir), dpi=100)
    plt.show()







                    








