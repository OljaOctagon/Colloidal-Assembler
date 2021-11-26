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


def read_config(val):
    pos_i = np.fromfile("positions_{}.bin".format(val))
    pos_i = np.reshape(pos_i, (-1,3))
    pos_i = pos_i[:,:2]

    orient_i = np.fromfile("orientations_{}.bin".format(val))
    orient_i = np.reshape(orient_i, (-1,5))[:,4]


    return pos_i, orient_i


def get_triangle_vertices(pos_i, orient_i):
    vertices = np.zeros((3,2))

    l2 = 0.5
    h2_small = np.sqrt(3)/6
    h2_large = 2*h2_small

    rotmat_i = rotation_matrix(orient_i[i])
    ax0 = np.array([[-l2, l2, 0], [-h2_small, -h2_small, h2_large]])

    #ax0 = np.array([
    #    [-l2,-h2_small],
    #    [l2,-h2_small,
    #    [0, h2_large]])

    ax_n = get_orient(ax0, rotmat_i)
    ax_n = ax_n.transpose()

    vertices[0] = pos_i[i] + ax_n[0]
    vertices[1] = pos_i[i] + ax_n[1]
    vertices[2] = pos_i[i] + ax_n[2]


    return vertices



def get_rhombi_vertices(pos_i,orient_i):

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

    polygon_color = '#C17DCB'

    if not os.path.isdir("./frames"):
        os.mkdir("./frames")

    for j,val in enumerate(check_point_values[-6:]):

        pos_i, orient_i = read_config(val)
        N_particles = len(pos_i)


        parser = argparse.ArgumentParser(description='Particle drawing methods')
        parser.add_argument('-polygon_type', type=str, choices=['rhombus', 'triangle'])

        args = parser.parse_args()

        vertex_functions = {
        'triangle': get_triangle_vertices,
        'rhombus': get_rhombi_vertices}

        get_vertices = vertex_functions[args.polygon_type]

        fig,ax = plt.subplots()
        ax.set_aspect('equal', 'box')

        for i in range(N_particles):

            vertices = get_vertices(pos_i, orient_i)

            polygon = patches.Polygon(vertices, 
                linewidth=0.1, 
                edgecolor='k',
                facecolor=polygon_color, alpha=0.7)
            
            ax.add_patch(polygon)

        ax.scatter(pos_i[:,0], pos_i[:,1],s=1)
        ax.set_title("Frame {}".format(j))
            
        plt.axis("equal")
        plt.axis('off')

        plt.axis("equal")
        plt.savefig("./frames/frame_{}.png".format(j), dpi=1000)
        
        plt.cla()
        plt.clf()
        plt.close('all')

