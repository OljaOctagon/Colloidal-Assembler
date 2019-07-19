import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import os.path
import glob
from collections import defaultdict

def parse_clusters(fname):
    '''return clusters as dictionary of (list of lists) for every checkpoint''' 
    clusters = defaultdict(list)

    with open(fname, 'r') as f:
        for line in f:
            line_list = line.split(" ")
            l_line = len(line_list)
            if l_line == 3:
                cp_i = int(line_list[0])
                clusters[cp_i].append([])

            if l_line == 1:
                clusters[cp_i][-1].append(int(line))

    return clusters

def parse_config(directory):
    """returns positions, box lengths
    as dictionary of arrays. Every dictionary entry is a checkpoint."""

    filelist = glob.glob(directory+"/positions_*.bin")
    g = lambda x: x.split("_")[-1].split(".")[0]

    checkpoints = [ int(g(l)) for l in filelist]
    max_check_point = max(checkpoints)

    start=max(8000000, min(checkpoints))
    freq=20000
    end= max_check_point

    pos = defaultdict(np.array)
    box_l = defaultdict(np.array)

    for point in range(start,end,freq):
        pos_i = np.fromfile("{}/{}_{}.bin".format(directory, 'positions', point))
        pos_i = np.reshape(pos_i, (-1,3))
        pos[point] = pos_i

        lbox_i = np.fromfile("{}/{}_{}.bin".format(directory, 'Box', point))
        box_l_i = np.array([lbox_i[3], lbox_i[4]])
        box_l[point] = box_l_i

    return pos,box_l

def get_pdist(vpos,box):
    l = len(vpos)
    pdist = np.zeros((l,l,3))

    for i in range(3):
        p = np.reshape(vpos[:,i], (l,1))
        p = p - p.transpose()
        pdist[:,:,i] = p - box[i]*np.rint(p/box[i])
        
    N_pdist = np.sqrt(
        np.power( pdist[:,:,0], 2)
        + np.power( pdist[:,:,1], 2)
        + np.power( pdist[:,:,2], 2) )

    return pdist, N_pdist

def calc_distance(cluster, pos_i, boxl_i):
    l = len(cluster)
    distance = np.zeros((l,l))
    pdist = np.zeros((l,l,2))
    for i in range(2):
       p = np.reshape(pos_i[cluster,i], (l,1))
       p = p - p.transpose()
       pdist[:,:,i] = p - boxl_i[i]*np.rint(p/boxl_i[i])

    distance = np.sqrt(
        np.power( pdist[:,:,0], 2)
        + np.power( pdist[:,:,1], 2))

    return distance

def is_loop(cluster, pos_i, boxl_i):
    is_loop = True
    distance = calc_distance(cluster, pos_i, boxl_i)
    limit = 1.5
    nneigh = np.array([ len(np.where(distance[i]<limit)[0]) for i in range(len(distance))])
    nneigh_lt_3 = np.where(nneigh < 3)[0].tolist()
    #print("cluster")
    if nneigh_lt_3:
        #print('no loop')
        is_loop = False

    return is_loop

def calculate_center_of_mass(checkpoint_i, pos_i, boxl_i, clusters_i, directory):

    # take only length five and six clusters
    six_p_clusters = [ cluster for cluster in clusters_i if len(cluster) == 6 or len(cluster) == 5 ]
    six_p_clusters  = [ cluster for cluster in clusters_i if is_loop(cluster, pos_i, boxl_i)]
    # calculate center of mass for all length six clusters and
    # write to file
    with open(directory+"/center_of_mass_stars.dat", 'a') as f:
        for c_id, cluster in enumerate(six_p_clusters):

            theta_x = (pos_i[cluster,0]/boxl_i[0])*2*np.pi
            theta_y = (pos_i[cluster,1]/boxl_i[1])*2*np.pi

            xp1 = np.mean( (boxl_i[0]/(2*np.pi))*np.cos(theta_x) )
            xp2 = np.mean( (boxl_i[0]/(2*np.pi))*np.sin(theta_x) )

            yp1 = np.mean( (boxl_i[1]/(2*np.pi))*np.cos(theta_y) )
            yp2 = np.mean( (boxl_i[1]/(2*np.pi))*np.sin(theta_y) )

            theta_av_x = np.arctan2(-xp2, -xp1) + np.pi
            theta_av_y = np.arctan2(-yp2, -yp1) + np.pi

            center_mass_x = (boxl_i[0]/(2*np.pi))*theta_av_x
            center_mass_y = (boxl_i[1]/(2*np.pi))*theta_av_y

            f.write("{} {} {} {} {} {} {} \n".format(checkpoint_i,
                                                    len(six_p_clusters),
                                                    boxl_i[0],
                                                    boxl_i[1],
                                                    c_id,
                                                    center_mass_x,
                                                    center_mass_y))

def run(directory):
    clusters =  parse_clusters(directory+"/All_Clusters_info.dat")
    pos,box_l = parse_config(directory)

    # calculate center of mass of all particles for every checkpoint 
    for key in pos:
        calculate_center_of_mass(key, pos[key], box_l[key], clusters[key], directory)
    #except:
    #    print('Warning: can not open file')

if __name__ == "__main__":
    dirlist = glob.glob("mu_0.4Energy_9.2symm_patchpos_0.*_Pressure_*")

    for dir_n in dirlist:
        print(dir_n)
        run(dir_n)
