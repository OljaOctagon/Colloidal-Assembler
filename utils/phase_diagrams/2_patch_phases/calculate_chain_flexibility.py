
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import os.path
from collections import defaultdict 
import glob
def files_exist(filen):
    try:
        open(filen,'r')
        return True
    except IOError:
        print("Error: file {} doesn't seem to exist.".format(filen))
        return False

def parse_chains(fname):
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
    min_check_point = min(checkpoints)

    freq=1000

    pos = defaultdict(np.array)
    box_l = defaultdict(np.array)
    orient = defaultdict(np.array)

    for point in range(min_check_point+freq,max_check_point+freq,freq):
        try:
            pos_i = np.fromfile("{}/{}_{}.bin".format(directory, 'positions', point))
            pos_i = np.reshape(pos_i, (-1,3))
            pos[point] = pos_i

            orient_i = np.fromfile("{}/{}_{}.bin".format(directory, 'orientations', point))
            orient_i = np.reshape(orient_i, (-1,5))[:,4]
            orient[point] = orient_i

            lbox_i = np.fromfile("{}/{}_{}.bin".format(directory, 'Box', point))
            box_l_i = np.array([lbox_i[3], lbox_i[4]])
            box_l[point] = box_l_i
            max_check_point = point
        except:
            max_check_point = point

    return pos,orient, box_l


def rotate_to_particle_frame(w, orient_vi):
    orient_vi = orient_vi -  np.pi/2.

    rot_mat = np.zeros((2,2))
   
    rot_mat[0,0] = np.cos(orient_vi)
    rot_mat[0,1] = -np.sin(orient_vi)

    rot_mat[1,0] = np.sin(orient_vi)
    rot_mat[1,1] = np.cos(orient_vi)

    wprime = np.dot(rot_mat,w)
    return wprime

def calculate_distances_from_v(v,w, box_l, orient_vi):

    delta_pos_x = w[:,0] - v[0]
    delta_pos_y = w[:,1] - v[1]
    delta_pos_x = delta_pos_x - box_l[0]*np.rint(delta_pos_x/box_l[0])
    delta_pos_y = delta_pos_y - box_l[1]*np.rint(delta_pos_y/box_l[1])

    w  = np.reshape(np.array([delta_pos_x, delta_pos_y]), (2,-1))
    
    arr = rotate_to_particle_frame(w,orient_vi)
    #np.sqrt(np.power(delta_pos_x,2) + np.power(delta_pos_y,2))
    return arr

def order_recursion(chain,pos_i, orient_i, elem_0, rest_chain, ordered_chain, box_l_i):
    limit=1.1
    limit_low=0.8
    if len(chain): 
        dist_0 = calculate_distances_from_v(pos_i[elem_0], pos_i[rest_chain], box_l_i, orient_i[elem_0])
        dist_0_max = np.maximum(np.fabs(dist_0[0,:]), np.fabs(dist_0[1,:]))
        arg_closest = np.where((dist_0_max < limit) & ( dist_0_max > limit_low))[0].astype(int)
        closest = np.array(rest_chain)[arg_closest]
        # first element is end of chain

    else:
        closest = 0

    if len(closest) == 1:
        ordered_chain.append(elem_0)
        rest_chain.remove(closest[0])
        ordered_chain, rest_chain = order_recursion(chain, pos_i, orient_i, closest[0],rest_chain, ordered_chain, box_l_i)
        return ordered_chain, rest_chain

    # element somewhere in the middle 
    elif len(closest) == 2:
        ordered_chain.append(elem_0)
        rest_chain.remove(closest[0])
        rest_chain.remove(closest[1])

        ordered_chain, rest_chain = order_recursion(chain, pos_i, orient_i, closest[0], rest_chain, ordered_chain, box_l_i)
        ordered_chain.reverse()
        ordered_chain, rest_chain = order_recursion(chain, pos_i, orient_i, closest[1], rest_chain, ordered_chain, box_l_i)

        return ordered_chain, rest_chain

    # last element
    elif len(closest) == 0:
        ordered_chain.append(elem_0)
        return ordered_chain, rest_chain

    else:
        ordered_chain = []
        rest_chain = []
        print("Info: too many next neighbours: ", len(closest), dist_0[:,arg_closest])
        return ordered_chain, rest_chain 

def order_chains(chains_i, pos_i, orient_i, box_l_i):
    ordered_chains = []
    for chain in chains_i:
        ordered_chain = []
        rest_chain = chain[1:]
        ordered_chain, _ = order_recursion(chain,pos_i, orient_i, chain[0], rest_chain, ordered_chain, box_l_i )
        if len(ordered_chain) > 0:
            ordered_chains.append(ordered_chain)

    return ordered_chains 

def calculate_flex_and_kink(chains, orient):

    def flex_and_kink(chain):

        rel_delta_angle = orient[chain[:-1]] - orient[chain[1:]]
        delta_angle = np.fabs(orient[chain[:-1]] - orient[chain[1:]])

        parallel_angles=np.array([0, np.pi, 2*np.pi,0])
        non_parallel_angles = np.array([np.pi/3, 2*np.pi/3, 4*np.pi/3, 5*np.pi/3])
        all_angles = np.concatenate((parallel_angles, non_parallel_angles)) 
        chain_bonds = np.zeros(len(delta_angle))
        rel_angle = np.zeros(len(delta_angle))
        for i, angle in enumerate(delta_angle):
            #abs
            parallel_min = np.min(np.fabs(parallel_angles - angle))
            non_parallel_min = np.min(np.fabs(non_parallel_angles - angle))
            pnp_min = np.array([parallel_min, non_parallel_min])
            delta_angle[i] = np.min(pnp_min)
            chain_bonds[i] = np.argmin(pnp_min)


            #rel
            rel_v = angle - all_angles
            rel_angle[i]=rel_v[np.argmin(np.fabs(rel_v))]

        flex = np.average(delta_angle)
        kink = np.average(chain_bonds)
        rel_flex = np.fabs(np.average(rel_angle))

        return flex, rel_flex, kink


    flex_kink = np.array([ np.array(flex_and_kink(chain)) for chain in chains if len(chain)>1]) 
    chain_length = [ len(chain) for chain in chains if len(chain)>1]
    if flex_kink.shape[0] > 0:
        return flex_kink[:,0], flex_kink[:,1], flex_kink[:,2], chain_length

    else:
        return None, None, None, None

def aggregate_chain_parameters(directory):

    df_chain=None
    chains =  parse_chains("{}/All_Clusters_info.dat".format(directory))
    pos,orient,box_l  = parse_config(directory)

    flex = np.array([])
    kink = np.array([])
    rel_flex = np.array([])
    chain_length = np.array([])

    for i in chains.keys():
        # Order chains
        ordered_chains = order_chains(chains[i],pos[i], orient[i], box_l[i])
        uoc = [l for l in chains[i] if len(l)>1]
        # List of flexibilities for every chain
        flex_i, rel_flex_i, kink_i, chain_length_i = calculate_flex_and_kink(ordered_chains, orient[i])
        # List of kinkyness for every chain
        #kink = calculate_kinkyness(ordered_chains, arr_orient)
        oc = [ l for l in ordered_chains if len(l)>1 ]
        flex = np.concatenate((flex,flex_i))
        kink = np.concatenate((kink, kink_i))
        rel_flex = np.concatenate((rel_flex, rel_flex_i))
        chain_length = np.concatenate((chain_length, chain_length_i))

    if flex is not None:
        df_chain = pd.DataFrame({'flex': flex, 'rel_flex': rel_flex, 'kink': kink, 'length': chain_length})
    return df_chain


def calculate_chain_parameters_MONO(patch_kind):
    phi_binary=0
    Chemical_Potentials = [0.3]
    Energies = [8.2]
    Patch_Positions = [0.2,0.3,0.4,0.5,0.6,0.7,0.8] 
    Nruns=8

    is_first = True

    for patch_pos in Patch_Positions:
        for mu in Chemical_Potentials:
            for energy in Energies:
                for run_i in range(1,Nruns+1):
                    print(patch_kind, patch_pos, mu, energy, run_i)
                    directory = "mu_{}Energy_{}{}_patchpos_{}_{}".format(mu,
                                                                            energy,
                                                                            patch_kind,
                                                                            patch_pos,
                                                                            run_i)

                    if files_exist("{}/All_Clusters_info.dat".format(directory)):
                        df_chain = aggregate_chain_parameters(directory)

                        if df_chain is not None:
                            df_chain['mu'] = mu
                            df_chain['patch_pos_A'] = patch_pos
                            df_chain['patch_pos_B'] = patch_pos
                            df_chain['energy'] = energy
                            df_chain['run'] = run_i
                            df_chain['topology'] = patch_kind
                            df_chain['phi_binary'] = phi_binary

                            if is_first == False:
                                df_all_chains = pd.concat([df_all_chains, df_chain])

                            if is_first == True:
                                df_all_chains = df_chain
                                is_first = False


    return df_all_chains

def calculate_chain_parameters_BINARY(patch_kind):
    Phi_Binary=[0.25, 0.5,0.75]
    Chemical_Potentials = [0.3]
    Energies = [8.2]
    Patch_Positions = [0.3 ,0.5, 0.7] 
    Nruns=8

    is_first = True

    for phi_binary in Phi_Binary:
        for patch_pos_A in Patch_Positions:
            patch_pos_B = patch_pos_A
            #for patch_pos_B in Patch_Positions:
            for mu in Chemical_Potentials:
                for energy in Energies:
                    for run_i in range(1,Nruns+1):
                        print(patch_kind, phi_binary, patch_pos_A, patch_pos_B, mu, energy, run_i)

                        directory = "mu_{}Energy_{}{}_patch_posA_{}posB_{}bin_{}_{}".format(mu,
                                                                                        energy,
                                                                                        patch_kind,
                                                                                        patch_pos_A,
                                                                                        patch_pos_B,
                                                                                        phi_binary,
                                                                                        run_i)

                        if files_exist("{}/All_Clusters_info.dat".format(directory)):
                            df_chain = aggregate_chain_parameters(directory)
                            if df_chain is not None:
                                df_chain['mu'] = mu
                                df_chain['patch_pos_A'] = patch_pos_A
                                df_chain['patch_pos_B'] = patch_pos_B
                                df_chain['energy'] = energy
                                df_chain['run'] = run_i
                                df_chain['topology'] = patch_kind
                                df_chain['phi_binary'] = phi_binary

                                if is_first == False:
                                    df_all_chains = pd.concat([df_all_chains, df_chain])

                                if is_first == True:
                                    df_all_chains = df_chain
                                    is_first = False


    return df_all_chains

import argparse 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-binary', type=str, default='off')
    parser.add_argument('-patch', type=str, default='symm')
    args = parser.parse_args()

    if args.binary == 'off':
        df_all_chains = calculate_chain_parameters_MONO(args.patch)

    if args.binary == 'on':
        df_all_chains = calculate_chain_parameters_BINARY(args.patch)

    outfile="/Users/ada/Documents/Code_Development_2016/2D_patchy/chain_analysis/nchains.csv"
    with open(outfile, 'a') as f:
        df_all_chains.to_csv(f, header=None)

