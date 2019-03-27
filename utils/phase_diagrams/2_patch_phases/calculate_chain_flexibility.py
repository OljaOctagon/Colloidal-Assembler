
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import os.path


def files_exist(f_data):
    try:
        for filen in f_data:
            open(f_data[filen],'r')
        return True
    except IOError:
        #print("Error: file {} doesn't seem to exist.".format(f_data[filen]))
        return False

def parse_chains(fname):
    '''return chain input as list of lists''' 
    chains = []
    with open(fname, 'r') as f:
        for line in f:
            l_line = len(line.split(" "))
            if l_line == 3:
                chains.append([])

            if l_line == 1:
                chains[-1].append(int(line))

    return chains

def parse_config(f_data):
    """returns positions, orientations and box lengths
    as arrays per particle from binary data"""

    pos = np.fromfile(f_data['f_positions'])
    pos = np.reshape(pos, (-1,3))

    orient = np.fromfile(f_data['f_orientations'])
    orient = np.reshape(orient, (-1,5))[:,4]

    lbox = np.fromfile(f_data['f_box'])
    box_l = np.array([lbox[3], lbox[4]])

    return pos, orient, box_l


def calculate_distances_from_v(v,w, box_l):
    delta_pos_x = w[:,0] - v[0]
    delta_pos_y = w[:,1] - v[1]
    delta_pos_x = delta_pos_x - box_l[0]*np.rint(delta_pos_x/box_l[0])
    delta_pos_y = delta_pos_y - box_l[1]*np.rint(delta_pos_y/box_l[1])
    return  np.array([delta_pos_x, delta_pos_y])

def order_recursion(chain,pos, elem_0, rest_chain, ordered_chain, box_l):
    limit=1.11
    dist_0 = calculate_distances_from_v(pos[elem_0], pos[rest_chain], box_l)
    dist_0_max = np.maximum(np.fabs(dist_0[0,:]), np.fabs(dist_0[1,:]))
    arg_closest = np.where(dist_0_max < limit)[0].astype(int)
    closest = np.array(rest_chain)[arg_closest]
    # first element is end of chain

    if len(closest) == 1:
        ordered_chain.append(elem_0)
        rest_chain.remove(closest[0])
        ordered_chain, rest_chain = order_recursion(chain, pos, closest[0],rest_chain, ordered_chain, box_l)
        return ordered_chain, rest_chain

    # element somewhere in the middle 
    elif len(closest) == 2:
        ordered_chain.append(elem_0)
        rest_chain.remove(closest[0])
        rest_chain.remove(closest[1])

        ordered_chain, rest_chain = order_recursion(chain, pos, closest[0], rest_chain, ordered_chain, box_l)
        ordered_chain.reverse()
        ordered_chain, rest_chain = order_recursion(chain, pos, closest[1], rest_chain, ordered_chain, box_l)

        return ordered_chain, rest_chain

    # last element
    elif len(closest) == 0:
        ordered_chain.append(elem_0)
        return ordered_chain, rest_chain

    else:
        ordered_chain = []
        rest_chain = []
        print("Info: too many next neighbours: ", len(closest))
        return ordered_chain, rest_chain 

def order_chains(chains, pos, box_l):
    ordered_chains = []
    for chain in chains:
        ordered_chain = []
        rest_chain = chain[1:]
        ordered_chain, _ = order_recursion(chain,pos,chain[0], rest_chain, ordered_chain, box_l )

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

    f_data = {
        'f_cluster' : directory+"/All_Clusters_info.dat",
        'f_positions' : directory+"/positions.bin",
        'f_orientations' : directory+"/orientations.bin",
        'f_box' : directory+"/Box.bin" }

    chains = []
    df_chain=None
    if files_exist(f_data):
        chains =  parse_chains(f_data['f_cluster'])
        pos,orient,box_l = parse_config(f_data)
        # Order chains
        ordered_chains = order_chains(chains,pos, box_l)
        # List of flexibilities for every chain
        flex, rel_flex, kink, chain_length = calculate_flex_and_kink(ordered_chains, orient)
        # List of kinkyness for every chain
        #kink = calculate_kinkyness(ordered_chains, arr_orient)
        if flex is not None:
            df_chain = pd.DataFrame({'flex': flex, 'rel_flex': rel_flex, 'kink': kink, 'length': chain_length})
    return df_chain


def calculate_chain_parameters_MONO(patch_kind):
    phi_binary=0
    Chemical_Potentials = [0.3,0.4]
    Energies = [ 5.2,6.2,7.2,8.2]
    Patch_Positions = [0.2, 0.3 ,0.4 ,0.5, 0.6, 0.7, 0.8] 
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
                    directory = directory+'/analysis'
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
    Energies = [ 5.2,6.2,7.2,8.2]
    Patch_Positions = [0.3 ,0.5, 0.7] 
    Nruns=8

    is_first = True

    for phi_binary in Phi_Binary:
        for patch_pos_A in Patch_Positions:
            for patch_pos_B in Patch_Positions:
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
                            directory = directory+'/analysis'
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

    outfile="/Users/ada/Documents/Code_Development_2016/2D_patchy/chain_analysis/chains.csv"
    with open(outfile, 'a') as f:
        df_all_chains.to_csv(f, header=None)

