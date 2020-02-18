import numpy as np
import pandas as pd
import glob
import re
import os
import argparse
import networkx as nx
import itertools

def check_file(fname):
    try:
        open(fname,"r")
        return 1
    except IOError:
        print("Error: File doesn't seem to exist.")
        return 0

cluster_type_dict = {
    '6-star': 6,
    '5-star': 5,
    'box': 3,
    'open-box':3,
    'open-6-star':6,
    'open-5-star': 5,
    'open-jenga':4,
    'l-shape': 2,
    'line-shape':2,
    'brick-shape':2}

patch_dict = {
    'mouse': [0,2],
    'manta': [0,1]
}

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


def filter(arr, particle_type):
    '''
    Filters out all non-parallel mouse bonds.

    Parameters
    ----------
    arr: array-like
        encodes the type of bonding,
        4 columns (id1,id2,pid1,pid2)

    Returns
    -------
    arr
       filtered array

    '''

    patch = patch_dict[particle_type]

    new_arr1 = arr[np.logical_and(arr[:,2]==patch[0], arr[:,3]==patch[1])]
    new_arr2 = arr[np.logical_and(arr[:,2]==patch[1], arr[:,3]==patch[0])]
    new_arr = np.concatenate((new_arr1,new_arr2))

    return new_arr

def calculate_yield(arr, cluster_type, N_particles):
    '''
    Calculates yield of specified cluster type.

    Parameters
    ----------
    arr: array like
        encodes the type of bonding,
        4 columns (id1,id2,pid1,pid2)

    cluster_type: str
        The cluster type of which the yield should be calculated.
        can be: "6-stars", "5-stars", boxes"

    particle_type: str
        name of the particle type.
        can be: 'manta', 'mouse

    N_particles: int
        Total number of particles in the system.

    Returns
    -------
    float:
         the yield of the specified cluster type

    Raises
    ------
    ValueError:
        when cluster type is not known.

    '''

    loop_shapes = ['5-star', '6-star', 'box', 'open-box', 'open-6-star', 'open-5-star', 'open-jenga']
    dimer_shapes = ['l-shape', 'brick-shape', 'line-shape']

    if cluster_type in loop_shapes:

        # Filter out bonds and construct graph
        m_arr = filter(arr, args.ptype)
        G_m = nx.Graph()
        G_m.add_edges_from(m_arr[:,:2])
        DG = nx.DiGraph(G_m)
        loops = list(nx.simple_cycles(DG))


        clusters = [ loop for loop in loops if len(loop)==cluster_type_dict[cluster_type] ]
        clusters = [ sorted(item) for item in clusters]
        clusters = list(clusters for clusters,_ in itertools.groupby(clusters))
        N_clusters=len(clusters)
        N_stars = N_clusters*cluster_type_dict[cluster_type]
        p_yield = N_stars/N_particles
        return p_yield

    if cluster_type in dimer_shapes:
        print("dimer!")
        labels = [ (i,j) for i,j in arr[:,2:]]

        patch_dict_2p = {
        'l-shape': (0,1),
        'line-shape': (1,1),
        'brick-shape': (0,0)
        }

        G_m = nx.Graph()
        G_m.add_edges_from(arr[:,:2])
        nx.set_edge_attributes(G_m, labels, 'patch_ids')

        dimers = [ cluster for cluster in nx.connected_components(G_m) if len(cluster) == 2]

        #dimer_attributes =
        N_dimer = 0
        for dimer in dimers:
        #    print()
            dimer = list(dimer)
            id1 = dimer[0]
            id2 = dimer[1]
            patch_tuple  = G_m.edges[id1,id2]['patch_ids']
            print("ptuple", patch_tuple, patch_dict_2p[cluster_type])
            if sorted(patch_tuple) == patch_dict_2p[cluster_type]:
                N_dimer +=2

        #N_dimer = len([ dimer for dimer in dimers if G_m.edges[dimer[0],dimer[1]]['patch_ids'] == patch_dict[cluster_type]])
        p_yield = N_dimer/N_particles
        return p_yield


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', type=str, default='yield_test.csv',
               help='output file name')
    parser.add_argument('-ctype', type=str, default='6-star',
               choices=['5-star', '6-star', 'box', 'open-box', 'open-6-star', 'open-5-star', 'open-jenga', 'l-shape', 'brick-shape', 'line-shape'],
               help='cluster type')
    parser.add_argument('-ptype', type=str, default='mouse',
                        choices=['mouse', 'manta'],
                        help = 'particle type')
    args = parser.parse_args()

    #-----------------------------

    filen = "patch_network.dat"
    filenpt = 'NPT_OUT.txt'
    dirlist = glob.glob("mu_*")
    pwd = os.getcwd()
    df = pd.DataFrame(columns=['mu','energy', 'topology', 'delta', 'delta_energy', 'mu2', 'p_yield','cluster_type'])

    for dir in dirlist:
        print(dir)
        numbers = re.findall(r"[-+]?\d*\.\d+|\d+", dir)
        mu = numbers[0]
        energy = numbers[1]
        delta = numbers[2]
        delta_energy = numbers[3]
        mu2 = numbers[4]
        topology = 'symm'

        if check_file('{}/{}'.format(dir,filen)) and check_file('{}/{}'.format(dir,filenpt)):
            arr = read_bonds('{}/{}'.format(dir, filen))
            bond_arr = arr[-1]

            N_particles = pd.read_csv('{}/{}'.format(dir, filenpt),
                                        header=None,
                                        delim_whitespace=True).values[-1,2]

            p_yield = calculate_yield(bond_arr,args.ctype,N_particles)

            df = pd.concat([df,pd.DataFrame({'mu': [mu],
                                                 'energy': [energy],
                                                 'topology': [topology],
                                                 'delta': [delta],
                                                 'delta_energy': [delta_energy],
                                                 'mu2' : [mu2],
                                                 'p_yield': [p_yield],
                                                 'cluster_type': [args.ctype]})])
            #except:
            #print("warning: couldn't calculate for {}".format(dir))


    yield_mean_std = df.groupby(['mu',
                                 'energy',
                                 'topology',
                                 'delta_energy',
                                 'mu2'
                                 'delta', 'cluster_type']).p_yield.agg(['mean','std']).reset_index()

    with open(args.o, 'a') as f:
        yield_mean_std.to_csv(f, header=None)
