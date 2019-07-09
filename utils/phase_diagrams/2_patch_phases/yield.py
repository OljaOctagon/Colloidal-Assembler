import numpy as np
import pandas as pd
import glob
import re
import os
import argparse 
import networkx as nx
import itertools


cluster_type_dict = {
    '6-star': 6,
    '5-star': 5,
    'box': 3,
    'open-box':3,
    'open-6-star':6,
    'open-5-star': 5}

patch_dict = {
    'mouse': [0,2],
    'manta': [0,1]
}

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

    # Filter out mouse np bonds and construct graph
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


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', type=str, default='yield.csv',
               help='output file name')
    parser.add_argument('-ctype', type=str, default='6-star',
               choices=['5-star', '6-star', 'box', 'open-box', 'open-6-star', 'open-5-star'],
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
    df = pd.DataFrame(columns=['mu','energy', 'topology', 'delta', 'p_yield','cluster_type'])

    for dir in dirlist:
        print(dir)
        numbers = re.findall(r"[-+]?\d*\.\d+|\d+", dir)
        mu = numbers[0]
        energy = numbers[1]
        delta = numbers[2]
        topology = 'symm'
        p = 'Asymm'
        if re.search(p, dir):
            topology = 'Asymm'

        try:
            bond_arr = pd.read_csv('{}/{}'.format(dir, filen),
                                header=None,
                                delim_whitespace=True).values

            N_particles = pd.read_csv('{}/{}'.format(dir, filenpt),
                                    header=None,
                                    delim_whitespace=True).values[-1,2]

            print(N_particles)

            p_yield = calculate_yield(bond_arr,args.ctype,N_particles)


            df = pd.concat([df,pd.DataFrame({'mu': [mu],
                                             'energy': [energy],
                                             'topology': [topology],
                                             'delta': [delta],
                                             'p_yield': [p_yield],
                                             'cluster_type': [args.ctype]})])


        except:
          print("warning: couldn't calculate for {}".format(dir))


    yield_mean_std = df.groupby(['mu',
                                 'energy',
                                 'topology',
                                 'delta', 'cluster_type']).p_yield.agg(['mean','std']).reset_index()

    print(args.o)
    with open(args.o, 'a') as f:
        yield_mean_std.to_csv(f, header=None)
    
