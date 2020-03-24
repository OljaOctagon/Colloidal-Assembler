
import numpy as np
import pandas as pd
import argparse
import gel_tools as gt
import networkx as nx
import glob
import os

# TODO:
# Find a way to find "crystal loop communities"
# find more measures for network topology

def sub_domain_size(connections_j, bond_dict, indicator):
    G_sub=nx.Graph()
    sub_connections_j = [ item[:2] for item in connections_j if bond_dict[tuple(sorted(item[2:]))] == indicator  ]
    G_sub.add_edges_from(sub_connections_j)
    sub_domain_lengths = gt.get_domain_lengths(G_sub)
    largest_sub_domain = np.max(sub_domain_lengths)

    return sub_domain_lengths, largest_sub_domain 

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str)
    parser.add_argument('-f', type=str)
    args = parser.parse_args()

    N_particles = 1500 
    valence=4
    T=0.15
    phi=0.5
    run_id="dma-as1_r001"

    # 1: paralelel
    # 0: non parallel
    bond_dict = {(0,0): 1, (0,1): 0, (1,1):1,
                (2,2): 1, (2,3): 0, (3,3):1,
                (0,2): 0, (0,3): 1, (1,2):1,
                (1,3):0}

    # 1: satisfied
    # 0: not satisfied 
    satisfied_bond_dict = {(0,0): 1, (0,1): 1, (1,1):1,
                (2,2): 1, (2,3): 1, (3,3):1,
                (0,2): 0, (0,3): 0, (1,2):0,
                (1,3):0}


    # read in connections, format: particle id1, particle id2, patch id1, patch id2
    dir_name = args.d
    f_name = args.f

    # get all check point values and sort them
    checkpoints= glob.glob("{}Box*.bin".format(dir_name))
    check_point_values = np.sort(
    [ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])

    pn_file = dir_name+f_name
    # network_arr format: network_arr.shape = ( frame_i, bond_rows_frame_i )
    connections = gt.read_bonds(pn_file)

    fe_name="Energy.dat"
    energy_file = dir_name+fe_name
    df_energy = pd.read_csv(energy_file, header=None, delim_whitespace=True)
    df_energy.columns = ['time','energy']

    results_dir = 'gel_observables'
    # make frame directory if it doesn't exist 
    if not os.path.isdir(dir_name+results_dir):
        os.mkdir(dir_name+results_dir)

    # TODO  simplify graph by combining all 3 loops
    columns = ['run_id',
               'topology',
               'delta','T',
               'phi','N_particles',
               'time',
               'domain_lengths',
               'p_domain_lengths',
               'np_domain_lengths',
               'degrees',
               'mean_degree',
               'std_degree',
               'pb',
               'node_connectivity',
               'N3_loop_percent',
               'largest_domain',
               'largest_pdomain',
               'largest_np_domain',
               'pbond_percent',
               'npbond_percent',
               'energy']

    df = pd.DataFrame(columns=columns)

    for j,val in enumerate(check_point_values[1:]):

         print("Checkpoint ....",val)
         # Initialize new results
         new_results = {}

         new_results['run_id'] = '0001'
         new_results['topology'] = 'dma-as1'
         new_results['delta'] = 0.2
         new_results['T'] = 0.15
         new_results['time'] = val
         new_results['N_particles'] = 1500
         new_results['energy'] = df_energy[df_energy.time==val].energy.values[0]


         # Make a graph for time j: 
         connections_j = connections[j]
         G=nx.Graph()
         G.add_edges_from(connections_j[:,:2])

         # Calculate domain lengths 
         domain_lengths = gt.get_domain_lengths(G)
         largest_domain = np.max(domain_lengths)

         new_results['domain_lengths'] = domain_lengths
         new_results['largest_domain'] = largest_domain

         # Calculate degree
         degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence

         new_results['degrees'] = degree_sequence
         new_results['mean_degree'] = np.mean(degree_sequence)
         new_results['std_degree']  = np.std(degree_sequence)

         # Caclulate fraction of number of bonds pbonds
         Nb_max = (N_particles * valence)/2
         Nbonds = G.number_of_edges()
         pbonds = Nbonds/Nb_max

         new_results['pb'] = pbonds

         # Node connectivity
         new_results['node_connectivity'] = nx.node_connectivity(G)

         # Percent pbonds/npbonds
         N_pbonds = np.sum([ bond_dict[tuple(sorted(l))] for l in connections_j[:,2:] ])
         pbond_percent = N_pbonds/Nbonds
         npbond_percent = 1 - pbond_percent

         new_results['pbond_percent'] = pbond_percent
         new_results['npbond_percent'] = npbond_percent

         # Caclulate sizes of p/np bonded domains
         p_domain_lengths, largest_p_domain = sub_domain_size(connections_j,bond_dict, 1)
         np_domain_lengths, largest_np_domain = sub_domain_size(connections_j, bond_dict, 0)

         new_results['p_domain_lengths'] = p_domain_lengths
         new_results['np_domain_lengths'] = np_domain_lengths
         new_results['largest_p_domain'] = largest_p_domain
         new_results['largest_np_domain'] = largest_np_domain

         N3_loops = gt.find_cliques(G)
         N3_loop_percent = (3*N3_loops)/N_particles
         new_results['N3_loop_percent'] = N3_loop_percent

         df = df.append(new_results, ignore_index=True)

    df.to_pickle("{}{}/network_data.pickle".format(dir_name, results_dir))
