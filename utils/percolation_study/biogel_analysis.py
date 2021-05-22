import glob
import pandas as pd
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.style as style
import matplotlib as mpl
import networkx as nx
import gel_tools as gt
import argparse
import random
import pickle


style.use('seaborn-ticks')
mpl.rcParams['font.family'] = "sans-serif"
plt.rcParams['axes.axisbelow'] = True

plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['figure.titlesize'] = 15

blue_c ='#9999FF'
red_c ='#FF9999'
purple_c='#C17DCB'
green_c='#61DCB7'

from collections import defaultdict
from functools import partial


def calculate_corr_blcs(arr, nbonds, nparticles, Ltime, bond_sequences, file_id, df_bond_life_cluster_size):
   
    a = np.arange(1, nparticles + 1)
    b = np.arange(1, nbonds + 1)
    f_to_number = lambda tup: (tup[0] - 1) * nbonds + (tup[1] - 1)
    arr_combinations = np.array([np.meshgrid(a, b)]).T.reshape(-1, 2)
    tuple_combinations = [(i, j) for [i, j] in arr_combinations]
    map_to_number = list(map(f_to_number, tuple_combinations))
    dict_to_number = dict(zip(tuple_combinations, map_to_number))

    list_keys = list(bond_sequences.keys())

    Npairs=10000
    Pairs = [ random.choice(list_keys) for i in range(Npairs)]
    # get bond life times

    from itertools import groupby

    def len_iter(items):
        return sum(1 for _ in items)

    def consecutive_one(data):
        return [len_iter(run) for val, run in groupby(data) if val]

    life_times = defaultdict(list)
    for pair in Pairs:
        life_times[pair] = consecutive_one(bond_sequences[pair])

    cluster_sizes = defaultdict(list)
    average_cluster_sizes = defaultdict(list)

    for time_i in range(1,Ltime+1):

        arr_i = arr[arr[:, 0] == time_i]

        connections = np.array([
            [float(dict_to_number[(pid1, lid1)]),
             float(dict_to_number[(pid2, lid2)])] for [pid1, lid1, pid2, lid2] in arr_i[:, 1:5]])

        G = nx.Graph()
        G.add_edges_from(connections)

        for pair_i in Pairs:
            # Add new cluster size if bonded at time_i
            if bond_sequences[pair_i][time_i] == 1:
                entry = len(list(nx.node_connected_component(G, pair_i[0])))
                cluster_sizes[pair_i].append(entry)

            # Bond breaks: calculate average cluster size
            if bond_sequences[pair_i][time_i] == 0 and bond_sequences[pair_i][time_i-1] == 1:
                av_size = np.mean(cluster_sizes[pair_i])
                average_cluster_sizes[pair_i].append(av_size)

            # Add the last size
            if time_i == Ltime and bond_sequences[pair_i][time_i] == 1:
                av_size = np.mean(cluster_sizes[pair_i])
                average_cluster_sizes[pair_i].append(av_size)

    # Combine life_times and average_cluster_sizes to one 2D array
    arr_combined = []
    for pair_i in Pairs:
        for ni in range(len(life_times[pair_i])):
            arr_combined.append([life_times[pair_i][ni], average_cluster_sizes[pair_i][ni]])

    arr_combined = np.array(arr_combined)

    new_results = pd.DataFrame(arr_combined,columns=['life_time','av_cluster_size'])
    new_results['id'] = np.repeat(file_id,len(arr_combined))
    df_bond_life_cluster_size = df_bond_life_cluster_size.append(new_results, ignore_index=True)
    return df_bond_life_cluster_size


def calculate_crosslinker_network_observables(arr, Ltime, file_id, nbonds, nparticles, df_crosslinker):

    for time_i in range(1, Ltime + 1):
        # get connections per polymer per time_i
        arr_i = arr[arr[:, 0] == time_i]

        a = np.arange(1, nparticles + 1)
        b = np.arange(1, nbonds + 1)
        f_to_number = lambda tup: (tup[0] - 1) * nbonds + (tup[1] - 1)
        arr_combinations = np.array([np.meshgrid(a, b)]).T.reshape(-1, 2)
        tuple_combinations = [(i, j) for [i, j] in arr_combinations]
        map_to_number = list(map(f_to_number, tuple_combinations))
        dict_to_number = dict(zip(tuple_combinations, map_to_number))

        connections = np.array([
            [dict_to_number[(pid1, lid1)],
             dict_to_number[(pid2, lid2)]] for [pid1, lid1, pid2, lid2] in arr_i[:, 1:5]])

        G = nx.Graph()
        G.add_edges_from(connections)

        # initalize new results
        new_results = {}
        # Calculate domain lengths
        new_results['id'] = file_id
        new_results['time'] = time_i
        domain_lengths = gt.get_domain_lengths(G)
        new_results['domain_lengths'] = domain_lengths
        new_results['largest_domain'] = np.max(domain_lengths)

        # Calculate degree
        degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
        new_results['degree_sequence'] = degree_sequence
        new_results['mean_degree'] = np.mean(degree_sequence)
        new_results['std_degree'] = np.std(degree_sequence)
        new_results['node_connectivity'] = nx.node_connectivity(G)

        df_crosslinker = df_crosslinker.append(new_results, ignore_index=True)

    return df_crosslinker

def calculate_polymer_network_observables(arr, Ltime, file_id, df_polymer):

    for time_i in range(1, Ltime + 1):
        # get connections per polymer per time_i
        arr_i = arr[arr[:, 0] == time_i]

        connections = np.column_stack((arr_i[:, 1], arr_i[:, 3]))
        G = nx.Graph()
        G.add_edges_from(connections)

        # initalize new results
        new_results = {}
        # Calculate domain lengths
        new_results['id'] = file_id
        new_results['time'] = time_i
        domain_lengths = gt.get_domain_lengths(G)
        new_results['domain_lengths'] = domain_lengths
        new_results['largest_domain'] = np.max(domain_lengths)

        # Calculate degree
        degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
        new_results['degree_sequence'] = degree_sequence
        new_results['mean_degree'] = np.mean(degree_sequence)
        new_results['std_degree'] = np.std(degree_sequence)
        new_results['node_connectivity'] = nx.node_connectivity(G)

        df_polymer = df_polymer.append(new_results, ignore_index=True)

    return df_polymer



def calculate_bond_lifetime(arr, nbonds, nparticles,Ltime, file_id, df_bl_histogram):

    a = np.arange(1, nparticles + 1)
    b = np.arange(1, nbonds + 1)

    f_to_number = lambda tup: (tup[0] - 1) * nbonds + (tup[1] - 1)
    #f_to_tuple = lambda c: ((c // nbonds) + 1, (c % nbonds) + 1)

    arr_combinations = np.array([np.meshgrid(a, b)]).T.reshape(-1, 2)
    tuple_combinations = [(i, j) for [i, j] in arr_combinations]
    map_to_number = list(map(f_to_number, tuple_combinations))
    dict_to_number = dict(zip(tuple_combinations, map_to_number))
    #dict_to_tuple = dict(zip(map_to_number, tuple_combinations))

    arr_id = np.array([
        [time, dict_to_number[(pid1, lid1)],
         dict_to_number[(pid2, lid2)]] for [time, pid1, lid1, pid2, lid2] in arr[:, :5]])

    connections = [arr_id[arr_id[:, 0] == time][:, 1:] for time in range(1, Ltime + 1)]
    # initialize bond sequences, keys: bond pairs, values: time series bonded states.
    # 1: is bonded 0 not bonded
    bond_sequences = defaultdict(partial(np.zeros, Ltime+1))
    print("get bond info")

    for time in range(Ltime):
        for [i, j] in connections[time]:
            bond_sequences[(i, j)][time+1] = 1

    from itertools import groupby

    def len_iter(items):
        return sum(1 for _ in items)

    def consecutive_one(data):
        return [len_iter(run) for val, run in groupby(data) if val]

    print("calculate lifetime")
    # calculate bond life times
    lifetimes = []
    arr_formed = []
    arr_reformed = []
    auto_corr_heat = np.empty((0,2))

    fig,ax = plt.subplots()

    for bond in bond_sequences.keys():
        i = bond[0]
        j = bond[1]
        # get life time
        seq = bond_sequences[i, j]
        lifetimes.extend(consecutive_one(seq))

        diffx = seq[1:]-seq[:-1]
        # get reformed bonds
        nformed = np.count_nonzero( diffx == 1 )
        arr_formed.append(nformed + seq[0])
        arr_reformed.append(nformed - 1)

        # get normalized autocorr of seq
        def autocorr(x):
            x = x - np.mean(x)
            result = np.correlate(x, x, mode='full')
            result = result[int(np.floor(result.size/2)):]
            rmax=np.max(result)
            result = result/rmax
            return result

        ac_seq_y = autocorr(seq)
        ac_seq_x = np.linspace(0,Ltime,Ltime+1)
        ac_seq = np.column_stack((ac_seq_x, ac_seq_y))
        a=np.random.randint(0,500)
        if a==0:
            ax.plot(ac_seq_x,ac_seq_y,lw=1)

        auto_corr_heat = np.concatenate((auto_corr_heat,ac_seq))

    plt.xlabel("Time")
    plt.ylabel("Autocorrelation of bond life time")
    plt.savefig("data_ana/{}/auto_corr_heat_seq.pdf".format(file_id))
    lifetimes = np.array(lifetimes)
    max_lifetime = np.max(lifetimes)

    Nformed = np.sum(np.array(arr_formed))
    Nreformed = np.sum(np.array(arr_reformed))
    percent_reformed = Nreformed/Nformed

    print("make reformed bond histogram")

    rval, bins, patches = plt.hist(arr_reformed,
                                   edgecolor=purple_c, bins=(np.linspace(0, Ltime, Ltime + 1) + 0.1),
                                   lw=2, alpha=0.5,
                                   density=True)
    fig, ax = plt.subplots()
    plt.yscale("log")
    x = bins[1:]

    # fit exponential
    def monoExp(x, l, t, ):
        return l * np.exp(-t * x)

    p0 = (1, 1)  # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(monoExp, x, rval, p0)
    l, t = params
    plt.plot(x, rval, c=red_c, linestyle='dotted', label='data')
    plt.plot(x, monoExp(x, l, t), c=blue_c, linestyle='--', label='biexp. fit')
    plt.legend(loc='best')

    plt.xlabel("reformed bonds")
    plt.ylabel("P")
    plt.ylim(0.001, 0.2)
    plt.xlim(0,150)
    plt.tight_layout()
    plt.savefig("data_ana/{}/reformed_bonds_fit.pdf".format(file_id))
    arr_data = np.array([l, t, percent_reformed])
    np.savetxt("data_ana/{}/param_reformed_bonds.dat".format(file_id), arr_data, delimiter=',', newline='\n')
    np.savetxt("data_ana/{}/reformed_bonds.dat".format(file_id), rval, newline='\n')


    print("make life time histogram")
    fig, ax = plt.subplots()
    plt.xlabel("bond life time")
    plt.ylabel("P")
    nval, bins, patches = plt.hist(lifetimes,
                                edgecolor='k', bins=(np.linspace(0, Ltime, Ltime + 1) + 0.1),
                                lw=2, alpha=1,
                                density=True)

    fig, ax = plt.subplots()
    plt.yscale("log")
    x = bins[1:]

    # fit exponential
    def monoExp(x, l, t, ):
        return l * np.exp(-t * x)

    p0 = (1, 1)  # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(monoExp, x, nval, p0)
    l, t = params

    plt.plot(x, nval, c=purple_c, linestyle='dotted', label='data')
    plt.plot(x, monoExp(x, l, t), c=blue_c, linestyle='--', label='biexp. fit')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.ylim(0.00001, 1)
    plt.savefig("data_ana/{}/bond_lifetime_fit.pdf".format(file_id))

    def cumulative_monoExp(xs, l, t):
        return 1 - l * np.exp(-t * xs)

    #H,x_edges, y_edges = np.histogram2d(auto_corr_heat[:,0], auto_corr_heat[:,1],
    #                                        bins=[500,100], density=True)

    np.savetxt("data_ana/{}/autocorr_lifetimes.dat".format(file_id),auto_corr_heat,newline='\n')
    #fig,ax  = plt.subplots()
    #sns.heatmap(H)
    #plt.savefig("data_ana/{}/autocorr_lifetimes_heatmap.pdf".format(file_id))

    arr_data = np.array([l, t, max_lifetime, percent_reformed])
    np.savetxt("data_ana/{}/param_bond_lifetime.dat".format(file_id), arr_data, delimiter=',', newline='\n')
    np.savetxt("data_ana/{}/bond_lifetime.dat".format(file_id), nval, newline='\n')

    new_results = pd.DataFrame(nval, columns=['p_life_time'])
    new_results['id'] = np.repeat(file_id, len(nval))
    new_results['time'] = np.linspace(1,2500,2500)

    df_bl_histogram = df_bl_histogram.append(new_results, ignore_index=True)

    return bond_sequences, df_bl_histogram

if __name__ == '__main__':

    files_list = glob.glob("*/pdf1/*link.dat")
    df_params = pd.read_csv("parameters.txt", delim_whitespace=True,dtype=str)
    print(df_params.head(5))
    nparticles = 340
    Ltime = 2500

    parser = argparse.ArgumentParser()
    parser.add_argument("-bond_life_time", type=int, help="bond life time", default=0)
    parser.add_argument("-polymer_network", type=int, help="calculate polymer network properties", default=0)
    parser.add_argument("-crosslinker_network", type=int, help="calculate crosslinker network properties", default=0)

    args = parser.parse_args()

    columns = ['id', 'time', 'mean_degree', 'std_degree',
               'degree_sequence', 'largest_domain', 'domain_lengths', 'node_connectivity']
    df_crosslinker = pd.DataFrame(columns=columns)

    columns = ['id', 'time', 'mean_degree', 'std_degree',
               'degree_sequence', 'largest_domain', 'domain_lengths', 'node_connectivity']
    df_polymer = pd.DataFrame(columns=columns)

    columns = ['id','life_time','av_cluster_size']
    df_bond_life_cluster_size = pd.DataFrame(columns=columns)

    columns = ['id','time','p_life_time']
    df_bl_histogram = pd.DataFrame(columns=columns)

    for file_i in files_list:
        print("read gel")
        file_id = file_i.split('/')[0]
        print(file_id)
        arr = pd.read_csv(file_i, delim_whitespace=True).values
        nbonds = int(df_params[df_params.ID == file_id].N.values[0])

        print("prepare data")
        # throw away links between crosslinkers on same polymer
        arr = arr[arr[:, 1] != arr[:, 3]]

        print(args.bond_life_time)
        if int(args.bond_life_time) == 0:
            print("calculate bond life time")
            bond_sequences, df_bl_histogram = calculate_bond_lifetime(arr, nbonds, nparticles, Ltime, file_id, df_bl_histogram)

            #
            #with open('bond_sequences.pickle', 'wb') as handle:
            #    pickle.dump(bond_sequences, handle, protocol=pickle.HIGHEST_PROTOCOL)

            #with open('bond_sequences.pickle', 'rb') as handle:
            #    bond_sequences = pickle.load(handle)

            print("calculate bond life time cluster size correlation")
            df_bond_life_cluster_size = calculate_corr_blcs(arr, nbonds, nparticles, Ltime,
                                                            bond_sequences, file_id,
                                                            df_bond_life_cluster_size)

        if int(args.polymer_network) == 0:
            print("calculate polymer network observalbes per time ")
            df_polymer = calculate_polymer_network_observables(arr, Ltime, file_id, df_polymer)


        if int(args.crosslinker_network) == 0:
            print("calculate crosslinker network observalbes per time")

            df_crosslinker = calculate_crosslinker_network_observables(arr,Ltime,
                                                                       file_id, nbonds,
                                                                       nparticles, df_crosslinker)

    df_bl_histogram.to_pickle("data_ana/bond_life_time.pickle")
    df_crosslinker.to_pickle("data_ana/network_data_crosslinker.pickle")
    df_polymer.to_pickle("data_ana/network_data_polymer.pickle")
    df_bond_life_cluster_size.to_pickle("data_ana/corr_bond_life_cluster_size.pickle")