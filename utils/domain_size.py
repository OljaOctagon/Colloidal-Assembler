import networkx as nx
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
# This script calculates the domain size of parallel and non-parallel bonds of the largest cluster in the rhombi system.
# We interpret the largest cluster as two networks. One with all bondings of 1 (parallel) and one with all bondings of -1 (non-parallel).
# The number of connected components gives the domain sizes. We then calculate the distribution, the mean, and the variance. 

# This is ana array of the format particle i, particle j, bonding (-1 or 1), cluster size

def calculate_network_domains(arr, p_index):
    G = nx.Graph()
    arr_p = arr[arr[:,2]==p_index]
    G.add_edges_from(arr_p[:,:2])
    domains = list(nx.connected_components(G))
    domain_sizes = [ len(domain) for domain in domains]
    plt.hist(domain_sizes)
    plt.savefig('domain_sizes_{}.png'.format('p_index'),dpi=300)
    mean=np.mean(domain_sizes)
    std = np.std(domain_sizes)
    cluster_size=arr[0,3]
    return mean, std, cluster_size

def domain_size(file):
    # This is an array of the format particle i, particle j, bonding (-1 or 1), cluster size
    arr = pd.read_csv(file, delim_whitespace=True, header=None).values
    # calculate for parallel
    p_mean, p_std, cs = calculate_network_domains(arr, 1)
    # calculate for non-parallel
    np_mean, np_std, cs = calculate_network_domains(arr, -1)

    with open("domain_sizes_results.dat", 'w') as fhandle:
        fhandle.write("{} {} {} {}\n".format(p_mean, p_std, 1, cs))
        fhandle.write("{} {} {} {}\n".format(np_mean, np_std, -1, cs))

if __name__ == "__main__":
    filen="bond_domains.dat"
    domain_size(filen)
