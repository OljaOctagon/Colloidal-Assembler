import networkx as nx
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
# This script calculates the domain size of parallel and non-parallel bonds of the largest cluster in the rhombi system.
# We interpret the largest cluster as two networks. One with all bondings of 1 (parallel) and one with all bondings of -1 (non-parallel).
# Calculating the connected components gives the size of all domains and number of domains.
# We then calculate the domain size distribution, the mean, and the variance. 

def calculate_network_domains(arr, p_index):
    G = nx.Graph()
    arr_p = arr[arr[:,2]==p_index]
    percentage=len(arr_p)/len(arr)

    G.add_edges_from(arr_p[:,:2])
    domains = list(nx.connected_components(G))
    domain_sizes = [ len(domain) for domain in domains]
    N_domains=len(domain_sizes)
    fig, ax = plt.subplots()
    plt.hist(domain_sizes, bins=N_domains, alpha=0.3)
    plt.savefig('domain_sizes_{}.png'.format(p_index),dpi=300)
    mean=np.mean(domain_sizes)
    std = np.std(domain_sizes)
    cluster_size=arr[0,3]
    return mean, std, N_domains, percentage, cluster_size

def domain_size(file):
    # This is an array of the format particle i, particle j, bonding (-1 or 1), cluster size
    arr = pd.read_csv(file, delim_whitespace=True, header=None).values
    # calculate for parallel
    p_mean, p_std, Ndp, pp, cs = calculate_network_domains(arr, 1)
    # calculate for non-parallel
    np_mean, np_std, Ndnp, pnp, cs = calculate_network_domains(arr, -1)

    Np_frac = Ndp/(Ndp+Ndnp)
    Nnp_frac = Ndnp/(Ndp+Ndnp)

    with open("domain_sizes_results.dat", 'w') as fhandle:
        fhandle.write("{} {} {} {} {} {}\n".format(p_mean, p_std, Np_frac, pp, 1, cs))
        fhandle.write("{} {} {} {} {} {}\n".format(np_mean, np_std, Nnp_frac, pnp, -1, cs))

if __name__ == "__main__":

    # bond_domains.dat contains information about bonding in the largest cluster
    # format: particle i, particle j, bonding (-1 for non-parallel, 1 for parallel)
    filen="bond_domains.dat"
    domain_size(filen)
