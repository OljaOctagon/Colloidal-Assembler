import numpy as np
import networkx as nx
import pandas as pd
import glob
from collections import defaultdict

# get particle ids of largest connected domain (largest cluster)
def calculate_max_domain(arr):
    G=nx.Graph()
    G.add_edges_from(arr[:,:2])
    domains = list(nx.connected_components(G))
    domain_length = np.array([ len(domain) for domain in domains ])
    d_id = np.argmax(domain_length)
    particles_max_domain = np.array(list(domains[d_id]))
    return particles_max_domain


if __name__ == "__main__":

    print("hello")
    dirlists=glob.glob("mu_0.25*")
    delta_values = [0.2, 0.3, 0.4]
    fn="patch_network.dat"
    print(dirlists)
    counts= defaultdict(list)

    for directory in dirlists:
        print(directory)
        try:
            arr = pd.read_csv(directory+"/"+fn, header=None, delim_whitespace=True).values
            delta = float(directory.split("_")[-2])
            print(delta)
            l_particles = calculate_max_domain(arr)
            l_arr = np.array([ arr_i for arr_i in arr if arr_i[0] in l_particles or arr_i[1] in l_particles])

            for elem_i in l_arr[:,2:]:
                elem = list(elem_i)
                # parallel bonds, signature: 0
                if elem == [3,0] or elem == [0,3] or elem == [2,1] or elem == [1,2]:
                    counts[delta].append(0)
                # parallel off edge bonds, signature: 1
                if elem == [3,3] or elem == [2,2] or elem == [1,1] or elem == [0,0]:
                    counts[delta].append(1)
                # non-parallel edge bonds, signature: 2
                if elem == [2,0] or elem == [0,2] or elem == [3,1] or elem == [1,3]:
                    counts[delta].append(2)
                # parallel off edge bonds, signature: 3
                if elem == [0,1] or elem == [1,0] or elem == [3,2] or elem == [2,3]:
                    counts[delta].append(3)
        except:
            print("Couldn't read {}".format(directory))

        print(counts)

    import matplotlib.pyplot as plt
    from cycler import cycler

    fig,axes=plt.subplots(1,3, sharex=True, sharey=True, figsize=(20,5))
    cmap=plt.cm.viridis_r
    c = cycler('color', cmap(np.linspace(0,1,3)) )
    bins = np.arange(-0.5,4.5,1)
    
    axes[0].set_ylabel("percentage in largest cluster", size=14)
    for i,key, c in zip(range(len(counts)), delta_values, c):
        axes[i].hist(counts[key],bins, normed=True, facecolor=c['color'], label='$\Delta={}$'.format(key))
        axes[i].set_xticks([0,1,2,3])
        axes[i].set_xticklabels(['2-p','2-p-off', '2-np','2-np-off'], fontsize=12)
        axes[i].legend(fontsize=15)

    fig.text(0.5, 0.02, "bond type", ha='center', size=14)
    #plt.tight_layout()
    plt.savefig("bond_histogram.png", dpi=300)
    plt.show()
