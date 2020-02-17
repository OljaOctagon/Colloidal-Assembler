import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os.path
import networkx as nx
import argparse
import glob
import re

def check_file(fname):
    try:
        open(fname,"r")
        return 1
    except IOError:
        print("Error: File doesn't seem to exist.")
        return 0


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

def plot_histogram_horizontal(dg):

    colors = ["#236AB9", "#ff0066", "#9900cc","#ffff00"]

    error_config = dict(ecolor='black', lw=2, capsize=6, capthick=2, alpha=0.8)
    width = 0.1

    delta_range = dg.delta.unique()
    ind = np.arange(len(delta_range)) + 0.2  - width/2.
    lw=3
    plt.xticks(ind+width/2.)

    min_delta = np.min(delta_range)
    max_delta = np.max(delta_range)

    for mu in dg.mu.unique():
        for topo in dg.topology.unique():
            for energy in dg.energy.unique():
                df_sub = dg[(dg.mu == mu) & (dg.topology == topo) & (dg.energy == energy)]
                arr = df_sub[['n_liquid', 'n_micell_1', 'n_micell_2','n_chain_and_loop']].values
                brr = df_sub[['packing_fraction']].values

                fig, ax =  plt.subplots(figsize=(15,11))

                plt.xlim((min_delta-width/2.,max_delta+width/2.))
                plt.ylim((0,1))
                plt.tick_params(axis='both', which='major', labelsize=28)
                plt.xlabel("$\\Delta $", size=38)
                plt.ylabel("$P_{c}$", size=38)


                plt.bar(ind+0.5, arr[:,0],
                        width, lw=lw,yerr=arr[:,1],
                        error_kw=error_config, color=colors[0], alpha=0.7)

                plt.bar(ind+0.5, arr[:,2],
                        width, lw=lw,yerr=arr[:,3],
                        error_kw=error_config, color=colors[1]
                        ,bottom=arr[:,0],alpha=0.7)

                plt.bar(ind+0.5, arr[:,4],
                        width, lw=lw,yerr=arr[:,5],
                        error_kw=error_config, color=colors[2],
                        bottom=arr[:,0]+arr[:,2], alpha=0.9)

                plt.bar(ind+0.5, arr[:,6], width,
                        lw=lw,yerr=arr[:,7],
                        error_kw=error_config, color=colors[3],
                        bottom=arr[:,0]+arr[:,2]+arr[:,4], alpha=0.7)

                for ei in range(len(delta_range)):
                    ax.text(ind[ei], 1.07, "$\\psi= {}$".format(round(brr[ei,0],3)), fontsize=24)
                    ax.text(ind[ei], 1.03, "$\pm {}$".format(round(brr[ei,1],3)), fontsize=24)

                plt.savefig("results_horizontal/mu_{}_{}_pachpos_{}_percent.pdf".format(mu,topo,delta))

def plot_histogram(dg):
    #color1 = #8da19e
    colors = ["#061e67", "#fd4949", "#fd4949","#ffc81d"]

    error_config = dict(ecolor='black', lw=2, capsize=6, capthick=2, alpha=0.8)
    width = 1

    energy_range = dg.energy.unique()
    ind = np.arange(len(energy_range)) + 5 - width/2.
    lw=3
    plt.xticks(ind+width/2.)

    min_energy = np.min(energy_range)
    max_energy = np.max(energy_range) #061e67


    for mu in dg.mu.unique():
        for topo in dg.topology.unique():
            for delta in dg.delta.unique():
                df_sub = dg[(dg.mu == mu) & (dg.topology == topo) & (dg.delta == delta)]
                arr = df_sub[['n_liquid', 'n_micell_1', 'n_micel

    if not os.path.isdir("./phase_histograms"):
        os.mkdir("./phase_histograms")l_2','n_chain_and_loop']].values
                brr = df_sub[['packing_fraction']].values

                fig, ax =  plt.subplots(figsize=(15,11))

                plt.xlim((min_energy-width/2.,max_energy+width/2.))
                plt.ylim((0,1))
                plt.tick_params(axis='both', which='major', labelsize=28)
                plt.xlabel("$\\epsilon [k_{B}T]$", size=38)
                plt.ylabel("$P_{c}$", size=38)


                plt.bar(ind+0.5, arr[:,0],
                        width, lw=lw,yerr=arr[:,1],
                        error_kw=error_config, color=colors[0], alpha=0.7)

                plt.bar(ind+0.5, arr[:,2],
                        width, lw=lw,yerr=arr[:,3],
                        error_kw=error_config, color=colors[1]
                        ,bottom=arr[:,0],alpha=0.7)

                plt.bar(ind+0.5, arr[:,4],
                        width, lw=lw,yerr=arr[:,5],
                        error_kw=error_config, color=colors[2],
                        bottom=arr[:,0]+arr[:,2], alpha=0.9)

                plt.bar(ind+0.5, arr[:,6], width,
                        lw=lw,yerr=arr[:,7],
                        error_kw=error_config, color=colors[3],
                        bottom=arr[:,0]+arr[:,2]+arr[:,4], alpha=0.7)

                for ei in range(len(energy_range)):
                    ax.text(ind[ei], 1.07, "$\\psi= {}$".format(round(brr[ei,0],3)), fontsize=24)
                    ax.text(ind[ei], 1.03, "$\pm {}$".format(round(brr[ei,1],3)), fontsize=24)

                plt.savefig("phase_histograms/mu_{}_{}_pachpos_{}_percent.pdf".format(mu,topo,delta))


if __name__ == '__main__':

    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-tsizes', nargs='+')
    parser.add_argument('-o', type=str, default='cluster_distribution_with_micell.csv')
    args  = parser.parse_args()


    args.tsizes = list(map(int,args.tsizes))

    assert (len(args.tsizes) <= 2), "too many cluster types, can be 1 or 2."


    filen = "patch_network.dat"
    filenpt = 'NPT_OUT.txt'
    dirlist = glob.glob("mu_*")
    pwd = os.getcwd()

    features = ['mu', 'energy', 'topology',
                'delta', 'packing_fraction',
                'n_liquid', 'n_micell_1', 'n_micell_2', 'n_chain_and_loop']

    df = pd.DataFrame(columns=features)

    if not os.path.isdir("./phase_histograms"):
        os.mkdir("./phase_histograms")

    for dir in dirlist:
        numbers = re.findall(r"[-+]?\d*\.\d+|\d+", dir)
        mu = numbers[0]
        energy = numbers[1]
        delta = numbers[2]
        topology = 'symm'

        if check_file('{}/{}'.format(dir,filen)) and check_file('{}/{}'.format(dir,filenpt)):
            print(dir,filen)
            network_arr = read_bonds('{}/{}'.format(dir, filen))
            arr = network_arr[-1]

            arr_npt = pd.read_csv('{}/{}'.format(dir, filenpt),
                                    header=None,
                                    delim_whitespace=True).values

            N_particles = arr_npt[-1,2]
            packing_fraction = arr_npt[-1,4]


            G = nx.Graph()
            G.add_edges_from(arr[:,:2])
            clusters = list(nx.connected_components(G))
            cluster_sizes = np.array([ len(cluster) for cluster in clusters ])


            # how many clusters are there
            N_clusters = len(cluster_sizes)

            # how many liquid clusters are there
            # liquid = all clusters smaller than tsizes[0]

            n_ci = lambda x,k:len(x[x==k])
            N_liquid=0
            for i in range(1, args.tsizes[0]):
                N_liquid += i*n_ci(cluster_sizes,i)

            n_liquid = N_liquid/N_particles

            # how many micells are there
            loops = list(nx.cycle_basis(G))
            loop_sizes = np.array([ len(loop) for loop in loops])

            N_micell = {}
            n_micell = []

            for mtype in args.tsizes:
                N_micell[mtype] = n_ci(loop_sizes,mtype)
                n_micell.append((mtype*N_micell[mtype])/N_particles)

            # how many chains are there
            N_chain_and_loop = 0

            for i in range(args.tsizes[0], np.max(cluster_sizes)+1):
                if i in args.tsizes:
                    N_chain_and_loop += i*(n_ci(cluster_sizes,i)- N_micell[i])
                else:
                    N_chain_and_loop += i*(n_ci(cluster_sizes,i))

            n_chain_and_loop = N_chain_and_loop/N_particles


            n_micell_2 = n_micell[0]
            n_micell_1 = 0.00001
            if len(n_micell)==2:
                n_micell_1 = n_micell[1]

            n_particles_size_1 = 1 - n_liquid - n_chain_and_loop - n_micell_1 - n_micell_2
            n_liquid = n_liquid + n_particles_size_1

            df = pd.concat(
                [df,pd.DataFrame({'mu': [float(mu)],
                                'energy': [float(energy)],
                                'topology': [topology],
                                'delta': [float(delta)],
                                'packing_fraction': [packing_fraction],
                                'n_liquid': [n_liquid],
                                'n_micell_1': [n_micell_1],
                                'n_micell_2': [n_micell_2],
                                'n_chain_and_loop' : [n_chain_and_loop]})])


    # aggregate mean and std
    group_para  = ['mu','energy','delta','topology']
    df_grouped = df.groupby(group_para).agg(['mean','std']).reset_index()

    dg = df.groupby(group_para).agg('mean').reset_index()
    # write dataframe
    dg.to_csv("phase_histograms/{}".format(args.o), index=False)
    # plot histogram
    plot_histogram(df_grouped)
    #plot_histogram_horizontal(df_grouped)
