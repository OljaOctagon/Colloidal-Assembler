import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df_crosslinker = pd.read_pickle("network_data_crosslinker.pickle")
df_polymer = pd.read_pickle("network_data_polymer.pickle")
df_bond_life_time = pd.read_pickle("bond_life_time.pickle")
df_corr_blcs = pd.read_pickle("corr_bond_life_cluster_size.pickle")

df_parameters = pd.read_csv("parameters.txt", delim_whitespace=True)
# Index(['ID', 'N', 'kb', 'rho', 'plink', 'epsilon', 'pstart'], dtype='object')

# Plot all crosslinker observables



# Plot all polymer observables
# columns df: id, time, domain_lengths, largest_domain, degree_sequence, std_degree, node_connectivity

# percolation averaged over time: fraction of largest domain.
# by epsilon, rho

# first: combine df_param with df_crosslinker

df_parameters.rename(columns = {'ID':'id'}, inplace=True)
df = pd.merge(df_crosslinker, df_parameters, on ='id')
nparticles = 340
fig,ax=plt.subplots()

mean_percolation = np.empty([])

for epsilon in Epsilons:
    df_1 = f[df.epsilon == epsilon]
    mean_epsi = np.mean(df_1.largest_domain.values)/nparticles


