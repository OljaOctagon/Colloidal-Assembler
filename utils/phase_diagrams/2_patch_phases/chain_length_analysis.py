# coding: utf-8
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os.path
import seaborn as sns
sns.set(style='whitegrid')

def check_file(fname):
    try:
        open(fname,"r")
        return 1
    except IOError:
        print("Error: File doesn't seem to exist.")
        return 0

Chemical_Potentials = [0.3,0.4]
Energies = [8.2, 7.2, 6.2, 5.2]
Patch_Topologies = ['Asymm']
patch_A = [0.3,0.5,0.7] 
patch_B = [0.3,0.5,0.7] 
phi_binary = [0.5,0.25,0.75]

Nruns=8

marker = ['o','s','^','d']
color = ['#d46a00','#00d46a', '#42fffc', '#6a00d4']

for mu in Chemical_Potentials:
    for patch_kind in Patch_Topologies:
        for pta in patch_A:
            for ptb in patch_B:
                for phi_b in phi_binary:
                    fig, ax = plt.subplots(3, 1, sharey=True, squeeze=False)
                    fig.set_canvas(plt.gcf().canvas)
                    
                    for i, energy in enumerate(Energies):   
                        df_equi = pd.DataFrame()
                        for nrun in range(1,Nruns+1):
                            dir_trace = 'mu_{}Energy_{}{}_patch_posA_{}posB_{}bin_{}_{}'.format(mu,energy,patch_kind,pta,ptb,phi_b,nrun)
                            fname = dir_trace+"/All_Clusters.dat"
                            f_exist = check_file(fname)
                            nexist=1
                            if f_exist == 1:              
                                dg = pd.read_csv(fname, header=None, delim_whitespace=True)
                                dg.columns=['time', 'cluster_size']
                                dg_equi = dg[dg.time > 100000 ]

                                if nexist == 2:
                                    df_equi = df.append(dg_equi)
                                    dg_max = dg_equi.groupby(['time']).aggregate('max').values
                                    df_max = df_max.append(dg_max)
                                    
                                if nexist == 1 and len(dg_equi)>0:
                                    df_equi = dg_equi
                                    df_max = df_equi.groupby(['time']).aggregate('max').values
                                    nexist+=1
                               

                        if len(df_equi)>0:
                            print(len(df_equi))
                            cluster_size = df_equi.cluster_size.values
                            # with weights
                            hist_p , bin_edges_p = np.histogram(cluster_size, 
                                bins=np.arange(1,np.max(cluster_size)+1), 
                                normed=True, 
                                weights=cluster_size)
                            # without weigts
                            hist, bin_edges = np.histogram(cluster_size, 
                                bins=np.arange(1,np.max(cluster_size)+1), 
                                normed=True)
                            # distribtion of maxima
                            hist_max, bin_edges_max = np.histogram(df_max, 
                                bins=np.arange(1,np.max(df_max)+1), 
                                normed=True)

                            ax[0,0].plot(bin_edges[:-1], hist, lw=2, marker=marker[i], color=color[i]) 
                            ax[0,0].set_xlim([1,18])
                            ax[0,0].set_ylim([-0.1,1.1])
                            ax[0,0].set_ylabel('P')
                            ax[0,0].set_xlabel('cluster size')
                           
                            ax[1,0].plot(bin_edges_p[:-1], hist_p, lw=2, marker=marker[i], color=color[i])
                            ax[1,0].set_xlim([1,18])
                            ax[1,0].set_ylim([-0.1,1.1])
                            ax[1,0].set_ylabel('P')
                            ax[1,0].set_xlabel('cluster size')
                                         
                            ax[2,0].plot(bin_edges_max[:-1], hist_max, lw=2, marker=marker[i], color=color[i])
                            ax[2,0].set_xlim([1,18])
                            ax[2,0].set_ylim([-0.1,1.1])
                            ax[2,0].set_ylabel('P')
                            ax[2,0].set_xlabel('maximum cluster size')
                    
                    plt.tight_layout()
                    plt.savefig("mu_{}{}_patch_posA_{}posB_{}bin{}.pdf".format(mu, patch_kind, pta,ptb,phi_b))
                    plt.close()
                
                       


"""
pd.read_csv("mu_0.3Energy_8.2symm_patch_posA_0.3posB_0.3bin_0.25_1/All_Clusters.dat", header=None, delim_whitespace=True)
df = pd.read_csv("mu_0.3Energy_8.2symm_patch_posA_0.3posB_0.3bin_0.25_1/All_Clusters.dat", header=None, delim_whitespace=True)
df
df = pd.read_csv("mu_0.3Energy_8.2symm_patch_posA_0.3posB_0.3bin_0.25_2/All_Clusters.dat", header=None, delim_whitespace=True)
df
df = pd.read_csv("mu_0.3Energy_8.2symm_patch_posA_0.3posB_0.3bin_0.25_3/All_Clusters.dat", header=None, delim_whitespace=True)
df
df.columns=['time', 'cluster_size']
df
df_equi = df[df.time > 100000 ]
df
df_equi.groupby('time')
df_equi.groupby(['time']).hist()
import matplotlib.pyplot()
import matplotlib.pyplot as plt
plt.show()
df_equi = df[df.time > 500000 ]
fig, ax = plt.subplots()
df_equi.groupby(['time']).hist()
plt.show()
df_equi = df[df.time > 1000000 ]
df_equi.groupby(['time']).hist()
plt.show()
df = pd.read_csv("mu_0.3Energy_10.2symm_patch_posA_0.3posB_0.3bin_0.25_3/All_Clusters.dat", header=None, delim_whitespace=True)
df = pd.read_csv("mu_0.3Energy_10.2symm_patch_posA_0.3posB_0.3bin_0.25_4/All_Clusters.dat", header=None, delim_whitespace=True)
df = pd.read_csv("mu_0.3Energy_8.2symm_patch_posA_0.3posB_0.3bin_0.25_4/All_Clusters.dat", header=None, delim_whitespace=True)
df_equi = df[df.time > 1000000 ]
df.columns=['time', 'cluster_size']
df_equi = df[df.time > 1000000 ]
df_equi.groupby(['time']).hist()
plt.show()
df_equi.groupby(['time']).aggregate('max')
df_max = df_equi.groupby(['time']).aggregate('max')
df_max
df_max.plot()
plot.show()
plt.show()
df_max.plot(kind='hist')
plt.show()
df_equi = df[df.time > 100000 ]
df_max = df_equi.groupby(['time']).aggregate('max')
df_max.plot(kind='hist')
plt.show()
df_max = df_equi.groupby(['time']).aggregate('max')
df_max
arr = df_max.values
arr
ar.shape
arr.shape
arr = df_max.values.flatten()
arr
df = pd.read_csv("mu_0.3Energy_5.2symm_patch_posA_0.3posB_0.3bin_0.25_4/All_Clusters.dat", header=None, delim_whitespace=True)
df_equi = df[df.time > 100000 ]
df.columns=['time', 'cluster_size']
df_equi = df[df.time > 100000 ]
df_max = df_equi.groupby(['time']).aggregate('max')
brr = df_max.values.flatten()
arr
brr
df_equi
df_equi.values
cluster_size = df_equi.cluster_size.values
counts, bins, patches = np.hist(cluster_size, bins=len(np.max(cluster_size)))
counts, bins, patches = np.hist(cluster_size, bins=len(np.max(cluster_size)), normed=True, wweights=np.arange(1,np.max(cluster_size)+1))
counts, bins, patches = np.histogram(cluster_size, bins=len(np.max(cluster_size)), normed=True, wweights=np.arange(1,np.max(cluster_size)+1))
counts, bins, patches = np.histogram(cluster_size, bins=np.max(cluster_size), normed=True, wweights=np.arange(1,np.max(cluster_size)+1))
counts, bins, patches = np.histogram(cluster_size, bins=np.max(cluster_size), normed=True, weights=np.arange(1,np.max(cluster_size)+1))
counts, bins, patches = np.histogram(cluster_size, bins=np.max(cluster_size), normed=True, weights=cluster_size)
hist , bin_edges = np.histogram(cluster_size, bins=np.max(cluster_size), normed=True, weights=cluster_size)
hist
np.max(cluster_size)
bin_edges
df = pd.read_csv("mu_0.3Energy_8.2symm_patch_posA_0.3posB_0.3bin_0.25_4/All_Clusters.dat", header=None, delim_whitespace=True)
df_equi = df[df.time > 100000 ]
df.columns=['time', 'cluster_size']
df_equi = df[df.time > 100000 ]
cluster_size = df_equi.cluster_size.values
hist , bin_edges = np.histogram(cluster_size, bins=np.max(cluster_size), normed=True, weights=cluster_size)
hist
cluster_size
bin_edges
hist , bin_edges = np.histogram(cluster_size, bins=np.arange(1,np.max(cluster_size)+1), normed=True, weights=cluster_size)
hist
bin_edges
df = pd.read_csv("mu_0.3Energy_8.2symm_patch_posA_0.5posB_0.5bin_0.5_4/All_Clusters.dat", header=None, delim_whitespace=True)
df = pd.read_csv("mu_0.3Energy_8.2symm_patch_posA_0.5posB_0.5bin_0.5_1/All_Clusters.dat", header=None, delim_whitespace=True)
df.columns=['time', 'cluster_size']
df_equi = df[df.time > 100000 ]
cluster_size = df_equi.cluster_size.values
hist , bin_edges = np.histogram(cluster_size, bins=np.arange(1,np.max(cluster_size)+1), normed=True, weights=cluster_size)
hist
hist , bin_edges = np.histogram(cluster_size, bins=np.arange(1,np.max(cluster_size)+1), normed=True)
hist
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('cd', '../')
get_ipython().run_line_magic('ls', '')
df = pd.read_csv("mu_0.3Energy_8.2symm_patch_posA_0.3posB_0.5bin_0.25_1/All_Clusters.dat", header=None, delim_whitespace=True)
df = pd.read_csv("mu_0.3Energy_8.2symm_patch_posA_0.3posB_0.5bin_0.25_1/All_Clusters.dat", header=None, delim_whitespace=True)
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('cd', 'lege_batterie_rl_chain/')
get_ipython().run_line_magic('ls', '')
df = pd.read_csv("mu_0.3Energy_8.2symm_patch_posA_0.3posB_0.5bin_0.25_1/All_Clusters.dat", header=None, delim_whitespace=True)
df.columns=['time', 'cluster_size']
df_equi = df[df.time > 100000 ]
cluster_size = df_equi.cluster_size.values
hist , bin_edges = np.histogram(cluster_size, bins=np.arange(1,np.max(cluster_size)+1), normed=True)
hist
hist , bin_edges = np.histogram(cluster_size, bins=np.arange(1,np.max(cluster_size)+1), normed=True, weights=cluster_size)
hist
get_ipython().run_line_magic('cd', '..')
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('cd', 'lege_batterie_rl_chain_Assym/')
get_ipython().run_line_magic('ls', '')
df = pd.read_csv("mu_0.3Energy_8.2Asymm_patch_posA_0.3posB_0.3bin_0.5_1/All_Clusters.dat", header=None, delim_whitespace=True)
df.columns=['time', 'cluster_size']
df_equi = df[df.time > 100000 ]
cluster_size = df_equi.cluster_size.values
hist , bin_edges = np.histogram(cluster_size, bins=np.arange(1,np.max(cluster_size)+1), normed=True, weights=cluster_size)
hist
np.max(cluster_size)
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('cd', '..')
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('cd', '..')
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('cd', '..')
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('cd', 'rr_chain/lege_batterie/')
get_ipython().run_line_magic('ls', '')
df = pd.read_csv("mu_0.3Energy_8.2Asymm_patch_pos_0.5_2/All_Clusters.dat", header=None, delim_whitespace=True)
df = pd.read_csv("mu_0.3Energy_8.2symm_patch_pos_0.5_2/All_Clusters.dat", header=None, delim_whitespace=True)
get_ipython().run_line_magic('ls', '')
df = pd.read_csv("mu_0.3Energy_8.2symm_patchpos_0.5_2/All_Clusters.dat", header=None, delim_whitespace=True)
df = pd.read_csv("mu_0.3Energy_8.2symm_patchpos_0.5_2/All_Clusters.dat", header=None, delim_whitespace=True)
get_ipython().run_line_magic('save', 'mysession 0-137')

"""
