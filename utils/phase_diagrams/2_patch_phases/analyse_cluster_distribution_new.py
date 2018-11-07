import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import os.path

def check_file(fname):
    try:
        open(fname,"r")
        return 1
    except IOError:
        print("Error: File doesn't seem to exist.")
        return 0

Chemical_Potentials = [0.3,0.4]
Energies = [ 5.2,6.2,7.2,8.2,9.2,10.2]
Patch_Toplogies = ['Asymm']
Patch_Positions = [0.2, 0.3 ,0.4 ,0.5, 0.6, 0.7, 0.8] 

Nruns=8
# stars
bins = [0,5,6,7,1000]

# boxes 
#bins = [0,3,4,1000]
n_bars = 4
cluster_sizes = np.zeros(n_bars)
    
f = open("cluster_distribution.csv", 'w')
f.write("chemical_potential,energy,phi_mean,phi_std,patch_topoloy,patch_position,Nleq5,N5,N6,Ngeq6\n")

for patch_kind in Patch_Toplogies:
    for mu in Chemical_Potentials:
        for patch_pos in Patch_Positions: 
            cluster_sizes = np.zeros((len(Energies), n_bars))
            std_cluster_sizes = np.zeros((len(Energies), n_bars))
            N_times = np.zeros(len(Energies))
            phi = np.zeros(len(Energies))
            phi_std = np.zeros(len(Energies))
            for ei, energy in enumerate(Energies):
                for nrun in range(1,Nruns+1):
                    directory = "mu_"+str(mu)+"Energy_"+str(energy)+patch_kind+"_patchpos_"+str(patch_pos)+"_"+str(nrun)
                    fname = directory+"/All_Clusters.dat"
                    f_exist = check_file(fname)

                    if f_exist == 1: 
                        # Calculate packing fraction 
                        if nrun == 1:
                            dg = pd.read_csv(directory+"/NPT_OUT.txt", header=None, delim_whitespace=True)
                        else:
                            de = pd.read_csv(directory+"/NPT_OUT.txt", header=None, delim_whitespace=True)
                            dg = dg.append(de)

                        # Read cluster data 
                        df = pd.read_csv(fname, delim_whitespace=True, header=None)                                                                                                                                                    
                        df.columns=['time','size']
                        unique_times = np.unique(df.time.values)
                        for time in unique_times:
                            cluster_size = df[df.time==time]['size'].values  
                        
                            # with percentage of particle in cluster 
                            cluster_size = np.repeat(cluster_size, cluster_size)
                            cluster_size_percent, bin_edges = np.histogram(cluster_size, 
                                bins=bins)

                            cluster_size_percent = cluster_size_percent/len(cluster_size)
                            std_cluster = cluster_size_percent * cluster_size_percent

                            cluster_sizes[ei] = cluster_sizes[ei] + cluster_size_percent 
                            std_cluster_sizes[ei] = std_cluster_sizes[ei] + std_cluster


                        N_times[ei] += len(unique_times)

                cluster_sizes[ei] = cluster_sizes[ei]/N_times[ei] 
                std_cluster_sizes[ei] = np.sqrt( - np.power(cluster_sizes[ei],2) + (std_cluster_sizes[ei]/N_times[ei]))  

                arr = dg.values[:,4]
                phi[ei] = np.round(np.mean(arr), decimals=3)
                phi_std[ei] = np.round(np.std(arr), decimals=3)

                f.write(str(mu)+
                    ","+str(energy)+
                    ","+str(phi[ei])+
                    ","+str(phi_std[ei])+
                    ","+patch_kind+
                    ","+str(patch_pos)+
                    ","+str(cluster_sizes[ei,0])+
                    ","+str(cluster_sizes[ei,1])+
                    ","+str(cluster_sizes[ei,2])+
                    ","+str(cluster_sizes[ei,3])+
                    "\n") 
     
            colors = ["#236AB9", "#9900cc", "#ff0066", "#ffff00"]

            error_config = dict(ecolor='black', lw=2, capsize=6, capthick=2, alpha=0.8)
            fig, ax =  plt.subplots(figsize=(15,11))
            width = 1
            ind = np.arange(len(Energies)) + 5.2 - width/2.
            lw=3
            plt.xticks(ind+width/2.)
            plt.xlim((5.2-width/2.,10.2+width/2.))
            plt.ylim((0,1))
            plt.tick_params(axis='both', which='major', labelsize=28)
            plt.xlabel("$\\epsilon [k_{B}T]$", size=38)
            plt.ylabel("$\\phi [\\%]$", size=38)
                
            plt.bar(ind+0.5, cluster_sizes[:,0], width, lw=lw,yerr=std_cluster_sizes[:,0], error_kw=error_config, color=colors[0], alpha=0.7)
            plt.bar(ind+0.5, cluster_sizes[:,1], width, lw=lw,yerr=std_cluster_sizes[:,1], error_kw=error_config, color=colors[1],bottom=cluster_sizes[:,0],alpha=0.7)
            plt.bar(ind+0.5, cluster_sizes[:,2], width, lw=lw,yerr=std_cluster_sizes[:,2], error_kw=error_config, color=colors[2],bottom=cluster_sizes[:,0]+cluster_sizes[:,1], alpha=0.9)
            plt.bar(ind+0.5, cluster_sizes[:,3], width, lw=lw,yerr=std_cluster_sizes[:,3], error_kw=error_config, color=colors[3], bottom=cluster_sizes[:,0]+cluster_sizes[:,1]+cluster_sizes[:,2], alpha=0.7)
            
            for ei in range(len(Energies)):

                ax.text(ind[ei], 1.07, "$\\psi=$"+str(phi[ei]), fontsize=24)
                ax.text(ind[ei], 1.03, "$\pm$"+str(phi_std[ei]), fontsize=24)

            plt.savefig("mu_"+str(mu)+patch_kind+"_patchpos_"+str(patch_pos)+"_percent.pdf")
            

    












                        
                        