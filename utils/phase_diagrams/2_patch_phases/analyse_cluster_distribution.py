import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import os.path


def check_file(fname):
    try:
        open(fname,"r")
        return 1
    except IOError:
        print "Error: File doesn't seem to exist."
        return 0


Chemical_Potentials = [0.3,0.4]
Energies = [ 5.2,6.2,7.2,8.2,9.2,10.2]
Patch_Toplogies = ['symm']
Patch_Positions = [0.2, 0.3 ,0.4 ,0.5, 0.6, 0.7, 0.8] 

Nruns=8
N_vals = np.zeros(4)

NVAL = np.zeros((len(Energies),len(N_vals)))

N_vals_single = np.zeros((len(N_vals),Nruns))
N_VALS_std = np.zeros((len(Energies), len(N_vals)))

Monomers = np.zeros(len(Energies))
Dimers = np.zeros(len(Energies))
Trimers = np.zeros(len(Energies))
Multimers = np.zeros(len(Energies))

phi_mean=np.zeros(len(Energies))
phi_std = np.zeros(len(Energies))

"""
Monomers = [[] for k in xrange(len(Energies)) ]
Dimers = [[] for k in xrange(len(Energies)) ]
Trimers = [[] for k in xrange(len(Energies)) ]
Quadmers = [[] for k in xrange(len(Energies)) ]
Mulimers = [[] for k in xrange(len(Energies)) ]
"""

f = open("cluster_distribution.csv", 'w')
f.write("chemical_potential,energy,phi_mean, phi_std, patch_topoloy,patch_position,N1,N2,N3,N4,btN4\n")

for patch_kind in Patch_Toplogies:
    for patch_pos in Patch_Positions:
        for mu in Chemical_Potentials:
            epsi=0
            for energy in Energies:
                data=np.zeros(1)
                Phi=np.zeros(1)
                for nrun in xrange(1,Nruns+1):
                    print nrun
                    directory = "mu_"+str(mu)+"Energy_"+str(energy)+patch_kind+"_patchpos_"+str(patch_pos)+"_"+str(nrun)
                    fname = directory+"/All_Clusters.dat"
                    f_exist = check_file(fname)

                    if f_exist == 1: 
                        df = pd.read_csv(fname, delim_whitespace=True, header=None)                                                                                                                                                    
                        df.columns=['time','size']
                        
                        data_single = df['size'].values
                        data = np.concatenate((data, data_single))
                        
                        d = np.diff(np.unique(data_single)).min()
                        left_of_first_bin = data_single.min() - float(d)/2
                        right_of_last_bin = data_single.max() + float(d)/2

                        N_single, bins, patches = plt.hist(data_single, 
                        np.arange(left_of_first_bin, right_of_last_bin + d, d),
                        weights=data_single,
                        alpha=0.5,
                        normed = True)

                        for nl in xrange(min(3,len(N_single))):
                            N_vals_single[nl,nrun-1] = N_single[nl]

                        N_vals_single[3,nrun-1] = np.sum(N_single[3:])

                        fname2=directory+"/NPT_OUT.txt"
                        dg = pd.read_csv(fname2, delim_whitespace=True ,header=None)
                        dg.columns=['time', 'crap', 'N', 'crap2', 'phi']
                        Phi=np.concatenate((Phi, dg['phi'].values))
                    
                
                N_VALS_std[epsi,0] = np.sqrt((np.var(N_vals_single[0])+np.var(N_vals_single[1]))/2.)
                N_VALS_std[epsi,2] = np.std(N_vals_single[2])
                N_VALS_std[epsi,3] = np.std(N_vals_single[3])

                data = np.delete(data,0)
                Phi = np.delete(Phi,0)
                print epsi,Phi
                phi_mean[epsi] = np.mean(Phi)
                print phi_mean[epsi]
                phi_std[epsi] = np.std(Phi)

                d = np.diff(np.unique(data)).min()
                left_of_first_bin = data.min() - float(d)/2
                right_of_last_bin = data.max() + float(d)/2
                N_cluster, bins, patches = plt.hist(data, 
                    np.arange(left_of_first_bin, right_of_last_bin + d, d),
                    weights=data,
                    alpha=0.5,
                    normed = True)
                for nl in xrange(min(3,len(N_cluster))):
                    N_vals[nl] = N_cluster[nl]

                N_vals[3] = np.sum(N_cluster[3:])
                NVAL[epsi] = N_vals
                f.write(str(mu)+
                    ","+str(energy)+
                    ","+str(phi_mean[epsi])+
                    ","+str(phi_std[epsi])+
                    ","+patch_kind+
                    ","+str(patch_pos)+
                    ","+str(N_vals[0])+
                    ","+str(N_vals[1])+
                    ","+str(N_vals[2])+
                    ","+str(N_vals[3])+
                    "\n")
                epsi = epsi+1

            phi_mean = np.around(phi_mean, decimals=3)
            phi_std  = np.around(phi_std, decimals=3)

            N_epsilon = len(Energies)
            
            for epsi  in  xrange(N_epsilon):
                Monomers[epsi] = NVAL[epsi][0]
                Dimers[epsi]   = NVAL[epsi][1]
                Trimers[epsi]  = NVAL[epsi][2]
                Multimers[epsi]= NVAL[epsi][3]   

            print "plot data"
            error_config = dict(ecolor='black', lw=2, capsize=6, capthick=2, alpha=0.8)

            fig, ax =  plt.subplots(figsize=(15,11))
            width = 1
            ind = np.arange(N_epsilon) + 5.2 - width/2.
            lw=3
            plt.xticks(ind+width/2.)
            plt.xlim((5.2-width/2.,10.2+width/2.))
            plt.ylim((0,1))
            plt.tick_params(axis='both', which='major', labelsize=28)
            plt.xlabel("$\\epsilon [k_{B}T]$", size=38)
            plt.ylabel("$\\phi [\\%]$", size=38)
            plt.bar(ind, Monomers+Dimers, width, lw=lw, yerr=N_VALS_std[:,0], error_kw=error_config, color="#236AB9", alpha=0.7)
            #plt.bar(ind, Dimers, width, yerr=N_VALS_std[:,1], error_kw=error_config, color="#216869",bottom=Monomers, alpha=0.7)
            plt.bar(ind, Trimers, width, lw=lw, yerr=N_VALS_std[:,2], error_kw=error_config, color="#ff0066",bottom=Monomers+Dimers, alpha=0.9)
            plt.bar(ind, Multimers, width, lw=lw, yerr=N_VALS_std[:,3], error_kw=error_config, color="#ffff00", bottom=Monomers+Dimers+Trimers, alpha=0.7)
            
            for i in xrange(len(ind)):
                ax.text(ind[i], 1.07, "$\\psi=$"+str(phi_mean[i]), fontsize=24)
                ax.text(ind[i], 1.03, "$\pm$"+str(phi_std[i]), fontsize=24)


            plt.savefig("mu_"+str(mu)+patch_kind+"_patchpos_"+str(patch_pos)+".pdf")
            






