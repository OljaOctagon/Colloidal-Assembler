import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter

sns.set(style="white")
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def gr_per_run(data, Delta, R, delta_time):
    a=Delta
    G = np.zeros(int((R-a)/Delta))
    start_time=data['time'].values[0]
    end_time=data['time'].values[-1]
    norm_time=(end_time-start_time)/float(delta_time)

    for t in range(start_time, end_time, delta_time):
        print(t)
        is_time = data["time"] == t
        P_t = data[is_time]
        is_id0 = P_t["id"] == 0
        Np = P_t[is_id0]["N"].values
        lx = P_t[is_id0]["Lx"].values
        ly = P_t[is_id0]["Ly"].values

        center_x=P_t["X"]
        center_y=P_t["Y"]

        particle_distx = np.reshape(center_x, (len(center_x), 1))
        particle_distx = particle_distx - particle_distx.transpose()

        particle_disty = np.reshape(center_y, (len(center_y), 1))
        particle_disty = particle_disty - particle_disty.transpose()

        particle_distx = particle_distx - lx*np.rint(particle_distx/lx)
        particle_disty = particle_disty - ly*np.rint(particle_disty/ly)

        particle_dist=np.sqrt(np.power(particle_distx,2) + np.power(particle_disty,2))

        rho=float(Np)/(lx*ly)
        S=int((R-a)/Delta)

        for i in range(1,S):
            rb = a+i*Delta
            ra = rb-Delta        
            A = np.pi*(rb*rb - ra*ra)      
            normed = rho*A
            gr = 0
            for id in range(int(Np)):
                dx= particle_dist[id]
                na=len(dx[np.where((dx > ra) & (dx < rb))])      
                gr = gr + na/normed
                                  
            G[i] = G[i] + gr/float(Np)

    G=G/float(norm_time)
    return G

def check_file(fname):
    try:
        open(fname,"r")
        return 1
    except IOError:
        print("Error: File doesn't seem to exist.")
        return 0

Chemical_Potentials = [0.3]
Energies = [8.2]
Patch_Toplogies = ['symm']
Patch_Positions = [0.7] 
Pressure=[100]
Nruns=8
Delta=0.01
R=10
delta_time=100

for patch_kind in Patch_Toplogies:
    for patch_pos in Patch_Positions:
        for mu in Chemical_Potentials:
            for energy in Energies:
                for p in Pressure:
                    G = np.zeros(int((R-Delta)/Delta))
                    for nrun in range(1,Nruns+1):
                        print(nrun)
                        directory = "mu_"+str(mu)+"Energy_"+str(energy)+patch_kind+"_patchpos_"+str(patch_pos)+"_Pressure_"+str(p)+"_"+str(nrun)
                        print(directory)
                        fname = directory+"/cluster_center_of_mass.dat"
                        f_exist = check_file(fname)

                        if f_exist == 1: 
                            data = pd.read_csv(fname, delim_whitespace=True, header=None)                                                                                                                                                    
                            data.columns=["time","N", "Lx", "Ly", "id", "X", "Y"]
                            G = G + gr_per_run(data, Delta, R, delta_time)

                    G = G/Nruns  

                    np.savetxt("gr_cluster_mu_"+str(mu)+patch_kind+"_patchpos_"+str(patch_pos)+".dat",G,newline="\n", delimiter=" ")            
                    r=np.arange(Delta,R,Delta)/2.    
                    plt.plot(r,np.ones(len(r)), 'k--', lw=2, alpha=0.2)
                    plt.plot(r,G,c="#3399ff", lw=2, label="box phase at $ \displaystyle \phi \\approx 0.78$")
                    plt.xlabel("$\displaystyle \sigma$", size=30)
                    plt.tick_params(labelsize=22)
                    plt.ylabel("$\displaystyle g(\sigma)$", size=30)
                    plt.locator_params(nbins=6)
                    plt.legend(prop={'size':18})
                    plt.xlim([0,5])
                    plt.ylim([0,10])

                    plt.tight_layout()
                    plt.savefig("gr_cluster_mu_"+str(mu)+patch_kind+"_patchpos_"+str(patch_pos)+".pdf")
                    plt.show()      














                        