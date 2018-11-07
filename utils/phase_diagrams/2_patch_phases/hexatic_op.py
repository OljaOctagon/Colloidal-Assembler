import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib as mpl
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib import rc 

sns.set(style='white')
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

colormap = plt.get_cmap('hsv')
#norm = mpl.colors.Normalize(0,np.pi/3.)

def get_pdist(df):
    Lx = df.Lx.unique()[0]
    Ly = df.Ly.unique()[0]
    l = len(df.x)
    pdist = np.zeros((l,l,2))

    p = np.reshape(df.x, (l,1))
    p = p - p.transpose()
    pdist[:,:,0] = p - Lx*np.rint(p/Lx)

    p = np.reshape(df.y, (l,1))
    p = p - p.transpose()
    pdist[:,:,1] = p - Ly*np.rint(p/Ly)

    N_pdist = np.sqrt( np.power( pdist[:,:,0], 2) + np.power( pdist[:,:,1], 2))
    return pdist, N_pdist


def hexatic_i(i,theta):
    return np.mean(np.exp(6*1j*theta))


#Pressures = [5,6,8,10,11,12, 14, 16,17, 18, 19, 20, 25, 30, 35, 40, 50, 60, 80,100]
#patch_positions = [0.3, 0.5, 0.7]
#Nruns = 8 

Nruns=8
patch_positions = [0.7]
Pressures = [35]
phi = [ [] for i in range(Nruns) ]

filetrace = 'mu_0.3Energy_8.2symm_patchpos_'
for pos in patch_positions:
    eqs = []
    for p in Pressures:
        for nrun in range(1, Nruns+1):
            run_trace = filetrace+str(pos)+"_Pressure_"+str(p)+"_"+str(nrun)
            file = "/cluster_center_of_mass.dat"
            dfile = run_trace+file
            df = pd.read_csv(dfile, header=None, delim_whitespace=True, 
                names=['time','N','Lx', 'Ly', 'id', 'x', 'y'])

            # get positions 
            Times = df.time.values[-1:]
            for frame in Times: 
                df_sub = df[df.time == frame]
                pdist, N_pdist = get_pdist(df_sub)
                # take first 12 neighbors
                sorted_dist = np.sort(N_pdist)[:,1:12]
                N=int(df_sub.N.unique()[0])
                Angle = np.zeros(N)
                for id in range(N):
                    counts, bin_edges = np.histogram(sorted_dist[id,:])
                    nneigh = np.argmax(counts)+np.argmin(counts[np.argmax(counts):])
                    r_cutoff = bin_edges[nneigh+1]

                    nn_id = np.argwhere(((N_pdist[id,:]<r_cutoff) & (N_pdist[id,:]>0))).flatten()    
                    print(id, nn_id) 
                    theta = np.arccos(pdist[id, nn_id,0]/N_pdist[id,nn_id])
                    psi = hexatic_i(id, theta)
                    angle = np.arctan2(psi.real,psi.imag) + np.pi
                    print(angle)
                    if angle > np.pi/3.:
                        angle = angle - (np.pi/3.)

                    Angle[id] = angle

                fig,ax = plt.subplots()
                plt.scatter(df_sub.x.values, df_sub.y.values, c=Angle, cmap=colormap, label="\displaystyle P = "+str(p))
                plt.axis('equal')
                #plt.xlim([0,35])
                #plt.ylim([0,30])
                plt.tick_params(axis='both', which='major', labelsize=14)

                #plt.legend(prop={'size':15}, fancybox = True)
                plt.xlabel("$\displaystyle x $", size=26)
                plt.ylabel("$\displaystyle y $", size=26)
                plt.tight_layout()
                plt.savefig("hexatic_pressure_35_patch_0.7_"+str(nrun)+".pdf")


 