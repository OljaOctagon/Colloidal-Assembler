import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib as mpl
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib import rc 
import networkx as nx
sns.set(style='white')
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

colormap = plt.get_cmap('hsv')
#norm = mpl.colors.Normalize(0,np.pi/3.)

def get_pdist(df):
    Lx = df.Lx.unique()[0]
    Ly = df.Ly.unique()[0]
    
    df_x = df.x.values
    df_y = df.y.values
    l = len(df_x)
    print(df_x)
    pdist = np.zeros((l,l,2))

    p = np.reshape(df_x, (l,1))
    p = p - p.transpose()
    pdist[:,:,0] = p - Lx*np.rint(p/Lx)

    p = np.reshape(df_y, (l,1))
    p = p - p.transpose()
    pdist[:,:,1] = p - Ly*np.rint(p/Ly)

    N_pdist = np.sqrt( np.power( pdist[:,:,0], 2) + np.power( pdist[:,:,1], 2))
    return pdist, N_pdist


def hexatic_i(i,theta):
    return np.mean(np.exp(6*1j*theta))


def calc_max_domain_size(Angle,nn_list,TH):
    Np = len(Angle)
    edges = []
    for id in range(Np):
        for j in nn_list[id]:
            dangle = np.fabs(Angle[id] - Angle[j])
            if  dangle < TH:
                edges.append((id,j))


    G = nx.Graph()
    G.add_edges_from(edges)
    domains = list(nx.connected_components(G))
    domain_sizes = [len(domain) for domain in domains ]
    #print(domain_sizes,Np)
    max_domain = np.max(domain_sizes)/Np
    return max_domain

#Pressures = [5,6,8,10,11,12, 14, 16,17, 18, 19, 20, 25, 30, 35, 40, 50, 60, 80,100]
Pressures = [5,10,14,18,20, 25, 30, 35, 40, 50, 60, 80,100]
patch_positions = [0.3]
Nruns=8
phi = [ [] for i in range(Nruns) ]

filetrace = 'mu_0.3Energy_8.2Asymm_patchpos_'
fig,ax = plt.subplots()
eqs = []
domains=[]
for p in Pressures:
    print(p)
    max_domain=[]
    for pos in patch_positions:
        for nrun in range(1, Nruns+1):
            run_trace = filetrace+str(pos)+"_Pressure_"+str(p)+"_"+str(nrun)
            file = "/cluster_center_of_mass.dat"
            dfile = run_trace+file
            df = pd.read_csv(dfile, header=None, delim_whitespace=True, 
                names=['time','N','Lx', 'Ly', 'id', 'x', 'y'])

            # get positions 
            Times = df.time.values[-1:]
            for frame in Times:
                nn_list = []
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
                    nn_list.append(nn_id)
                    theta = np.arccos(pdist[id, nn_id,0]/N_pdist[id,nn_id])
                    psi = hexatic_i(id, theta)
                    angle = np.arctan2(psi.real,psi.imag) + np.pi
                    Angle[id] = angle 

                TH=15*(2*np.pi/360)
                size_l=calc_max_domain_size(Angle, nn_list,TH)
                max_domain.append(size_l)
                with open('all_domains.dat', 'a') as f:
                    f.write("{} {} {} {}\n ".format(pos, p, nrun, size_l))

                
                fig,ax = plt.subplots()
                plt.scatter(df_sub.x.values,
                            df_sub.y.values,
                            c=Angle, cmap=colormap,norm = mpl.colors.Normalize(vmin=0.,vmax=2*np.pi),
                            label="\displaystyle P = "+str(p))
                plt.axis('equal')
                plt.axis('off')
                plt.tight_layout()
                plt.savefig("hexatic_pressure_{}_patch_{}_{}.png".format(p,pos,nrun), dpi=300)
                plt.close()
                

    domains.append(np.array([np.mean(max_domain), np.std(max_domain)]))

domains = np.array(domains)

plt.errorbar(Pressures, domains[:,0], domains[:,1],
             capthick=3,
             capsize=5,
             ms=5, marker='o', c='k', alpha=0.5, markeredgecolor='k')

ax.tick_params(axis='both', which='major', labelsize=22)
plt.xlabel("$P$", size=30)
plt.ylabel("$ \langle S_{L}\\rangle $", size=30)
plt.tight_layout()
plt.savefig("max_domain_threshold_15.png", dpi=300)

Pressures = np.reshape(np.array(Pressures), (-1,1))
out_arr_domains = np.concatenate((Pressures,domains), axis=1)
np.savetxt("domains.dat", out_arr_domains, newline='\n')

