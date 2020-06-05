import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
import matplotlib as mpl
import matplotlib.style as style
from matplotlib import patches


style.use('seaborn-ticks') 
mpl.rcParams['font.family'] = "sans-serif"
#sns.set_context('poster')
plt.rcParams['axes.axisbelow'] = True

plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['figure.titlesize'] = 15 

blue_c ='#9999FF'
red_c ='#FF9999'
purple_c='#C17DCB'
green_c='#61DCB7'


dd = 0.05
dr = 0.01

# paths for reading and saving 
savepath = "../plots/varr/"
filelist = glob.glob("../analysis_data/varr_mp_merged/*")

# read pickled files 
df = pd.read_pickle(filelist[0])
for file in filelist[1:]:
    df = pd.concat([df,pd.read_pickle(file)])

# convert types
convert_dict = {'#state': int,
                '#cp': int,
                'mu': float, 
                'delta':float,
                'r_patch':float, 
                'percolation':float, 
                #'#p', 
                #'#np', 
                'psi':float, 
                #'max_cluster', 
                #'max_p_domain', 
                #'max_np_domain',
                'polybonds/bonds':float, 
                'avgpolysize':float,
                'max_bonds':float,
                #'theta_diff_list',
                #'g(r)',
                #'cluster_sizes',
                #'p_domain_sizes',
                #'np_domain_sizes',
                #'energy',
                #'density,'
               } 
  
df = df.astype(convert_dict) 
df['max_bonds'] = np.where(df['max_bonds'] < 1, 1, df['max_bonds'])


# Function to plot psi 
def plotpsi(data,std,radius,deltas):
    plt.errorbar(deltas, data, yerr = std,
                 label = f'$r$ = {radius}', linewidth=1, marker="o", ms=7, capsize=3, capthick=2)
    plt.legend(loc='upper right', bbox_to_anchor=(1, 0.95), ncol=2)
    plt.ylim(-1,1)
    plt.xlabel("$\Delta$")
    plt.ylabel("$\Psi}$")

def get_edge_points(pos_i,ax_n,sign_p):
    edge_n = np.zeros(2)
    edge_n = pos_i + sign_p[0]*ax_n[0]/2. + sign_p[1]*ax_n[1]/2.

    return edge_n


def plot_gr(g_av,delta,radius,mu):
    x = np.linspace(0,len_g,601)/100
    fig, ax = plt.subplots(figsize=(13,3))


    left, bottom, width, height = [0.8, 0.4, 0.06, 0.3]
    ax2 = fig.add_axes([left, bottom, width, height])

    sin60 = np.sin(np.pi/3.)
    cos60 = np.cos(np.pi/3.)

    ax_n = np.array([[1,cos60],[0,sin60]])
    pos = np.array([0,0])
    edges=np.zeros((4,2))

    edges[0] = get_edge_points(pos,ax_n,np.array([-1,-1]))
    edges[1] = get_edge_points(pos,ax_n,np.array([+1,-1]))
    edges[2] = get_edge_points(pos,ax_n,np.array([+1,+1]))
    edges[3] = get_edge_points(pos,ax_n,np.array([-1,+1]))

    pdelta = 0.5
    particle_patches = edges[0] + pdelta*(edges[3]-edges[0])
    #particle_patches[1] = edges[2] + pdelta[1]*(edges[3]-edges[2])
    #particle_patches[2] = edges[0] + pdelta[2]*(edges[1]-edges[0])
    #particle_patches[3] = edges[1] + pdelta[3]*(edges[2]-edges[1])

    rhombi = patches.Polygon(edges,
                                linewidth=0.5,
                                edgecolor='k',
                                facecolor='white',
                                alpha=1)

    pcolor = 'r' 
    patch = patches.Circle((particle_patches[0],particle_patches[1]),
                                radius=radius,
                            facecolor=pcolor)
    ax2.add_patch(patch)
    ax2.add_patch(rhombi)
    ax2.set_xlim((-1,1))
    ax2.set_ylim((-1,1))
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.axes.xaxis.set_visible(False)
    ax2.axes.yaxis.set_visible(False)


    ax.plot(x,g_av, label=f'$\Delta={delta},r={radius}$',  c='#147448', lw=3 )
    ax.legend(loc='best')
    ax.set_xlim((0.8,6.0))
    ax.set_ylim((0,10))
    ax.set_xlabel("l")
    ax.set_ylabel("g(l)")
    #ax.plot(5.7,3.8, marker='o', ms=radius*150, c='r')
    #plt.tight_layout()
    plt.savefig("../plots/gr_{}_delta_{}_mu_{}.png".format(radius,delta,mu), dpi=300)
    plt.close()


mus = np.sort(np.array(df.mu.unique()))
deltas = np.sort(np.array(df.delta.unique()))
radii = np.sort(np.array(df.r_patch.unique()))
len_g = len(df['g'][0])

mu_dict ={x:i for i,x in enumerate(mus)}
del_dict={x:i for i,x in enumerate(deltas)}
rad_dict={x:i for i,x in enumerate(radii)}

# initialize plot arrays 
bond_max  = np.zeros((len(mus), len(radii), len(deltas)))
bond_avg  = np.zeros((len(mus), len(radii), len(deltas)))
polybonds = np.zeros((len(mus), len(radii), len(deltas)))
polyratio = np.zeros((len(mus), len(radii), len(deltas)))
psi_avg   = np.zeros((len(mus), len(radii), len(deltas)))
psi_std   = np.zeros((len(mus), len(radii), len(deltas)))
bonds     = np.zeros((len(mus), len(radii), len(deltas), 2)) #rename
g_avg     = np.zeros((len(mus), len(radii), len(deltas), len_g))#160??


for properties, ensemble in df.groupby(['mu', 'delta', 'r_patch']):
    mu = properties[0]
    delta = properties[1]
    radius = properties[2]
    i_m = mu_dict[mu]
    i_d = del_dict[delta]
    i_r = rad_dict[radius]

    bond_max[i_m,i_r,i_d] = np.amax(ensemble['max_bonds'].values)
    bond_avg[i_m,i_r,i_d] = np.mean(ensemble['max_bonds'].values)

    psi_avg[i_m,i_r,i_d] = np.mean(ensemble['psi'].values)
    psi_std[i_m,i_r,i_d] = np.std(ensemble['psi'].values)

    g_avg[i_m,i_r,i_d] = np.average(ensemble['g'].values, axis=0)

    bonds[i_m,i_r,i_d,0] = ensemble['psi'][ensemble['psi'] > 0.3].count()
    bonds[i_m,i_r,i_d,1] = ensemble['psi'][ensemble['psi'] < -0.3].count()


#plot Psi
sns.set_palette(sns.cubehelix_palette(8, start=.5, rot=-.75))
r_to_plot = [0.05,0.08,0.10,0.15,0.2]
d_to_plot = [0.2,0.3,0.4,0.5]
for i_m in range(len(mus)):
    fig,ax=plt.subplots(figsize=(5,7))
    filename="psi_op_var_r_{}.pdf".format(mus[i_m])
    for radius in r_to_plot:
            data = psi_avg[i_m,rad_dict[radius]]
            std = psi_std[i_m,rad_dict[radius]]
            plotpsi(data,std,radius,deltas)
    plt.tight_layout()
    plt.savefig(savepath+filename)
    for radius in radii:  
        for delta in d_to_plot:
            data = g_avg[i_m,rad_dict[radius],del_dict[delta]]
            plot_gr(data,delta, radius, mus[i_m])
   

