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

##################


def plot_domain_sizes(arr_p, arr_np,radius,mu,delta):

    sns.set(style='white')

    fig, ax = plt.subplots()
    #nbins_p = int(np.ceil(np.max(arr_p)/2.))
    #nbins_np = int(np.ceil(np.max(arr_np)/2.))
    nbins = 100
    

    ax.hist(arr_p, alpha=0.4,bins=nbins, normed=True, facecolor='red', label='parallel')
    ax.hist(arr_np,alpha=0.4, bins=nbins, normed=True, facecolor='blue', label='non-parallel')
    ax.set_xlabel("domain sizes")
    ax.set_ylabel("P")
    #handles, labels = ax.get_legend_handles_labels()

    fig.legend(loc=(0.35,0), ncol=2)
    plt.tight_layout()
    plt.savefig("../plots/domain_sizes/domain_size_histogram_{}_{}_{}.png".format(mu,delta,radius), dpi=300)


#################
def get_edge_points(pos_i,ax_n,sign_p):
    edge_n = np.zeros(2)
    edge_n = pos_i + sign_p[0]*ax_n[0]/2. + sign_p[1]*ax_n[1]/2.

    return edge_n



mus = np.sort(np.array(df.mu.unique()))
deltas = np.sort(np.array(df.delta.unique()))
radii = np.sort(np.array(df.r_patch.unique()))

# initialize plot arrays 


p_distri = {}
np_distri = {}
theta = {}
packing = {} 

for properties, ensemble in df.groupby(['mu', 'delta', 'r_patch']):
    mu = properties[0]
    delta = properties[1]
    radius = properties[2]

    p_arr = ensemble['p_domain_sizes'].values
    np_arr = ensemble['np_domain_sizes'].values
    theta_arr =  ensemble['theta_diff_list'].values
    packing_arr = ensemble['density'].values

    p_distri[properties] = []
    np_distri[properties] = []
    theta[properties] = []
    packing[properties] = packing_arr



    for item in p_arr:
        p_distri[properties] = np.concatenate((p_distri[properties],item))
    for item in np_arr:
        np_distri[properties] = np.concatenate((np_distri[properties],item))

    for item in theta_arr:
        theta[properties] = np.concatenate((theta[properties], item))


    #plot_domain_sizes(p_distri,np_distri,radius,mu,delta)
    #plot_av_domain(p_distri,np_distri,radius,mu,delta)

radii = [0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2]
for mu in mus:
    for delta in [0.2,0.5,0.35]:
        bins = np.array([0,5,15,25,35,45,55,65,75,85,95,105,120,150,200,1000])
        np_distri_map = np.zeros((len(radii), (len(bins)-1)))
        fig,ax = plt.subplots(figsize=(10,10))
        for ri, radius in enumerate(radii):
            arr = np_distri[(mu,delta,radius)]
            arr_sum = np.sum(arr)
            for bi in range(len(bins)-1):
                np_distri_map[ri,bi] = np.sum(arr[(arr>=bins[bi]) & (arr<bins[bi+1])])/arr_sum

        cmap = sns.cubehelix_palette(10, start=.5, rot=-.75, as_cmap=True, reverse=True)
        #plt.imshow(np_distri_map, cmap='bone', interpolation='none')
        sns.heatmap(np_distri_map, cmap=cmap, vmin=0,vmax=0.8, linewidths=0.5, linecolor='k')
        plt.xlabel("domain sizes")
        plt.ylabel("radius")
        ax.set_xticks(np.arange(len(bins[:-1])))
        ax.set_yticks(np.arange(len(radii)))
        ax.set_xticklabels(bins[:-1])
        ax.set_yticklabels(radii)
        plt.yticks(rotation=0) 
        #plt.savefig("../plots/np_domain_sizes_{}_{}.pdf".format(mu,delta))

        ##############3
        bins = np.array([0,5,15,25,35,45,55,65,75,85,95,105,120,150,200,1000])
        p_distri_map = np.zeros((len(radii), (len(bins)-1)))
        fig,ax = plt.subplots(figsize=(10,10))
        for ri, radius in enumerate(radii):
            arr = p_distri[(mu,delta,radius)]
            arr_sum = np.sum(arr)
            for bi in range(len(bins)-1):
                p_distri_map[ri,bi] = np.sum(arr[(arr>=bins[bi]) & (arr<bins[bi+1])])/arr_sum

        cmap = sns.cubehelix_palette(10, start=0.2, light=1, as_cmap=True,reverse=True)
        sns.heatmap(p_distri_map, cmap=cmap,vmin=0, vmax=0.8, linewidths=0.5, linecolor='k')
        #plt.imshow(p_distri_map, cmap='pink', interpolation='none')
        plt.xlabel("domain sizes")
        plt.ylabel("radius")
        ax.set_xticks(np.arange(len(bins[:-1])))
        ax.set_yticks(np.arange(len(radii)))
        ax.set_xticklabels(bins[:-1])
        ax.set_yticklabels(radii)
        plt.yticks(rotation=0) 
        #plt.savefig("../plots/p_domain_sizes_{}_{}.pdf".format(mu,delta))

        sns.set_palette(sns.cubehelix_palette(16, start=0.5, rot=-0.75))

        for RI, RADIUS in enumerate(radii):

            fig,ax = plt.subplots(figsize=(10,10))
            for ri, radius in enumerate(radii):
                nbins = 200
                arr = theta[(mu,delta,radius)]
                if ri != RI:
                    n,bins,npatches = ax.hist(arr, histtype='step', lw=2, alpha=0.4,
                                            bins=nbins, density=True, label='r = {}'.format(radius))
            ax.set_xlabel("mutual orientations")
            ax.set_ylabel("P")
            ax.set_xticks((0,np.pi/3,2*np.pi/3,np.pi))
            ax.set_xticklabels(['0','$\pi$/3', '$2\pi/3$', '$\pi$'])
            ax.set_xlim((0,np.pi))
            ax.set_ylim((0,6))

            nbins = 200
            arr = theta[(mu,delta,RADIUS)]
            n,bins,npatches = ax.hist(arr, histtype='step', lw=3, alpha=1, color='red',
                                        bins=nbins, density=True, label='r = {}'.format(RADIUS))


            plt.legend()
            plt.tight_layout()

            left, bottom, width, height = [0.6, 0.75, 0.2, 0.2]
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

            pdelta = delta 
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
                                        radius=RADIUS,
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

            #plt.savefig("../plots/orientation/orientation_histogram_{}_{}_rad_{}.pdf".format(mu,delta,RADIUS))
            plt.close()

print("INFO len delta", len(deltas))
for mu in mus:
    sns.set_palette(sns.cubehelix_palette(10, start=.5, rot=-.75, reverse=True))
    fig,ax = plt.subplots(figsize=(10,7))
    for delta in deltas:
        arr_mean = np.zeros(len(radii))
        arr_std = np.zeros(len(radii))
        for i,radius in enumerate(radii):
            arr_mean[i] = np.mean(packing[(mu,delta,radius)])
            arr_std[i] = np.std(packing[(mu,delta,radius)])
        plt.errorbar(radii, arr_mean, yerr = arr_std,
                     label = f'$\Delta$ = {delta}', linewidth=2, mec='k', marker="o", ms=5, capsize=3, capthick=2)
    plt.legend(loc='upper right', bbox_to_anchor=(1, 0.3), ncol=2)
    plt.ylim((0.4,1))
    plt.xlabel("$r$")
    plt.ylabel("$\phi}$")
    plt.tight_layout()
    plt.savefig("../plots/packing/packing_mu_{}.pdf".format(mu))

