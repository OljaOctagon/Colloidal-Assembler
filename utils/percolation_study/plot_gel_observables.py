import argparse 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.style as style
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


def plot_domain_distribution(arr, dir_name):
    fig,ax = plt.subplots()
    plt.hist(arr,bins=range(1,np.max(arr)), fc=blue_c, edgecolor='gray', lw=1, alpha=0.7, density=True)
    plt.xlabel("cluster size")
    plt.ylabel("P")
    plt.locator_params(axis='y', nbins=6)
    plt.tight_layout()
    plt.savefig("{}cluster_size_distribution.pdf".format(dir_name))

def plot_degree_distribution(arr, dir_name):
    fig,ax = plt.subplots()
    plt.hist(arr,bins=[-0.5,0.5,1.5,2.5,3.5,4.5], fc=purple_c, edgecolor='gray', lw=1, alpha=0.7, density=True, rwidth=0.5)
    plt.xlabel("degree")
    plt.ylabel("P")
    plt.xticks([0,1,2,3,4])
    plt.xlim((-0.5,4.5))
    plt.locator_params(axis='y', nbins=6)
    plt.tight_layout()
    plt.savefig("{}degree_distribution.pdf".format(dir_name))


def plot_bond_probability(arr, dir_name):
    fig,ax = plt.subplots()

    plt.xlabel(" time [MC sweeps]")
    plt.ylabel("$p_{b}$")
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.plot(arr[:,0],arr[:,1],c=purple_c, lw=3,ms=3,marker='o')
    plt.tight_layout()
    plt.savefig("{}pbt.pdf".format(dir_name))


def plot_connectivity(arr, dir_name):
    fig,ax = plt.subplots()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    plt.xlabel("time [MC sweeps]")
    plt.ylabel("node connectivity")
    plt.plot(arr[:,0],arr[:,1],c=blue_c, lw=3,ms=3,marker='o')
    plt.tight_layout()
    plt.savefig("{}connectivity_time.pdf".format(dir_name))


def plot_largest_domain(arr, dir_name):
    fig,ax = plt.subplots()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    plt.xlabel("time [MC sweeps]")
    plt.ylabel("$L_{c}/N_{p}$")
    plt.plot(arr[:,0],arr[:,1]/1500,c=red_c, lw=3,ms=3,marker='o')
    plt.tight_layout()
    plt.savefig("{}largest_cluster_time.pdf".format(dir_name))

def plot_mean_degree(arr, dir_name):
    fig,ax = plt.subplots()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    plt.xlabel("time [MC sweeps]")
    plt.ylabel("average degree")
    plt.plot(arr[:,0],arr[:,1], c=purple_c, lw=3)
    plt.errorbar(arr[::5,0],arr[::5,1], yerr=arr[::5,2], c=purple_c, elinewidth=2, lw=0,ms=3,marker='o', capsize=5, capthick=3)
    plt.ylim((0,4))
    plt.yticks([0,1,2,3,4])
    plt.tight_layout()
    plt.savefig("{}average_degree_time.pdf".format(dir_name))

def plot_pnp_domain_distribution(arr,brr, dir_name):

    fig,ax=plt.subplots()
    # patterns = ('-', '+', 'x', '\\', '*', 'o', 'O', '.')
    plt.hist(arr,bins=range(1,np.max(arr)), fc=blue_c, edgecolor='gray', lw=1, alpha=0.7, density=True,label='largest p-domain', hatch='\\')
    plt.hist(brr,bins=range(1,np.max(arr)), fc=red_c, edgecolor='gray', lw=1, alpha=0.7, density=True,label='larget np-domain', hatch='//')
    plt.xlabel("largest domain")
    plt.ylabel("P")
    plt.locator_params(axis='y', nbins=6)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("{}pnp_largest_domain_size_distribution.pdf".format(dir_name))


def plot_largest_pnp_domain(arr,brr, dir_name):
    fig,ax = plt.subplots()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    plt.xlabel("time [MC sweeps]")
    plt.ylabel("domain size")
    plt.plot(arr[:,0],arr[:,1],c=blue_c, lw=3,ms=3,marker='o',label='p-domains')
    plt.plot(brr[:,0],brr[:,1],c=red_c, lw=3,ms=3,marker='o',label='np-domains')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("{}pnp_domain_size_time.pdf".format(dir_name))


def plot_bond_type_percent(arr,brr,dir_name):
    fig,ax = plt.subplots()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    plt.xlabel("time [MC sweeps]")
    plt.ylabel("domain size")
    plt.plot(arr[:,0],arr[:,1],c=blue_c, lw=3,ms=3,marker='o',label='p-bonds')
    plt.plot(brr[:,0],brr[:,1],c=red_c, lw=3,ms=3,marker='o',label='np-bonds')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("{}bond_type_time.pdf".format(dir_name))



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str)
    parser.add_argument('-f', type=str)
    args = parser.parse_args()

    file_name = args.f
    dir_name = args.d

    df = pd.read_pickle(dir_name+file_name)

    # Static plots

    # domain distribution at the end 
    tmax=df.time.max()
    arr= df[df.time==tmax].domain_lengths.values[0]
    plot_domain_distribution(arr, dir_name)


    # degree distribution at the end
    tmax=df.time.max()
    arr = df[df.time==tmax].degrees.values[0]
    plot_degree_distribution(arr, dir_name)


    # p and np domain lengths at the end
    tmax=df.time.max()
    arr =  df[df.time==tmax].p_domain_lengths.values[0]
    brr =  df[df.time==tmax].np_domain_lengths.values[0]
    plot_pnp_domain_distribution(arr,brr, dir_name)

    # Time dependent plots

    # plot pb(t)
    arr = df[['time','pb']].values
    plot_bond_probability(arr, dir_name)

    # plot connectiviy(t)
    arr = df[['time','node_connectivity']].values
    plot_connectivity(arr, dir_name)

    # plot largest domain (t)
    arr = df[['time','largest_domain']].values
    plot_largest_domain(arr, dir_name)

    # plot average degree (t)
    arr = df[['time', 'mean_degree','std_degree']].values
    plot_mean_degree(arr, dir_name)

    # plot largest p/np domain (t)
    arr = df[['time','largest_p_domain']].values
    brr = df[['time','largest_np_domain']].values
    plot_largest_pnp_domain(arr,brr, dir_name)

    # plot bond type percent (t)
    arr = df[['time','pbond_percent']].values
    brr = df[['time','npbond_percent']].values
    plot_bond_type_percent(arr,brr, dir_name)
