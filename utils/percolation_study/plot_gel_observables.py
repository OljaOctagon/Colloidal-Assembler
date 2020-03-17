import argparse 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.style as style
style.use('seaborn-poster') 
mpl.rcParams['font.family'] = "sans-serif"
sns.set_context('poster')
plt.rcParams['axes.axisbelow'] = True


def plot_domain_distribution(arr):
    pass

def plot_degree_distribution(arr):
    pass

def plot_bond_probablity(arr):
    pass

def plot_connectivity(arr):
    pass

def plot_largest_domain(arr):
    pass 

def plot_mean_degree(arr):
    pass

def plot_pnp_domain_distribution(arr,brr):
    pass

def plot_largest_pnp_domain(arr,brr):
    pass


if __name__ == "__main__":i

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str)
    parser.add_argument('-f', type=str)
    args = parser.parse_args()

    file_name = args.f
    dir_name = args.p

    df = pd.read_picke(dir_name+file_name)

    # Info: columns:
    columns = ['run_id',
            'topology',
            'delta','T',
            'phi','N_particles',
            'time',
            'domain_lengths',
            'p_domain_lengths',
            'np_domain_lengths',
            'degrees',
            'mean_degree',
            'std_degree',
            'pb',
            'node_connectivity',
            'N3_loop_percent',
            'largest_domain',
            'largest_pdomain',
            'largest_np_domain',
            'pbond_percent',
            'npbond_percent']


    # Static plots

    # domain distribution at the end 
    tmax=df.time.max()
    arr= df[df.time==tmax].domain_lengths.values
    plot_domain_distribution(arr)


    # degree distribution at the end
    tmax=df.time.max()
    arr = df[df.time==tmax].degrees.values
    plot_degree_distribution(arr)


    # p and np domain lengths at the end
    tmax=df.time.max()
    arr =  df[df.time==tmax].p_domain_lengths.values
    brr =  df[df.time==tmax].np_domain_lengths.values
    plot_pnp_domain_distribution(arr,brr)

    # Time dependent plots

    # plot pb(t)
    arr = df[['time','pb']].values
    plot_bond_probability(arr)

    # plot connectiviy(t)
    arr = df[['time','connectivity']].values
    plot_connectivity(arr)

    # plot largest domain (t)
    arr = df[['time','largest_domain']].values
    plot_largest_domain(arr)

    # plot average degree (t)
    arr = df[['time', 'mean_degree','std_degree']]
    plot_mean_degree(arr)

    # plot largest p/np domain (t)
    arr = df[['time','largest_p_domain']].values
    brr = df[['time','largest_np_domain']].values
    plot_largest_pnp_domain(arr,brr)
