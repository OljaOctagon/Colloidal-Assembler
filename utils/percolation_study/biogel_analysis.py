import glob
import pandas as pd
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.style as style
import matplotlib as mpl
import seaborn as sns

style.use('seaborn-ticks')
mpl.rcParams['font.family'] = "sans-serif"
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

from collections import defaultdict
from functools import partial

def calculate_bond_lifetime(arr, nbonds, nparticles,Ltime, file_id):

    a = np.arange(1, nparticles + 1)
    b = np.arange(1, nbonds + 1)

    f_to_number = lambda tup: (tup[0] - 1) * nbonds + (tup[1] - 1)
    f_to_tuple = lambda c: ((c // nbonds) + 1, (c % nbonds) + 1)

    arr_combinations = np.array([np.meshgrid(a, b)]).T.reshape(-1, 2)
    tuple_combinations = [(i, j) for [i, j] in arr_combinations]
    map_to_number = list(map(f_to_number, tuple_combinations))
    dict_to_number = dict(zip(tuple_combinations, map_to_number))
    dict_to_tuple = dict(zip(map_to_number, tuple_combinations))

    arr_id = np.array([
        [time, dict_to_number[(pid1, lid1)],
         dict_to_number[(pid2, lid2)]] for [time, pid1, lid1, pid2, lid2] in arr[:, :5]])

    connections = [arr_id[arr_id[:, 0] == time][:, 1:] for time in range(1, Ltime + 1)]
    # initialize bond sequences, keys: bond pairs, values: time series bonded states.
    # 1: is bonded 0 not bonded
    bond_sequences = defaultdict(partial(np.zeros, Ltime))
    print("get bond info")

    for time in range(Ltime):
        for [i, j] in connections[time]:
            bond_sequences[(i, j)][time] = 1

    # TODO: calculate number of reformed bonds over all formed bonds

    # TODO : here we can calculate recurrence: percentage of reunite in number of broken bonds
    # TODO :

    from itertools import groupby

    def len_iter(items):
        return sum(1 for _ in items)

    def consecutive_one(data):
        return [len_iter(run) for val, run in groupby(data) if val]

    print("calculate lifetime")
    # calculate bond life times
    lifetimes = []
    arr_formed = []
    arr_reformed = []
    auto_corr_heat = np.empty((0,2))

    for bond in bond_sequences.keys():
        i = bond[0]
        j = bond[1]
        # get life time
        seq = bond_sequences[i, j]
        lifetimes.extend(consecutive_one(seq))

        diffx = seq[1:]-seq[:-1]
        # get reformed bonds
        nformed = np.count_nonzero( diffx == 1 )
        arr_formed.extend(nformed + seq[0])
        arr_reformed.extend(nformed - 1)

        # get normalized autocorr of seq
        def autocorr(x):
            x = x - np.mean(x)
            result = np.correlate(x, x, mode='full')
            result = result[int(np.floor(result.size/2)):]
            rmax=np.max(result)
            result = result/rmax
            return result

        ac_seq_y = autocorr(seq)
        ac_seq_x = np.linspace(0,Ltime,Ltime)
        ac_seq = np.column_stack((ac_seq_x, ac_seq_y))
        auto_corr_heat = np.concatenate((auto_corr_heat,ac_seq))


    lifetimes = np.array(lifetimes)
    max_lifetime = np.max(lifetimes)

    Nformed = np.sum(np.array(arr_formed))
    Nreformed = np.sum(np.array(arr_reformed))
    percent_reformed = Nreformed/Nformed

    print("make histogram")
    fig, ax = plt.subplots()
    plt.xlabel("bond life time")
    plt.ylabel("P")
    nval, bins, patches = plt.hist(lifetimes,
                                edgecolor='k', bins=(np.linspace(0, Ltime, Ltime + 1) + 0.1),
                                lw=2, alpha=1,
                                density=True)

    fig, ax = plt.subplots()
    plt.yscale("log")
    x = bins[1:]

    # fit exponential
    def monoExp(x, l, t, ):
        return l * np.exp(-t * x)

    p0 = (1, 1)  # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(monoExp, x, nval, p0)
    l, t = params

    plt.plot(x, nval, c=purple_c, linestyle='dotted', label='data')
    plt.plot(x, monoExp(x, l, t), c=blue_c, linestyle='--', label='biexp. fit')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.ylim(0.00001, 1)
    plt.savefig("data_ana/{}/bond_lifetime_fit.pdf".format(file_id))

    def cumulative_monoExp(xs, l, t):
        return 1 - l * np.exp(-t * xs)

    H,x_edges, y_edges = np.histogram2d( auto_corr_heat[:,0], auto_corr_heat[:,1],
                                            bins=[500,100], density=True)

    fig,ax  = plt.subplots()
    sns.heatmap(H)
    plt.savefig("data_ana/{}/autocorr_lifetimes_heatmap.pdf")

    arr = np.array([l, t, max_lifetime, percent_reformed])
    np.savetxt("data_ana/{}/param_bond_lifetime.dat".format(file_id), arr, delimiter=',', newline='\n')
    np.savetxt("data_ana/{}/bond_lifetime.dat".format(file_id), nval, newline='\n')

if __name__ == '__main__':

    files_list = glob.glob("*/pdf1/*link.dat")
    df_params = pd.read_csv("parameters.txt", delim_whitespace=True, header=True)
    nparticles = 340
    Ltime = 2500

    for file_i in files_list:
        print("read gel")
        file_id = file_i.split('/')[0]
        arr = pd.read_csv(file_i, delim_whitespace=True).values
        nbonds = df_params['ID' == file_id]["N"]

        print("prepare data")
        # throw away links between crosslinkers on same polymer
        arr = arr[arr[:, 1] != arr[:, 3]]

        print("calculate bond life time")
        calculate_bond_lifetime(arr, nbonds, nparticles, Ltime, file_id)


