import numpy as np
import glob
import pandas as pd
import argparse
from spanning_func import get_spanning
import gel_tools as gt
import re
from scipy.optimize import curve_fit
import pandas as pd
import configparser
import matplotlib.pyplot as plt
from os.path import exists
import networkx as nx
from collections import defaultdict


def generator_from_fsys(fsys_iterator):

    for dir_i in fsys_iterator:
        config = configparser.ConfigParser()
        config.read('{}para.ini'.format(dir_i))

        N = int(config['System']['Number_of_Particles'])
        phi = float(config['System']['Packing_Fraction'])
        temperature = float(config['System']['Temperature'])
        ptype = config['Rhombus']['rhombus_type']
        delta = config['Rhombus']['patch_delta']
        patch_size = config['Rhombus']['patch_size']

        pos_files = glob.glob('{}positions_*.bin'.format(dir_i))

        # get the last value from the string
        def g(x): return int(re.findall(r'\d+', x)[-1])

        mc_times = list(map(g, pos_files))

        last_time = np.max(mc_times)

        pos_file = "{}positions_{}.bin".format(dir_i, last_time)
        pos = np.fromfile(pos_file)
        pos = np.reshape(pos, (-1, 3))
        pos = pos[:, :2]

        box_file = "{}Box_{}.bin".format(dir_i, last_time)
        box = np.fromfile(box_file)

        yield ptype, phi, temperature, delta, last_time, pos, box


if __name__ == '__main__':

    # read data either through files system via glob or via db
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, choices=['fsys'])
    parser.add_argument('-run_id', type=str)
    parser.add_argument('-spanning', type=bool, default=True)
    parser.add_argument('-degree', type=bool, default=True)
    parser.add_argument('-psi', type=bool, default=True)

    args = parser.parse_args()

    gen_fsys = generator_from_fsys(glob.glob("double*/double*/"))
    gen_dict = {'fsys': gen_fsys}
    columns = ['id', 'ptype', 'delta', 'phi', 'temperature',
               'current_time']

    print("Calculating for all double* files in dir ")
    print("spanning: ", args.spanning)
    print("degree: ", args.degree)
    print("psi", args.psi)

    if args.spanning == True:
        columns.append('frac_largest')
        columns.append('frac_largest_virtual')

    if args.degree == True:
        columns.append('frac_degree_0')
        columns.append('frac_degree_1')
        columns.append('frac_degree_2')
        columns.append('frac_degree_3')
        columns.append('frac_degree_3')

    if args.psi == True:
        columns.append("psi_all")
        columns.append("psi_largest")

    df = pd.DataFrame(columns=columns)
    gen = gen_dict[args.input]

    for ptype, phi, temperature, delta, last_time, pos, box in gen:
        new_results = {}
        temp_str = "{0:.2f}".format(temperature)
        dir_name = "{}/{}_phi_{}_delta_{}_temp_{}".format(
            ptype, ptype, phi, delta, temp_str)
        file_name = "{}/patch_network.dat".format(dir_name)

        print("Evaluating data of {}".format(dir_name))
        subdir_name = "{}_phi_{}_delta_{}_temp_{}".format(
            ptype, phi, delta, temperature)

        new_results['id'] = subdir_name
        new_results['ptype'] = ptype
        new_results['delta'] = delta
        new_results['phi'] = phi
        new_results['temperature'] = temperature
        new_results['current_time'] = last_time

        if exists(file_name):

            connections = gt.read_bonds(file_name)[-1]
            connections = connections[:, :2]
            G = nx.Graph()
            G.add_edges_from(connections)

            if args.spanning == True:
                frac_largest, virtual_frac_largest = get_spanning(pos,
                                                                  box,
                                                                  connections,
                                                                  G,)
                new_results['frac_largest'] = frac_largest
                new_results['frac_largest_virtual'] = virtual_frac_largest

            if args.degree == True:
                degree_sequence = sorted(
                    (d for n, d in G.degree()), reverse=True)
                dmax = max(degree_sequence)

                degrees, degrees_N = np.unique(degree_sequence,
                                               return_counts=True)

                N_particles = len(pos)
                degrees_f = defaultdict(lambda: 0)
                for di, dn in zip(degrees, degrees_N):
                    degrees_f[di] = dn/N_particles

                new_results['frac_degree_0'] = degrees_f[0]
                new_results['frac_degree_1'] = degrees_f[1]
                new_results['frac_degree_2'] = degrees_f[2]
                new_results['frac_degree_3'] = degrees_f[3]
                new_results['frac_degree_4'] = degrees_f[4]

            if args.psi == True:
                file_psi = "{}/psi_op.dat".format(dir_name)
                dg = pd.read_csv(file_psi, delim_whitespace=True,
                                 names=['psi_all', 'psi_largest', 'N_largest'])

                arr = dg.values
                new_results['psi_all'] = arr[-1, 0]
                new_results['psi_largest'] = arr[-1, 1]

        else:
            print("{}: doesn't exist".format(file_name))
            if args.spanning == True:
                new_results['frac_largest'] = np.nan
                new_results['frac_largest_virtual'] = np.nan

            if args.degree == True:
                new_results['frac_degree_0'] = np.nan
                new_results['frac_degree_1'] = np.nan
                new_results['frac_degree_2'] = np.nan
                new_results['frac_degree_3'] = np.nan
                new_results['frac_degree_4'] = np.nan

            if args.psi == True:
                new_results['psi_all'] = np.nan
                new_results['psi_largest'] = np.nan

        df = df.append(new_results, ignore_index=True)

    df.to_pickle("results_gel_{}.pickle".format(
        args.run_id))
