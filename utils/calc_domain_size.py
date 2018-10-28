import sys
import os
import re
import numpy as np

def calculate_last_frames(key_1, key_2, delta, filen):
    rootdir = os.getcwd()
    subdirs = os.listdir(rootdir)
    print(subdirs)
    for subdir in subdirs:
        rel_path = subdir.split('/')[-1]
        if rel_path.startswith(key_1) or rel_path.startswith(key_2):
            # clear previous output
            print()
            print('worked!', rel_path)
            print()
            os.system('rm {}/*1.bin'.format(rel_path))
            os.system('rm {}/{}'.format(rel_path, filen))
            files = os.listdir(rel_path) 
            # get last frame
            g = lambda filename: int(re.split('\.|\_',filename)[1])
            n_iter = max([ g(filename) for filename in files
                              if filename.endswith(".bin") and
                              filename.startswith("positions")])
            os.system("cp McPoly {} ".format(rel_path))
            os.system("cp para.ini.original para.ini")

            xdelta = float(rel_path.split('_')[-2])
            print(xdelta)
            os.system("gsed -i 's/xdelta/{}/g' para.ini".format(xdelta))
            os.system("cp para.ini {}".format(rel_path))
            f_patch_energy = 'patch_energy_{}.bin'.format(n_iter)
            os.system("cp patch_energy.bin {}/{}".format(rel_path, f_patch_energy))
            os.system("cd {}; ./McPoly -f {} {}; cd ..".format(rel_path, n_iter, n_iter+1))
            os.system("cp domain_size.py {}; cd {}; python domain_size.py".format(rel_path, rel_path))
            os.system('cat {}/domain_sizes_results.dat >> domain_sizes_results_{}.dat'.format(rel_path, key_1))


import pandas as pd
import json
def calculate_mean_domain_sizes(key_1, nmin, nmax):
    df = pd.read_csv("domain_sizes_results_{}.dat".format(key_1), delim_whitespace=True, header=None)
    df.columns = ['mean', 'std', 'bonding_type', 'csize']
    dg = df[df.csize>nmin]
    dg = dg[dg.csize<nmax]
    print(dg)
    return (len(dg), dg.mean.mean(), dg.std.std())

if __name__ == "__main__":
    nmin=200
    nmax=1000
    deltas=[0.5]
    domain_values = {}
    filen="bond_domains.dat"
    for delta in deltas:
        key_1 = 'mu_0.25Energy_-5.2Asymm_patchpos_{}'.format(delta)
        key_2 = 'mu_0.25Energy_-5.2Asymm_patchpos_{}'.format(1-delta)
        calculate_last_frames(key_1,key_2,delta, filen)
        domain_values[delta] = calculate_mean_domain_sizes(key_1, nmin, nmax)

    def default(o):
        if isinstance(o, np.int64): return int(o)  
        raise TypeError

    with open('domain_sizes_mean.json', 'w') as fp:
        json.dump(domain_values, fp, default=default)
