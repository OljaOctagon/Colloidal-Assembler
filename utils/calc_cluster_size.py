import sys
import os
import re
import numpy as np

def calculate_all_frames():
    rootdir = os.getcwd()
    for subdir, _, files in os.walk(rootdir):
        # clear previous output
        print()
        print(subdir)
        print()
        os.system('rm {}/*1.bin'.format(subdir))
        os.system('rm {}/psi_op.dat'.format(subdir))
        for filename in files:
            if filename.endswith(".bin") and filename.startswith("positions"):
                split_name = re.split('\.|\_',filename)
                n_iter=int(split_name[1])
                os.system("cp McPoly {} ".format(subdir))
                os.system("cp para.ini {}".format(subdir))
                f_patch_energy = 'patch_energy_{}.bin'.format(n_iter)
                os.system("cp patch_energy.bin {}/{}".format(subdir, f_patch_energy))

                print(n_iter, n_iter+1)
                os.system("cd {}; ./McPoly -f {} {}; cd ..".format(subdir,
                                                                n_iter,
                                                                n_iter+1))

    os.system("for dir in ./*/; do cat $dir/psi_op.dat >> psi_op.dat; done")

def calculate_last_frames(key_1, key_2, delta):
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
            os.system('rm {}/psi_op.dat'.format(rel_path))
            files = os.listdir(rel_path) 
            # get last frame
            g = lambda filename: int(re.split('\.|\_',filename)[1])
            n_iter = max([ g(filename) for filename in files
                              if filename.endswith(".bin") and
                              filename.startswith("positions")])
            print('n_iter!!!!', n_iter)
            os.system("cp McPoly {} ".format(rel_path))
            os.system("cp para.ini.original para.ini")

            xdelta = float(rel_path.split('_')[-2])
            print(xdelta)
            os.system("gsed -i 's/xdelta/{}/g' para.ini".format(xdelta))
            os.system("cp para.ini {}".format(rel_path))
            f_patch_energy = 'patch_energy_{}.bin'.format(n_iter)
            os.system("cp patch_energy.bin {}/{}".format(rel_path, f_patch_energy))

            os.system("cd {}; ./McPoly -f {} {}; cd ..".format(rel_path, n_iter, n_iter+1))
            os.system('cat {}/psi_op.dat >> psi_op_{}.dat'.format(rel_path, key_1))


import pandas as pd
import json
def calculate_mean_psi(key_1, nmin, nmax):
    df = pd.read_csv("psi_op_{}.dat".format(key_1), delim_whitespace=True, header=None)
    df.columns = ['global', 'local', 'csize']
    dg = df[df.csize>nmin]
    dg = dg[dg.csize<nmax]
    print(dg)
    return (len(dg), dg.local.mean(), dg.local.std())

if __name__ == "__main__":
    nmin=200
    nmax=1000
    deltas=[0.2, 0.3, 0.4, 0.5]
    psi_values = {}
    for delta in deltas:
        key_1 = 'mu_0.25Energy_-5.2Asymm_patchpos_{}'.format(delta)
        key_2 = 'mu_0.25Energy_-5.2Asymm_patchpos_{}'.format(1-delta)
        calculate_last_frames(key_1,key_2,delta)
        psi_values[delta] = calculate_mean_psi(key_1, nmin, nmax)

    def default(o):
        if isinstance(o, np.int64): return int(o)  
        raise TypeError

    with open('psi_mean.json', 'w') as fp:
        json.dump(psi_values, fp, default=default)
