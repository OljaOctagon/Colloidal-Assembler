import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter

sns.set(style="white")
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def gr_per_run(data, Delta, R, delta_time):
    a=Delta
    G = np.zeros(int((R-a)/Delta))
    start_time=data['time'].values[0]
    end_time=data['time'].values[-1]
    norm_time=(end_time-start_time)/float(delta_time)

    print(start_time)
    print(end_time)
    print('norm', norm_time)

    for t in range(start_time, end_time+delta_time, delta_time):
        print(t)
        is_time = data["time"] == t
        P_t = data[is_time]
        is_id0 = P_t["id"] == 0
        Np = P_t[is_id0]["N"].values
        lx = P_t[is_id0]["Lx"].values
        ly = P_t[is_id0]["Ly"].values

        center_x=P_t["X"]
        center_y=P_t["Y"]

        particle_distx = np.reshape(center_x, (len(center_x), 1))
        particle_distx = particle_distx - particle_distx.transpose()

        particle_disty = np.reshape(center_y, (len(center_y), 1))
        particle_disty = particle_disty - particle_disty.transpose()

        particle_distx = particle_distx - lx*np.rint(particle_distx/lx)
        particle_disty = particle_disty - ly*np.rint(particle_disty/ly)

        particle_dist=np.sqrt(np.power(particle_distx,2) + np.power(particle_disty,2))

        rho=float(Np)/(lx*ly)
        S=int((R-a)/Delta)

        for i in range(1,S):
            rb = a+i*Delta
            ra = rb-Delta        
            A = np.pi*(rb*rb - ra*ra)      
            normed = rho*A
            gr = 0
            for id in range(int(Np)):
                dx= particle_dist[id]
                na=len(dx[np.where((dx > ra) & (dx < rb))])      
                gr = gr + na/normed
            G[i] = G[i] + gr/float(Np)

    G=G/float(norm_time)
    return G

def check_file(fname):
    try:
        open(fname,"r")
        return 1
    except IOError:
        print("Error: File doesn't seem to exist.")
        return 0


import argparse
import glob 

if __name__ == '__main__':

    Delta=0.01
    R=10
    delta_time=20000
    print('hi')
    parser = argparse.ArgumentParser()
    parser.add_argument('-dirs', default='mu_0.4Energy_9.2symm_patchpos_*')
    parser.add_argument('-fname', default='center_of_mass_stars.dat')
    args = parser.parse_args()

    dirs = glob.glob("{}*".format(args.dirs))
    print(dirs)
    G = np.zeros(int((R-Delta)/Delta))
    for dir_i in dirs:
        print(dir_i)
        filen = "{}/{}".format(dir_i,args.fname)
        f_exist = check_file(filen)
        if f_exist == 1:
            df = pd.read_csv(filen, delim_whitespace=True, header=None)
            df.columns=["time","N", "Lx", "Ly", "id", "X", "Y"]
            G = G + gr_per_run(df, Delta, R, delta_time)

    G = G/len(dirs)
    print(len(dirs))
    np.savetxt("{}_gr2d.dat".format(args.dirs), G, delimiter=' ', newline='\n', fmt='%.5f')            
