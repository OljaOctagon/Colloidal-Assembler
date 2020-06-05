import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-radius", type=float)
parser.add_argument("-delta", type=float)

p = parser.parse_args()

radius=p.radius
patch_delta=p.delta

df = pd.read_csv("patch_state.dat", header=None, delim_whitespace=True, 
	names=['dx', 'dy', 'dphi'])

df = df.drop_duplicates()

def pbc(phi):
    if phi > np.pi:
        phi = 2*np.pi - phi
    if phi < -np.pi:
        phi = -2*np.pi - phi

    return phi

df['dphi'] = df.apply(lambda x: pbc(x.dphi), axis=1)
df['dx'] = -1*df['dx']

data = df.values

Lx=1
maxr=0.1
sigma = maxr*2*Lx
delta=0.001*Lx
bins_px = np.arange(0, sigma+delta, delta)
bins_py = np.arange(-sigma,sigma+delta, delta)
bins_phi = np.arange(-np.pi,np.pi+delta,delta)

counts, edges = np.histogramdd(data, bins=[bins_px, bins_py, bins_phi])
flat_counts = counts.flatten()
N_nonzero = len(flat_counts[flat_counts>0])
volume = N_nonzero*(np.power(delta,3))
N_counts = len(flat_counts)
with open("../p_bonding_vol.dat", 'a') as f:
    f.write("{} {} {} {} {}\n".format(patch_delta,radius,volume, df.dphi.max(), df.dphi.min()))
