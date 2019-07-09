import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 

df = pd.read_csv("Acceptances.txt", header=None, delim_whitespace=True, 
	names=['time', 'nmoves', 'nbroken', 'phi_i', 'patch_ix', 'patch_iy', 'phi_j','patch_jx', 'patch_jy','wx', 'wy'])


df['dpatch_x'] = df.patch_ix - df.patch_jx
df['dpatch_y'] = df.patch_iy - df.patch_jy
df['dphi'] = np.arctan2(df.wx,df.wy)

data = df[['dpatch_x', 'dpatch_y', 'dphi']].values

print(np.max(data[:,0]))
print(np.min(data[:,0]))
print(np.max(data[:,1]))
print(np.min(data[:,1]))

print(np.min(data[:,2]))
print(np.max(data[:,2]))

Lx=1
sigma = 0.1*Lx
delta=0.0001*Lx
bins_px = np.arange(0, sigma+delta, delta)
bins_py = np.arange(-sigma,sigma+delta, delta)
bins_phi = np.arange(-np.pi,np.pi+delta,delta)

counts, edges = np.histogramdd(data, bins=[bins_px, bins_py, bins_phi])
flat_counts = counts.flatten()
N_nonzero = len(flat_counts[flat_counts>0])

volume = N_nonzero*(np.power(delta,3))

N_counts = len(flat_counts)
print( N_nonzero, N_counts, N_nonzero/N_counts, volume)
