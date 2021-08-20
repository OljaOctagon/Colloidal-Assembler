import numpy as np 
import matplotlib.pyplot as plt
from glob import glob 
import argparse 
import pandas as pd 

parser=argparse.ArgumentParser()
parser.add_argument('-ptype',type=str, choices=['dma-as1', 'dmo-s1', 'dmo-s2', 'dmo-as1'])
args = parser.parse_args()


pdict = {'dma-as1': 'double_manta_asymm_1',
'dmo-s1': 'double_mouse_symm_1',
'dmo-s2':  'double_mouse_symm_2',
'dmo-as1': 'double_mouse_asymm_1'}

N=1500
Phi = [0.01,0.05,0.1,0.2,0.3,0.4,0.5]
Delta = [0.2,0.3,0.4,0.5]

len_phi=len(Phi)
len_delta=len(Delta)

fig,ax=plt.subplots(len(Phi),len(Delta), figsize=(10,20))

for i,phi in enumerate(Phi):
	for j, delta in enumerate(Delta):
		axi=ax[i,j]
		axi.set_yscale("log")
		filelist=glob('{}_phi_{}_delta_{}_temp_*/Energy.dat'.format(pdict[args.ptype],phi,delta))

		for file_i in filelist:
			arr=pd.read_csv(file_i, delim_whitespace=True).values

			arr[:,1] = -1*arr[:,1]
			axi.plot(arr[:,0],arr[:,1],lw=0.1)

		axi.set_xlabel("MC sweeps")
		axi.set_ylabel("abs(U)")

plt.tight_layout()
plt.savefig("{}_energy.pdf".format(pdict[args.ptype]))