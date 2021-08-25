import numpy as np 
import matplotlib.pyplot as plt
from glob import glob 
import argparse 
import pandas as pd 
import re

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
Temps = ['0.01','0.04','0.05','0.07','0.09','0.10','0.11','0.12','0.13','0.14','0.15','0.16']


len_phi=len(Phi)
len_delta=len(Delta)

fig,ax=plt.subplots(len(Phi),len(Delta), figsize=(10,20))
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.viridis(np.linspace(0,1,12)))


for i,phi in enumerate(Phi):
	for j, delta in enumerate(Delta):
		axi=ax[i,j]
		axi.set_yscale("log")
	
		for temp in Temps:
			file_i=	'{}_phi_{}_delta_{}_temp_{}/Energy.dat'.format(pdict[args.ptype],phi,delta, temp)
					
			#temp = float(re.findall(r"[-+]?\d*\.\d+|\d+",file_i)[3])
			arr=pd.read_csv(file_i, delim_whitespace=True).values

			arr[:,1] = -1*arr[:,1]
			axi.plot(arr[:,0],arr[:,1],lw=0.3, label="$T = ${}".format(temp))

		axi.set_title("$\Delta$={},$\phi$={}".format(delta,phi))
		axi.set_xlabel("sweeps")
		axi.set_ylabel("abs(U)")
		box = axi.get_position()
		axi.set_position([box.x0, box.y0, box.width * 0.65, box.height*0.8])
#plt.tight_layout()
plt.legend(loc=1, bbox_to_anchor=(0.2,-0.3),ncol=6)

plt.savefig("{}_energy.pdf".format(pdict[args.ptype]))