import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser()
parser.add_argument('-ptype', type=str,
                    choices=['dma-as1', 'dmo-s1', 'dmo-s2', 'dmo-as1'])
args = parser.parse_args()

pdict = {'dma-as1': 'double_manta_asymm_1',
         'dmo-s1': 'double_mouse_symm_1',
         'dmo-s2':  'double_mouse_symm_2',
         'dmo-as1': 'double_mouse_asymm_1'}

N = 1500
Phi = [0.01, 0.03, 0.05, 0.075, 0.08, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225,
       0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525]
Delta = [0.2, 0.3, 0.4]
Temps = ['0.01', '0.02', '0.03', '0.04, 0.05', '0.06', '0.07', '0.08', '0.09', '0.10',
         '0.11', '0.12', '0.13', '0.14', '0.15', '0.16']


len_delta = len(Delta)
len_temps = len(Temps)

fig, ax = plt.subplots(len_temps, len_delta, figsize=(10, 20))
plt.rcParams["axes.prop_cycle"] = plt.cycler(
    "color", plt.cm.viridis(np.linspace(0, 1, 23)))

for i, ti in enumerate(Temps):
    for j, delta in enumerate(Delta):
        axi = ax[i, j]
        axi.set_yscale("log")

        for phi in Phi:
            file_i = '{}_phi_{}_delta_{}_temp_{}/Energy.dat'.format(
                pdict[args.ptype], phi, delta, ti)

            #temp = float(re.findall(r"[-+]?\d*\.\d+|\d+",file_i)[3])
            arr = pd.read_csv(file_i, delim_whitespace=True).values

            arr[:, 1] = -1*arr[:, 1]
            axi.plot(arr[:, 0], arr[:, 1], lw=0.3,
                     label="$\\phi = ${}".format(phi))

        axi.set_title("$\Delta$={},T={}".format(delta, ti))
        axi.set_xlabel("sweeps")
        axi.set_ylabel("abs(U)")
        box = axi.get_position()
        axi.set_position([box.x0, box.y0, box.width * 0.65, box.height*0.8])
# plt.tight_layout()
plt.legend(loc=1, bbox_to_anchor=(0.2, -0.3), ncol=6)

plt.savefig("{}_energy.pdf".format(pdict[args.ptype]))
