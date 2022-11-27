import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import argparse
import pandas as pd
import re
from matplotlib.ticker import NullFormatter

parser = argparse.ArgumentParser()
parser.add_argument('-ptype', type=str,
                    choices=['dma-as1', 'dmo-s1', 'dmo-s2', 'dmo-as1'])

parser.add_argument(
    "-b", type=str, choices=["batch1", "batch2", "batch3", "batch4"])
args = parser.parse_args()

pdict = {'dma-as1': 'double_manta_asymm_1',
         'dmo-s1': 'double_mouse_symm_1',
         'dmo-s2':  'double_mouse_symm_2',
         'dmo-as1': 'double_mouse_asymm_1'}

N = 1500
dict_phi = {"batch1": [0.01, 0.03, 0.05, 0.075, 0.08],
            "batch2": [0.1, 0.125, 0.15, 0.175, 0.2, 0.225],
            "batch3": [0.25, 0.275, 0.3, 0.325, 0.35, 0.375],
            "batch4": [0.4, 0.425, 0.45, 0.475, 0.5, 0.525]}

Phi = dict_phi[args.b]

Delta = [0.2, 0.3, 0.4]
Temps = ['0.01', '0.02', '0.03', '0.04', '0.05', '0.06', '0.07', '0.08', '0.09', '0.10',
         '0.11', '0.12', '0.13', '0.14', '0.15', '0.16']

len_phi = len(Phi)
len_temps = len(Temps)

fig, ax = plt.subplots(len_temps, len_phi, figsize=(10, 20))
plt.rcParams["axes.prop_cycle"] = plt.cycler(
    "color", plt.cm.viridis(np.linspace(0, 1, 23)))

for delta in Delta:
    for i, ti in enumerate(Temps):
        for j, phi in enumerate(Phi):
            axi = ax[i, j]
            axi.set_yscale("log")
            for ri in range(1, 9):
                print(ri)
                file_i = '{}_phi_{}_delta_{}_temp_{}_run_{}/Energy.dat'.format(
                    pdict[args.ptype], phi, delta, ti, ri)

                print(file_i)
                arr = pd.read_csv(file_i, delim_whitespace=True).values
                arr[:, 1] = -1*arr[:, 1]
                axi.plot(arr[:, 0], arr[:, 1]/1000, lw=0.3,
                         label="run id = {}".format(ri))

            axi.set_title("$\phi$={},T={}".format(phi, ti), size=5)
            axi.set_xlabel("sweeps", size=5)
            axi.set_ylabel("abs(U)", size=5)
            axi.tick_params(axis='both', which='major', labelsize=5)

            axi.yaxis.set_major_formatter(NullFormatter())
            axi.yaxis.set_minor_formatter(NullFormatter())

            axi.set_ylim((0, 5))
            box = axi.get_position()
            axi.set_position(
                [box.x0, box.y0, box.width * 0.65, box.height*0.8])
    plt.tight_layout()
    plt.legend(loc=1, bbox_to_anchor=(0.2, -0.3), ncol=6)
    plt.savefig("{}_{}_energy.pdf".format(delta, pdict[args.ptype]))
