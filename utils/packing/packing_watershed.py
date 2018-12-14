import numpy as np 
import cv2
from matplotlib import pyplot as plt
import argparse 
def calculate_packing_pixels(filename):
    img = cv2.imread(filename, cv2.IMREAD_GRAYSCALE)
    n_white_pix = np.sum(img == 255)
    size = img.size
    packing = 1 - (n_white_pix/size)    
    return packing

import os
filelist = next(os.walk('.'))[2]

parser = argparse.ArgumentParser()
parser.add_argument('-key')
args = parser.parse_args()

filesplits = []
patchpos = [0.2,0.3,0.4,0.5]

for patch in patchpos:
    filesplits.append(list(filter(lambda s: s.startswith("{}_{}".format(args.key,patch)), filelist)))

with open("measured_packing_{}.dat".format(args.key), 'w') as fopen:
    for files, patch in zip(filesplits, patchpos): 
        packing = np.array([ calculate_packing_pixels(fn) for fn in files ])
        packing_mean = np.round(np.mean(packing), decimals=3)
        packing_std = np.round(np.std(packing), decimals=3)

        fopen.write("{} {} {}\n".format(patch, packing_mean, packing_std))
        fopen.write("{} {} {}\n".format(1-patch, packing_mean, packing_std))
