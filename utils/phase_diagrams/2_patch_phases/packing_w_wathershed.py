import numpy as np 
import cv2
from matplotlib import pyplot as plt

def calculate_packing_pixels(filename):
    img = cv2.imread(filename, cv2.IMREAD_GRAYSCALE)
    n_white_pix = np.sum(img == 255)
    size = img.size
    packing = 1 - (n_white_pix/size)    
    return packing

import os
filelist = next(os.walk('.'))[2]

starp = []
starp.append(list(filter(lambda s: s.startswith('starbox_02'), filelist)))
starp.append(list(filter(lambda s: s.startswith('starbox_03'), filelist)))
starp.append(list(filter(lambda s: s.startswith('starbox_04'), filelist)))
#starp.append( list(filter(lambda s: s.startswith('starbox_05'), filelist))) 
starp.append( list(filter(lambda s: s.startswith('starbox_06'), filelist))) 
starp.append( list(filter(lambda s: s.startswith('starbox_07'), filelist))) 
starp.append( list(filter(lambda s: s.startswith('starbox_08'), filelist))) 

patchpos = [0.2,0.3,0.4,0.6,0.7,0.8]

for i in range(len(starp)):
    packing = np.array([ calculate_packing_pixels(fn) for fn in starp[i]])
    packing_mean = np.round(np.mean(packing), decimals=3)
    packing_std = np.round(np.std(packing), decimals=3)
    
    print(patchpos[i], packing_mean, packing_std)
