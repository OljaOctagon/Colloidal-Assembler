import numpy as np 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-time', type=int)
args = parser.parse_args()

filen = 'patch_energy_{}.bin'.format(args.time)

arr = np.fromfile(filen, dtype=np.int32)
arr = np.reshape(arr,(-1,4))

b = [ 0 if item.tolist() == [0,1,1,0]  else 1 for item in arr]
with open('binary_op.dat', 'w') as f:
    f.write("1000\n")
    f.write("Particles of frame\n")    
    for item in b:
        f.write(str(item)+"\n")


	

