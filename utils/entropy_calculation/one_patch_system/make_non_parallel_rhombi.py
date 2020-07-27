import argparse  
import numpy as np 


parser = argparse.ArgumentParser()
parser.add_argument("-delta", type=float)

p = parser.parse_args()

delta=p.delta

Lx=1.0
alpha=np.pi/3.
ax=Lx*np.cos(alpha)
h=Lx*np.sin(alpha)

# positions
d = Lx;
A=np.zeros((2,3))


shift = np.array([2,2,0.05])


var_delta = 1-2*delta
#shift_sign = np.sign(var_delta)
#var_delta = np.abs(var_delta) 

delta_arr = np.array([0,var_delta,0])
A[0] = np.array([d*(np.sqrt(3)/2.),0,0]) + shift +  delta_arr 
A[1] = np.array([0,0,0]) + shift

A.tofile("positions_1.bin")

# orientations
Phi=np.zeros((2,5))
#p=np.pi
p1=  3*np.pi/6
p2 = np.pi/6  
Phi[0] = np.array([0,0,0,0,1*p1])
Phi[1] = np.array([0,0,0,0,1*p2])


Phi.tofile("orientations_1.bin")

patch_e = np.zeros(3000).astype(np.int32)
patch_e[:8] = np.array([1,0,1,1,1,1,1,0])
patch_e.tofile("patch_energy_1.bin")


