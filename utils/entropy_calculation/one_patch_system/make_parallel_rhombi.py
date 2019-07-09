
import numpy as np 

Lx=1.0
alpha=np.pi/3.
ax=Lx*np.cos(alpha)
h=Lx*np.sin(alpha)

# positions
d = Lx;
A=np.zeros((2,3))
#shift=np.array([5.47723,4.74342,0.05])
shift = np.array([2,2,0.05])

A[0] = np.array([d*(np.sqrt(3)/2.),d*0.5,0]) + shift
A[1] = np.array([0,0,0]) + shift

A.tofile("positions_0.bin")

# orientations
Phi=np.zeros((2,5))
#p=np.pi
p1=np.pi/6.
p2 = np.pi/6. + np.pi

Phi[0] = np.array([0,0,0,0,1*p1])
Phi[1] = np.array([0,0,0,0,1*p2])


Phi.tofile("orientations_0.bin")

