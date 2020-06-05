
import numpy as np 


Lx=1.0
alpha=np.pi/3.
ax=Lx*np.cos(alpha)
h=Lx*np.sin(alpha)

# positions

d = np.sqrt(np.power(Lx+ax,2) + np.power(h,2))/2.
c = Lx/2.
b = np.sqrt(np.power(d,2)-np.power(c,2))

A=np.zeros((6,3))
shift=np.array([5.47723,4.74342,0.05])

A[0] = np.array([0,d,0]) + shift
A[1] = np.array([b,c,0]) + shift
A[2] = np.array([b,-c,0]) + shift
A[3] = np.array([0,-d,0]) + shift
A[4] = np.array([-b,-c,0]) + shift
A[5] = np.array([-b,c,0]) + shift


A.tofile("positions_0.bin")


# orientations
offset=np.pi/3.0
Phi=np.zeros((6,5))
p=np.pi/3.0

Phi[0] = np.array([0,0,0,0,offset - 0*p])
Phi[1] = np.array([0,0,0,0,offset - 1*p])
Phi[2] = np.array([0,0,0,0,offset - 2*p])
Phi[3] = np.array([0,0,0,0,offset - 3*p])
Phi[4] = np.array([0,0,0,0,offset - 4*p])
Phi[5] = np.array([0,0,0,0,offset - 5*p])

Phi.tofile("orientations_0.bin")

