# coding: utf-8
import numpy as np
import pandas as pd 
df = pd.read_csv("Acceptances.txt", header=None, delim_whitespace=True, 
    names=['time', 'nmoves', 'nbroken', 'phi_i', 'patch_ix', 'patch_iy', 'phi_j','patch_jx', 'patch_jy','wx', 'wy'])
df['dpatch_x'] = df.patch_ix - df.patch_jx
df['dpatch_y'] = df.patch_iy - df.patch_jy
df['dphi'] = np.arctan2(df.wx,df.wy)

data = df[['dpatch_x', 'dpatch_y', 'dphi']].values


def omega_y_max_pos(y,delta):
    sigma=0.1
    return np.arctan2(np.sqrt(np.power(sigma,2) - np.power(y,2)), delta + y)

def omega_y_max_neg(y,delta):
    sigma=0.1
    return np.arctan2(np.sqrt(np.power(sigma,2) - np.power(y,2)), delta - y)


def omega_x_max(x,delta):
    return np.arctan2(x,delta)



import matplotlib.pyplot as plt 
import seaborn as sns 

fig,ax = plt.subplots()
plt.figure(figsize=(7,14))  

plt.xlabel("patch dx", size=16)
plt.ylabel("patch dy", size=16)

plt.xlim([0,0.1])
plt.ylim([-0.1,0.1])

plt.hist2d(data[:,0], data[:,1], bins=[np.arange(0,0.1+0.001,0.001), np.arange(-0.1,0.1+0.001,0.001)])
ax.tick_params(labelsize=12)
plt.xticks([0,0.025,0.05,0.075,0.1])
plt.yticks([-0.1,-0.05,0,0.05,0.1])
ax.tick_params(labelsize=30)
plt.tight_layout()
plt.savefig("patch_d_histo_xy.pdf")


fig,ax = plt.subplots()
plt.xlabel("patch dx", size=16)
plt.ylabel("$\omega$", size=16)

plt.xlim([0,0.1])

plt.hist2d(data[:,0], data[:,2], bins=[np.arange(0,0.1+0.001,0.001), np.arange(-0.6,0.6,0.001)])
ax.tick_params(labelsize=12)
plt.xticks([0,0.025,0.05,0.075,0.1])
ax.tick_params(labelsize=12)
plt.tight_layout()
plt.savefig("patch_d_histo_xw.pdf")

fig,ax = plt.subplots()
plt.xlabel("patch dy", size=16)
plt.ylabel("$\omega$", size=16)

plt.xlim([-0.1,0.1])
plt.hist2d(data[:,1], data[:,2], bins=[np.arange(-0.1,0.1+0.001,0.001), np.arange(-0.6,0.6,0.001)])

#y1=np.arange(0,0.1,0.001)
#y2=np.arange(-0.1,0,0.001)
#plt.plot(y1, omega_y_max_pos(y1,0.5), c='red', lw=3)
#plt.plot(y2, omega_y_max_pos(y2,0.5), c='yellow', lw=3)
#plt.plot(y1, -omega_y_max_pos(-y1,0.5), c='blue', lw=3)
#plt.plot(y2, -omega_y_max_pos(-y2,0.5), c='yellow', lw=3)

ax.tick_params(labelsize=12)
plt.xticks([-0.1,-0.05,0,0.05,0.1])
ax.tick_params(labelsize=12)
plt.tight_layout()
plt.savefig("patch_d_histo_yw.pdf")






