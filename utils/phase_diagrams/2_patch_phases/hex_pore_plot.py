import numpy as np
import matplotlib.pyplot as plt 


def area_hexagon(Lx, patch_pos):
    ''' The bigger pores are hexagon and the smaller triangles of 
    l=(fabs(patch_pos - (1-patch_pos)).'''
    pore_p = np.fabs(2*patch_pos - 1)
    pore_size = pore_p*Lx
    pore_area = (3*np.sqrt(3)/2)*np.power(pore_size,2)
    return pore_area


def area_triangle(Lx, patch_pos):
    pore_p = np.fabs(2*patch_pos - 1)
    pore_size = pore_p*Lx
    pore_area = (np.sqrt(3)/(6)*np.power(pore_size,2))
    return pore_area



def area_rhombi(Lx,patch_pos):
    ''' the unit cell has one rhombic pore. 
    Its area is a function of the patch size. '''
    pore_p = np.fabs(2*patch_pos - 1)
    pore_size = pore_p*Lx
    pore_area = pore_size*np.sin(np.pi/3.)
    return pore_area


import matplotlib.pyplot as plt 
from matplotlib import rc 
import seaborn as sns 

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
sns.set(style='whitegrid')

#-----------

Lx = 1
patch_pos = np.arange(0,1.05,0.1)
area_h = area_hexagon(Lx, patch_pos)
area_r = area_rhombi(Lx,patch_pos)
area_t = area_triangle(Lx, patch_pos)

fig, ax = plt.subplots()
ax.plot(patch_pos,area_h, lw=2, marker='h', c='#4d1ccb', 
    label='$\displaystyle \\textrm{hexagonal pore}$',
    alpha=0.6)

ax.plot(patch_pos,area_t, lw=2, marker='^', c='g', 
    label='$\displaystyle \\textrm{triangular pore}$',
    alpha=0.6)

ax.plot(patch_pos,area_r, lw=2, marker='D', c='r', 
    label='$\displaystyle \\textrm{rhombic pore}$',
    alpha=0.6)



ax.set_xlim([0,1])
#ax.set_ylim([0.0,1.03])
ax.legend(loc='upper center', fancybox=True, prop={'size': 13})
ax.tick_params(axis='both', which='major', labelsize=18)
ax.set_xlabel("$ \displaystyle \Delta$", size=22)
ax.set_ylabel("$ \displaystyle \\textrm{pore area} [1/l^2]$", size=22)
plt.tight_layout()
plt.savefig("pore_area_function.png", dpi=400)

