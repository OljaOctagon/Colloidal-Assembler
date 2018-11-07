import numpy as np 
import pandas as pd 

def area_particles(Lx):
    '''Are that particls inhabit is given by 6 full particles and 6 half ones.
    That results in 9 full 60 degree rhombi. '''
    area_rhombi = Lx*np.sin(np.pi/3.) 
    return 9*area_rhombi

def area_pores(Lx, patch_pos):
    ''' The pores are hexagon and triangles of l=(fabs(patch_pos - (1-patch_pos)).
    In total there is one full hexagon, 6 1/3 hexagon, and 6 triangles. 
    That equals the area of 4 full hexagons '''
    pore_p = 2*patch_pos - 1
    pore_length = pore_p*Lx
    hexagon = (3*np.sqrt(3)/2)*np.power(pore_length,2)
    return 4*hexagon


Lx = 1
area_pa = area_particles(Lx)
patch_pos = np.arange(0,1.05,0.05)
area_po = area_pores(Lx, patch_pos)
total_area = area_pa + area_po
packing = area_pa/total_area


exp_packing = pd.read_csv("measured_packing_starboxes.dat", header=None,
    delim_whitespace=True).values

import matplotlib.pyplot as plt 
from matplotlib import rc 
import seaborn as sns 

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

sns.set(style='whitegrid')
fig,ax = plt.subplots()
plt.plot(patch_pos,packing, lw=2, marker='s', c='#4d1ccb', 
    label='theoretical maximum packing',
    alpha=0.9)

plt.errorbar(exp_packing[:,0], exp_packing[:,1], yerr=exp_packing[:,2], 
    lw=0, c='r', elinewidth=1,
    capsize=3, capthick=2, 
    label='measured packing for $\displaystyle r_{p} = 0.1$',
    alpha=1.0)

plt.xlim([0,1])
plt.legend(loc='upper center', fancybox=True, prop={'size': 12}, bbox_to_anchor=(0.5, 1.1), ncol=2)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlabel("$ \displaystyle \Delta$", size=22)
plt.ylabel("$ \displaystyle \phi$", size=22)
plt.tight_layout()
plt.savefig("theoretical_packing_fraction_starboxes.pdf")
plt.show()

