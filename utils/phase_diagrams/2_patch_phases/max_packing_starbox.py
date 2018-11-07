import numpy as np 
import pandas as pd 

def area_particles(Lx):
    '''The unit cell is given by 6 full particles and 6 half ones.
    That results in 9 full 60 degree rhombi. '''
    area_rhombi = Lx*np.sin(np.pi/3.) 
    return 9*area_rhombi

def area_pores(Lx, patch_pos):
    ''' The pores are hexagon and triangles of 
    l=(fabs(patch_pos - (1-patch_pos)). In total there is 
    one full hexagon, 6 1/3 hexagon, and 6 triangles. 
    That equals the area of 4 full hexagons '''
    pore_p = 2*patch_pos - 1
    pore_length = pore_p*Lx
    hexagon = (3*np.sqrt(3)/2)*np.power(pore_length,2)
    return 4*hexagon


def area_particles_jm(Lx):
    ''' the particles of the unit cell are 2 full rhombi'''
    area_rhombi = Lx*np.sin(np.pi/3.) 
    return 2*area_rhombi

def area_pores_jm(Lx,patch_pos):
    ''' the unit cell has one rhombic pore. 
    Its area is a function of the patch size. '''
    pore_p = 1 - 2*np.minimum(patch_pos,1-patch_pos)
    pore_length = pore_p*Lx
    rhombi_p = pore_length*np.sin(np.pi/3.)
    return rhombi_p

# starboxes
#--------------------

Lx = 1
area_pa = area_particles(Lx)
patch_pos = np.arange(0,1.05,0.05)
area_po = area_pores(Lx, patch_pos)
total_area = area_pa + area_po
packing = area_pa/total_area

exp_file = "measures_packing_starboxes.dat"

exp_packing = pd.read_csv(exp_file, header=None,
    delim_whitespace=True).values

# janus mice open lattice

#----------------------------

Lx = 1
area_pa_jm = area_particles_jm(Lx)
patch_pos_jm = np.arange(0,1.05,0.05)
area_po_jm = area_pores_jm(Lx, patch_pos)
total_area_jm = area_pa_jm + area_po_jm
packing_jm = area_pa_jm/total_area_jm

#exp_file = "measured_packing_jm.dat"
#exp_packing = pd.read_csv(exp_file, header=None,
#    delim_whitespace=True).values

#--------------------------------------------------------

import matplotlib.pyplot as plt 
from matplotlib import rc 
import seaborn as sns 

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
sns.set(style='whitegrid')

#-----------

fig, ax = plt.subplots()
ax.plot(patch_pos,packing, lw=2, marker='s', c='#4d1ccb', 
    label='theoretical maximum packing janus manta',
    alpha=0.9)

ax.errorbar(exp_packing[:,0], exp_packing[:,1], yerr=exp_packing[:,2], 
    lw=0, c='r', elinewidth=1,
    capsize=3, capthick=2, 
    label='measured packing for $\displaystyle r_{p} = 0.1$',
    alpha=1.0)

ax.set_xlim([0,1])
ax.set_ylim([0.4,1.03])
ax.legend(loc='upper center', fancybox=True, prop={'size': 12}, bbox_to_anchor=(0.5, 1.1), ncol=2)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xlabel("$ \displaystyle \Delta$", size=22)
ax.set_ylabel("$ \displaystyle \phi$", size=22)

ax.plot(patch_pos,packing_jm, lw=2, marker='o', c='#36d154',
    label='theoretical maximum packing janus mice',
    alpha=0.9)


plt.tight_layout()
plt.savefig("theoretical_packing_fraction_starboxes.pdf")

#------------------------------------------

fig, ax = plt.subplots()
ax.plot(patch_pos,packing_jm, lw=2, marker='o', c='#36d154',
    label='theoretical maximum packing janus mice',
    alpha=0.9)

ax.set_xlim([0,1])
ax.set_ylim([0.4,1.03])
ax.legend(loc='upper center', fancybox=True, prop={'size': 12}, bbox_to_anchor=(0.5, 1.1), ncol=2)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xlabel("$ \displaystyle \Delta$", size=22)
ax.set_ylabel("$ \displaystyle \phi$", size=22)

plt.tight_layout()
plt.savefig("theoretical_packing_fraction_janus_mice.pdf")

plt.show()

