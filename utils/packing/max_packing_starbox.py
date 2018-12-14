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
patch_pos = np.arange(0,1.01,0.01)
area_po = area_pores(Lx, patch_pos)
total_area = area_pa + area_po
packing = area_pa/total_area

exp_file = "measures_packing_starboxes.dat"

exp_packing = pd.read_csv(exp_file, header=None,
    delim_whitespace=True).values

exp_file_center = "measures_packing_starboxes_center.dat"
exp_packing_center = pd.read_csv(exp_file_center, header=None,
    delim_whitespace=True).values

exp_file = "measured_packing_starboxes.dat"
dma = pd.read_csv(exp_file, header=None,
    delim_whitespace=True).values
# janus mice open lattice

#----------------------------

Lx = 1
area_pa_jm = area_particles_jm(Lx)
patch_pos_jm = np.arange(0,1.01,0.01)
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


# dma-as2
#-----------

fig, ax = plt.subplots()
ax.plot(patch_pos,packing, lw=2, c='#336699', 
    label='sticky limit, triangular-pore phase',
    alpha=0.9)

ax.errorbar(exp_packing[:,0], exp_packing[:,1], yerr=exp_packing[:,2], 
    lw=0, c='#003366', elinewidth=1,
    capsize=3, capthick=2, 
    label='$\displaystyle r_{p} = 0.1$, triangular-pore phase',
    alpha=1.0)

print(exp_packing_center)
ax.errorbar(exp_packing_center[:,0], exp_packing_center[:,1], yerr=exp_packing_center[:,2], 
    lw=0, c='#23f890', elinewidth=1,
    capsize=3, capthick=2, 
    label='$\displaystyle r_{p} = 0.1$, dma-center random phase',
    alpha=1.0)


#ax.plot(dma[:,0], dma[:,1], c='#003366', lw=2, alpha=0.4 )

# dmo-as2
#-----------------------------

ax.plot(patch_pos,packing_jm, lw=2, c='#cd5c5c',
    label='sticky limit, rhombic-pore phases',
    alpha=0.9)

#-------------------------------

ax.set_xlim([0,1])
ax.set_ylim([0.4,1.03])

ax.set_xticks([0.2,0.3,0.4,0.5,0.6,0.7,0.8])
ax.legend(loc='upper center', fancybox=True, prop={'size': 7}, bbox_to_anchor=(0.5, 1.2), ncol=2)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xlabel("$ \displaystyle \Delta$", size=22)
ax.set_ylabel("$ \displaystyle \phi$", size=22)

plt.tight_layout()
plt.savefig("packing_fraction.png", dpi=300)
plt.show()
#------------------------------------------
