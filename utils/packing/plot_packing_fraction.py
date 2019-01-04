import numpy as np 
import pandas as pd 

def area_particles(Lx):
    '''The unit cell is given by  full particles and 6 half ones.
    That results in 9 full 60 degree rhombi. '''
    area_rhombi = Lx*np.sin(np.pi/3.) 
    return 3*area_rhombi

def area_pores(Lx, patch_pos):
    ''' The pores are hexagon and triangles of 
    l=(fabs(patch_pos - (1-patch_pos)). In total there is 
    one full hexagon, 6 1/3 hexagon, and 6 triangles. 
    That equals the area of 4 full hexagons '''
    pore_p = 2*patch_pos - 1
    pore_length = pore_p*Lx
    hexagon = (3*np.sqrt(3)/2)*np.power(pore_length,2)
    triangle = hexagon*(1/6.)
    return hexagon + 2*triangle


def area_particles_jm(Lx):
    ''' the particles of the unit cell are one rhombi'''
    area_rhombi = np.power(Lx,2)*np.sin(np.pi/3.) 
    return area_rhombi

def area_pores_jm(Lx,patch_pos):
    ''' the unit cell has one rhombic pore. 
    Its area is a function of the patch size. '''
    pore_length = np.fabs(Lx - 2*patch_pos)
    rhombi_p = np.power(pore_length,2)*np.sin(np.pi/3.)
    return rhombi_p

# theoretical dma-as2
#--------------------

Lx = 1
area_pa = area_particles(Lx)
patch_pos = np.arange(0,1.01,0.01)
area_po = area_pores(Lx, patch_pos)
total_area = area_pa + area_po
packing = area_pa/total_area


# theoretical dmo-as2 

#----------------------------

Lx = 1
area_pa_jm = area_particles_jm(Lx)
patch_pos_jm = np.arange(0,1.01,0.01)
area_po_jm = area_pores_jm(Lx, patch_pos)
total_area_jm = area_pa_jm + area_po_jm
packing_jm = area_pa_jm/total_area_jm

#--------------------------------------------------------

import matplotlib.pyplot as plt 
from matplotlib import rc 
import seaborn as sns 
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
sns.set(style='whitegrid')

fig, ax = plt.subplots()

# dma-as2
#-----------

exp_file = "measured_packing_dma-as2.dat"

exp_packing = pd.read_csv(exp_file, header=None,
    delim_whitespace=True).values

print(exp_packing)
exp_packing_center = exp_packing[exp_packing[:,0] == 0.5]
exp_packing_off_center = exp_packing[exp_packing[:,0]!=0.5]


ax.plot(patch_pos,packing, lw=2, c='#336699', 
    label='sticky limit, triangular-pore phase',
    alpha=0.9)

ax.plot(exp_packing_off_center[:,0],
        exp_packing_off_center[:,1],
        lw=0,
        marker='s',
        c='#003366',
        label='$\displaystyle r_{p} = 0.05$, triangular-pore phase',
        alpha=1.0)

ax.plot(exp_packing_center[:,0],
        exp_packing_center[:,1],
        lw=0,
        c='#008365',
        marker='o',
        label='$\displaystyle r_{p} = 0.05$, dma-c random phase',
        alpha=1.0)


# dmo-as2
#-----------------------------

exp_file_2 = "measured_packing_dmo-as2.dat"

exp_packing_2 = pd.read_csv(exp_file, header=None,
    delim_whitespace=True).values

exp_packing_center_2 = exp_packing_2[exp_packing_2[:,0] == 0.5]
exp_packing_off_center_2 = exp_packing_2[exp_packing_2[:,0]!=0.5]

ax.plot(patch_pos,packing_jm, lw=2, c='#cd5c5c',
    label='sticky limit, rhombic-pore phases',
    alpha=0.9)


ax.plot(exp_packing_off_center_2[:,0],
            exp_packing_off_center_2[:,1],
            lw=0, c='r',
            marker='d',
            label='$\displaystyle r_{p} = 0.05$, rhombic-pore phase',
            alpha=1.0)


ax.plot(exp_packing_center_2[:,0],
            exp_packing_center_2[:,1],
            lw=0, c='#e77a0b',
            marker='*',
            label='$\displaystyle r_{p} = 0.05$, dmo-c random phase',
            alpha=1.0)

#-------------------------------

ax.set_xlim([0,1])
ax.set_ylim([0.4,1.03])

ax.set_xticks([0.2,0.3,0.4,0.5,0.6,0.7,0.8])
ax.set_yticks([0.4,0.5,0.6,0.7,0.8,0.9,1.0])
ax.legend(loc='upper center', fancybox=True, prop={'size':10}, bbox_to_anchor=(0.5, 1.4), ncol=2)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xlabel("$ \displaystyle \Delta$", size=22)
ax.set_ylabel("$ \displaystyle \phi$", size=22)

plt.tight_layout()
plt.savefig("packing_fraction.png", dpi=300)
plt.show()
#------------------------------------------
