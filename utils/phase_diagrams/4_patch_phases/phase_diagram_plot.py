import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib import colors 

patch_values=[0.2,0.3,0.4,0.5,0.6,0.7,0.8]
energy_values = [7.2,6.2,5.2,4.2]

l_p = len(patch_values)
l_e = len(energy_values)
gridsize=(l_p, l_e)

# load psi values for all phase states

# 4.2 0.2 value
# 4.2 0.3 value
# ....
# 7.2 0.8 value 

fig, ax = plt.subplots()
zvals = np.random.randint(-1,1, size=(l_e,l_p))
cmap = plt.get_cmap('jet')

img = plt.imshow(zvals,interpolation='nearest',
                    cmap = cmap)

x_og = np.arange(-.5,6.5,1)
y_og = np.arange(-.5,3.5,1)
ax.set_xticks(x_og)
ax.set_yticks(y_og)

plt.xticks(x_og, patch_values, rotation='horizontal')
plt.yticks(y_og, energy_values, rotation='horizontal')

plt.grid(b=True,
	which='major', 
	axis='both', 
	linestyle='-', 
	color='k', linewidth=2)

plt.xlabel("$\\Delta$", size=25)
plt.ylabel("$\\epsilon [ k_{B}T ]$", size=25)

plt.tight_layout()
plt.savefig("test_phase_diagram.pdf")
plt.show()

'''
# setup the plot
fig, ax = plt.subplots(1,1, figsize=(6,6))

# define the data
x = np.random.rand(20)
y = np.random.rand(20)
tag = np.random.randint(0,20,20)
tag[10:12] = 0 # make sure there are some 0 values to showup as grey

# define the colormap
cmap = plt.cm.jet
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
cmaplist[0] = (.5,.5,.5,1.0)
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

# define the bins and normalize
bounds = np.linspace(0,20,21)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# make the scatter
scat = ax.scatter(x,y,c=tag,s=np.random.randint(100,500,20),cmap=cmap, norm=norm)

# create a second axes for the colorbar
ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')

ax.set_title('Well defined discrete colors')
ax2.set_ylabel('Very custom cbar [-]', size=12)
'''
