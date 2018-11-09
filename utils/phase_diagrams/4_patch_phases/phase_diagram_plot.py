import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib import colors 
import json
import itertools

import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.transforms as transforms
def size_hexagon(Lx, patch_pos):
    ''' The bigger pores are hexagon and the smaller triangles of 
    l=(fabs(patch_pos - (1-patch_pos)).'''
    pore_p = np.fabs(2*patch_pos - 1)
    pore_size = pore_p*Lx
    return pore_size

def size_triangle(Lx, patch_pos):
    pore_p = np.fabs(2*patch_pos - 1)
    pore_size = pore_p*Lx
    return pore_size

def size_rhombi(Lx,patch_pos):
    ''' the unit cell has one rhombic pore. 
    Its area is a function of the patch size. '''
    pore_p = np.fabs(2*patch_pos - 1)
    pore_size = pore_p*Lx
    return pore_size

def list_append(lst, item):
    lst.append(item)
    return lst

def parse_json(filen):
    with open(filen) as fhandle:
        data = json.load(fhandle)

    data = np.array([[1-float(key), data[key][1]] for key in data ])
    data = data[data[:,0].argsort()]
    return data

files=['psi_mean_assym_manta.json', 'psi_mean_assym_mouse.json']

patch_values=[0.2,0.3,0.4,0.5,0.6,0.7,0.8]
energy_values = [7.2,6.2,5.2,4.2]

l_p = len(patch_values)
l_e = len(energy_values)


# define the colormap
cmap = plt.cm.rainbow
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
cmaplist[0] = (1,1,0,1)
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
# define the bins and normalize
liquid_val=-1.1
bounds = np.linspace(liquid_val,1,21)
norm = colors.BoundaryNorm(bounds, cmap.N)

pore_size_functions = [size_hexagon, size_rhombi]
pore_dict = dict(zip(files, pore_size_functions))

def draw_hexagon(grid_point, pore_size):
    return mpatches.RegularPolygon(grid_point, 6, pore_size)

def draw_rhombi(grid_point, pore_size):
    trans = transforms.Affine2D().skew_deg(0,15).rotate_deg(-15).translate(-0.5,0)
    return mpatches.Rectangle(grid_point, pore_size, pore_size, transform=trans)

poly_functions = [draw_hexagon, draw_rhombi]
poly_function_dict = dict(zip(files, poly_functions))


name_dict = dict(zip(files, ['dma', 'dmo']))
def draw_pores(patch_values, energy_values,pore_size_function, poly_function):
    patches = []
    l_p = len(patch_values)
    l_e = len(energy_values)
    grid = np.mgrid[0:l_p, 0:l_e].reshape(2, -1).T
    grid_keys = itertools.product(patch_values, energy_values)
    grid_dict = dict(zip(grid_keys, grid))
    for patch_delta in patch_values:
        pore_size = pore_size_function(0.5, patch_delta)
        for ei in energy_values[:-1]:
            grid_point = grid_dict[(patch_delta, ei)]
            polygon = poly_function(grid_point, pore_size)
            patches.append(polygon)

    collection = PatchCollection(patches, edgecolors='k', facecolors='none')
    ax.add_collection(collection)


for filen in files:
    parsed_data = parse_json(filen)[:,1]
    parsed_data_double = np.concatenate((parsed_data[::-1][:-1], parsed_data))
    data_4kbt = np.reshape(liquid_val*np.ones(l_p), (1,-1))
    data_5kbt = np.reshape(parsed_data_double, (1,-1))
    data_6kbt = data_5kbt
    data_7kbt = data_5kbt

    zvals = np.concatenate((data_7kbt, data_6kbt,
                            data_5kbt, data_4kbt), axis=0)

    print(zvals)

    fig, ax = plt.subplots()
    #zvals = np.random.randint(-1,1, size=(l_e,l_p))

    img = plt.imshow(zvals,interpolation='nearest',
                     cmap = cmap, norm = norm)

    pos_start = -0.5
    x_og = np.arange(pos_start,pos_start+l_p,1)
    y_og = np.arange(pos_start,pos_start+l_e,1)
    ax.set_xticks(x_og)
    ax.set_yticks(y_og)

    plt.xticks(x_og, patch_values, rotation='horizontal')
    plt.yticks(y_og, energy_values, rotation='horizontal')

    plt.grid(b=True,
             which='major', 
             axis='both', 
             linestyle='-', 
             color='k', linewidth=2)

    draw_pores(patch_values, energy_values, pore_dict[filen], poly_function_dict[filen])
    plt.xlabel("$\\Delta$", size=20)
    plt.ylabel("$\\epsilon [ k_{B}T ]$", size=20)
    cbar = plt.colorbar(ticks=[-1.1,-1,0,1])
    cbar.ax.set_yticklabels(['liquid', '-1 (np)','0 (r) $\psi$', '1 (p)'])
    cbar.ax.tick_params(labelsize=12)
    plt.tick_params(axis='both',labelsize=12)
    plt.tight_layout()
    plt.savefig("phase_diagram_{}.png".format(name_dict[filen]), dpi=300)

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
