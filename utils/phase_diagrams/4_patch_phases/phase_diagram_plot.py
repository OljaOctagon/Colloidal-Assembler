import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib import colors 
import json
import itertools
import matplotlib.patches as mpatches
import matplotlib.transforms as transforms
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

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

def size_zero(Lx, patch_pos):
    return 0

def list_append(lst, item):
    lst.append(item)
    return lst

#TODO remove when parse_json finished
def parse_json_old(filen):
    with open(filen) as fhandle:
        data = json.load(fhandle)

    data = np.array([[1-float(key), data[key][1]] for key in data ])
    data = data[data[:,0].argsort()]
    return data

def parse_json(filen):
    with open(filen) as fhandle:
        data=json.load(fhandle)

    arr = np.array([ data[key1][key2][0] for key1 in data for key2 in data[key1]])
    patch_values=[0.5,0.6,0.7,0.8]
    print(filen, list(data.keys()))
    energy_values = list(data['0.2'].keys())
    arr = np.reshape(arr, (len(patch_values), len(energy_values)))
    arr = arr[:,::-1]
    arr =np.concatenate((arr, arr[::-1][1:]))

    return patch_values, energy_values, arr 

def draw_hexagon(grid_point, pore_size):
    return mpatches.RegularPolygon(grid_point, 6, pore_size,
                                   edgecolor='k',
                                   facecolor='none')

def draw_rhombi(grid_point, pore_size):
    trans = transforms.Affine2D().skew(0,15).rotate(-15).translate(grid_point[0], grid_point[1]) + ax.transData

    return mpatches.Rectangle([0,0], pore_size, pore_size,
                              edgecolor='k',
                              facecolor='none',
                              transform=trans)

def draw_zero(grid_point, pore_size):
    trans = transforms.Affine2D().skew(0,15).rotate(-15).translate(grid_point[0], grid_point[1]) + ax.transData

    return mpatches.Rectangle([0,0], pore_size, pore_size,
                              edgecolor='k',
                              facecolor='none',
                              transform=trans)

def draw_pores(patch_values, energy_values,pore_size_function, poly_function):
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
            ax.add_patch(polygon)

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


#-----------------------------------------------------

if __name__ == "__main__":

    files=[
        'psi_mean_double_manta_Asymm_2.json',
        'psi_mean_double_mouse_Asymm_2.json',
        'psi_mean_checkers_Asymm_2.json',
        'psi_mean_double_manta_symm_1.json',
        'psi_mean_checkers_symm_1.json',
        'psi_mean_checkers_Asymm_1.json']
    name_dict = dict(zip(files, [
        'dma_as2',
        'dmo_as2',
        'checkers_as2',
        'dma_s1',
        'checkers_s1'
        'checkers_as1']))

    poly_functions = [
        draw_hexagon,
        draw_rhombi,
        draw_rhombi,
        draw_zero,
        draw_zero,
        draw_zero]
    
    poly_function_dict = dict(zip(files, poly_functions))
    pore_size_functions = [
        size_hexagon,
        size_rhombi,
        size_rhombi,
        size_zero,
        size_zero,
        size_zero]

    pore_dict = dict(zip(files, pore_size_functions))

    for filen in files:

        patch_values, energy_values, zvals = parse_json(filen)
        l_p = len(patch_values)
        l_e = len(energy_values)

        fig, ax = plt.subplots()
        img = plt.imshow(zvals,interpolation='nearest',
                        cmap = cmap, norm = norm)

        pos_start = 0
        x_og = np.arange(pos_start,pos_start+l_p,1)
        y_og = np.arange(pos_start,pos_start+l_e,1)
        ax.set_xticks(x_og)
        ax.set_yticks(y_og)

        plt.xticks(x_og, patch_values, rotation='horizontal')
        plt.yticks(y_og, energy_values, rotation='horizontal')

        # Minor ticks
        ax.set_xticks(x_og-0.5, minor=True);
        ax.set_yticks(y_og-0.5, minor=True);

        plt.grid(b=True,
                which='minor', 
                linestyle='-', 
                color='k',
                linewidth=2)

        draw_pores(patch_values,
                   energy_values,
                   pore_dict[filen],
                   poly_function_dict[filen])

        plt.xlabel("$\\Delta$", size=20)
        plt.ylabel("$\\epsilon [ k_{B}T ]$", size=20)
        #cbar = plt.colorbar(ticks=[-1.1,-1,0,1])
        #cbar.ax.set_yticklabels(['liquid', '-1 (np)','0 (r) $\psi$', '1 (p)'])
        cbar = plt.colorbar(boundaries=np.linspace(-1,1,20), ticks=[-1,0,1])
        cbar.ax.set_yticklabels(['-1 (np)', '0 (r)', '1 (p)'])
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label('$\psi$', size=24)
        rect1 = mpatches.Rectangle((5,4.3),0.35,0.35,
                                   clip_on=False,
                                   edgecolor='k',
                                   linewidth=1,
                                   facecolor='#ffff00')
        plt.text(5.5,4.6,"liquid", size=13)

        rect2 = mpatches.Rectangle((5,4.8),0.35,0.35,
                                   clip_on=False,
                                   edgecolor='k',
                                   linewidth=1,
                                   facecolor='#bdbdbd')

        plt.text(5.5,5.0,"clusters", size=13)

        # Add the patch to the Axes
        ax.add_patch(rect1)
        ax.add_patch(rect2)

        plt.tick_params(axis='both',labelsize=12)
        plt.tight_layout()
        plt.savefig("phase_diagram_{}.png".format(name_dict[filen]), dpi=300)

    plt.show()
