import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib import colors 
import json
import itertools
import matplotlib.patches as mpatches
import matplotlib.transforms as transforms
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
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

def parse_json(filen):
    with open(filen) as fhandle:
        data= json.load(fhandle)

    dsorted = lambda x: np.sort(np.array(list(x)).astype(float))

    arr = np.array([ data[str(key1)][str(key2)][1] for key1 in dsorted(data.keys()) for key2 in dsorted(data[str(key1)])])
    patch_values=[0.2,0.3,0.4,0.5,0.6,0.7,0.8]
    print(filen, list(data.keys()))
    energy_values = dsorted(data['0.2'].keys())
    arr = np.reshape(arr, (len(patch_values)-3, len(energy_values)))
    arr = arr[:,::-1]
    arr =np.concatenate((arr, arr[::-1][1:]))
    arr[:,-1] = -1.1
    arr = np.transpose(arr)
    return patch_values, energy_values, arr 

def draw_hexagon(grid_point, pore_size,ax):
    return mpatches.RegularPolygon(grid_point, 6, pore_size,
                                   edgecolor='k',
                                   facecolor='none')

def draw_rhombi(grid_point, pore_size,ax):
    trans = transforms.Affine2D().skew(0,15).rotate(-15).translate(grid_point[0], grid_point[1]) + ax.transData

    return mpatches.Rectangle([0,0], pore_size, pore_size,
                              edgecolor='k',
                              facecolor='none',
                              transform=trans)

def draw_zero(grid_point, pore_size,ax):
    trans = transforms.Affine2D().skew(0,15).rotate(-15).translate(grid_point[0], grid_point[1]) + ax.transData

    return mpatches.Rectangle([0,0], pore_size, pore_size,
                              edgecolor='k',
                              facecolor='none',
                              transform=trans)


def draw_pores(patch_values, energy_values,pore_size_function, poly_function,ax):
    l_p = len(patch_values)
    l_e = len(energy_values)
    grid = np.mgrid[0:l_p, 0:l_e].reshape(2, -1).T
    grid_keys = itertools.product(patch_values, energy_values)
    grid_dict = dict(zip(grid_keys, grid))
    for patch_delta in patch_values:
        pore_size = pore_size_function(0.5, patch_delta)
        for ei in energy_values[:-1]:
            grid_point = grid_dict[(patch_delta, ei)]
            polygon = poly_function(grid_point, pore_size,ax)
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
        'psi_mean_double_manta_symm_1.json',
        'psi_mean_double_manta_symm_2.json',
        'psi_mean_double_manta_Asymm_1.json',
        'psi_mean_double_manta_Asymm_2.json',
        'psi_mean_double_mouse_symm_1.json',
        'psi_mean_double_mouse_symm_2.json',
        'psi_mean_double_mouse_Asymm_1.json',
        'psi_mean_double_mouse_Asymm_2.json',
        'psi_mean_checkers_symm_1.json',
        'psi_mean_checkers_Asymm_1.json',
        'psi_mean_checkers_Asymm_2.json',
    ]


    dmo_to_replace = [
        'psi_mean_double_mouse_symm_1.json',
        'psi_mean_double_mouse_symm_2.json',
        'psi_mean_double_mouse_Asymm_1.json']

    dmo_replacement= 'psi_mean_double_mouse_Asymm_2.json'

    name_dict = dict(zip(files, [
        'dma-s1',
        'dma-s2',
        'dma-as1',
        'dma-as2',
        'dmo-s1',
        'dmo-s2',
        'dmo-as1',
        'dmo-as2',
        'checkers-s1',
        'checkers-as1',
        'checkers-as2',
    ]))

    psi_keys=[1,
              1,
              0,
              1,
              0,
              0,
              0,
              1,
              1,
              1,
              1]

    poly_functions = [
        draw_zero,
        draw_zero,
        draw_zero,
        draw_hexagon,
        draw_zero,
        draw_zero,
        draw_zero,
        draw_rhombi,
        draw_zero,
        draw_zero,
        draw_rhombi]


    poly_function_dict = dict(zip(files, poly_functions))
    pore_size_functions = [
        size_zero,
        size_zero,
        size_zero,
        size_hexagon,
        size_zero,
        size_zero,
        size_zero,
        size_rhombi,
        size_zero,
        size_zero,
        size_rhombi]


    pore_dict = dict(zip(files, pore_size_functions))


    def draw_basic_figure(l_p,l_e, patch_values, energy_values, zvals, fig,ax):
        img = ax.imshow(zvals,interpolation='nearest',
                        cmap = cmap, norm = norm)

        pos_start = 0
        x_og = np.arange(pos_start,pos_start+l_p,1)
        y_og = np.arange(pos_start,pos_start+l_e,1)

        #ax.set_xticks(x_og, patch_values, rotation='horizontal')
        #ax.set_yticks(y_og, energy_values[::-1], rotation='horizontal')

        ax.set_xticks(x_og)
        ax.set_yticks(y_og)

        ax.set_xticklabels(patch_values, minor=False)
        ax.set_yticklabels(energy_values[::-1], minor=False)

        # Minor ticks
        ax.set_xticks(x_og-0.5, minor=True);
        ax.set_yticks(y_og-0.5, minor=True);

        ax.xaxis.set_major_formatter(plt.NullFormatter())
        ax.yaxis.set_major_formatter(plt.NullFormatter())

        ax.grid(b=True,
                which='minor', 
                linestyle='-', 
                color='k',
                linewidth=1)

        draw_pores(patch_values,
                   energy_values,
                   pore_dict[filen],
                   poly_function_dict[filen],
                   ax)

        ax.set_xlabel("$\\Delta$", labelpad=-3, size=13)
        ax.set_ylabel("$\\epsilon$", labelpad=-3, size=13)

        ax.tick_params(axis='both',labelsize=5)
        return img
    def draw_cluster_format():
        pass

    nrow=3
    ncol=4
    fig, axes = plt.subplots(nrow,ncol)
    from itertools import product

    #fig.delaxes(axes[-1, -1])

    Coord = list(product(np.arange(nrow), np.arange(ncol)))[:-1]
    print(len(Coord))

    _,_,zvals_replacement = parse_json(dmo_replacement)

    l_p = 0
    l_e = 0
    for filen,psi_key,coord in zip(files, psi_keys,Coord):
        ax = axes[coord[0], coord[1]]
        patch_values, energy_values, zvals = parse_json(filen)
        l_p = len(patch_values)
        l_e = len(energy_values)

        if filen in dmo_to_replace:
            zvals[:,3] = zvals_replacement[:,3]

        img = draw_basic_figure(l_p,l_e, patch_values, energy_values, zvals,fig,ax)
        ax.set_title(name_dict[filen], size=10)

        if psi_key == 0:
            rect2 = mpatches.Rectangle((-0.5,-1.9),7.1,1.5,
                                        clip_on=False,
                                        edgecolor='None',
                                        linewidth=0,
                                        facecolor='k',
                                        alpha=0.3)
            ax.add_patch(rect2)

    cbar = fig.colorbar(img, ax=axes, orientation='horizontal', fraction=0.05, boundaries=np.linspace(-1,1,20), ticks=[-1,0,1])
    cbar.ax.set_xticklabels(['-1 (np)', '0 (r)', '1 (p)'])
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('$\psi$', size=12)

    rect1 = mpatches.Rectangle((-0.2,0),0.08,1.,
                                clip_on=False,
                                edgecolor='k',
                                linewidth=1,
                                facecolor='#ffff00')


    
    rect3 = mpatches.Rectangle((1.1,0),0.08,1.,
                                clip_on=False,
                                edgecolor='k',
                                linewidth=1,
                                facecolor='k',
                                alpha=0.3)

    cbar.ax.add_patch(rect1)
    cbar.ax.add_patch(rect3)

    fig.text(0.1,0.06,"liquid", size=12)
    fig.text(0.85,0.06,"no phase", size=12)
    fig.text(0.85,0.03,"off center", size=12)
    ax = axes[2,3]

    pos_start = 0
    deltax=1.0
    x_og = np.arange(pos_start,pos_start+l_p,deltax)
    y_og = np.arange(pos_start,pos_start+l_e,deltax)

    ax.set_xticks(x_og)
    ax.set_yticks(y_og)

    ax.set_xlim([pos_start,pos_start + (l_p-1)*deltax])
    ax.set_ylim([pos_start,pos_start + (l_e-1)*deltax])

    ax.set_xticklabels(patch_values, minor=False)
    ax.set_yticklabels(energy_values*-1, minor=False)

    # Minor ticks
    #ax.set_xticks(x_og-(deltax/2), minor=True);
    #ax.set_yticks(y_og-(deltax/2), minor=True);

    ax.set_xlabel("$\\Delta$", labelpad=0.1, size=10)
    ax.set_ylabel("$\\epsilon [ k_{B}T ]$", labelpad=0.1, size=10)
    ax.tick_params(axis='both',labelsize=7)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position(('data',-0.5))
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    #asp = np.diff(ax.get_xlim())[0] / np.diff(ax.get_ylim())[0]
    ax.set_aspect(0.9)

    plt.savefig("phase_diagram.png", dpi=300)
    plt.show()
