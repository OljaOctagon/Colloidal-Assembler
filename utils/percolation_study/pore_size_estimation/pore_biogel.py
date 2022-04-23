import numpy as np 
import pandas as pd 
import networkx as nx 
import matplotlib as pyplot 
from collections import defaultdict 
import pore_tool as pt 
from math import ceil
import yaml
from yaml.loader import SafeLoader


if __name__ == '__main__':

	print("read data...")
	fdir = "example_data/biogel/positions.xyz" 
	df_pos = pd.read_csv("{}/{}".format(fdir,fdata), delim_whitespace=True, names=["type,",'x','y','z'])
	pos = df_pos[['x','y','z']].values 

	print("initialize variables...")
    with open('param_pore_biogel.yaml') as f:
        param = yaml.load(f, Loader=SafeLoader)

    ndim = int(param['dimension'])

	# generate particles object 
	sigma= float(param['sigma'])
	particles = pt.Spheres(pos,ndim,sigma)

	# generate box object 
	lx = 30 
	origin = (-15)*np.ones(ndim)
	box = pt.Box(origin,lx,ndim)

	# initalize cells object
	# absolute lmin = particles.sigma 
	lmin = particles.sigma 
	cells = pt.Cells(lmin,box.lx,ndim)

	# initialize voxels object
	voxcels_per_cell = int(param['voxcels_per_cell'])
	vxl = cells.lx/voxcels_per_cell
	vxn = cells.nx*voxcels_per_cell 
	voxcels = pt.Voxcels(vxl,vxn, box.origin,ndim)
	
	print("generate cell lists...")
	# initalize cells object 
	cells.list_generate(particles.pos,voxcels.coords, voxcels.pos, box.origin)

	print("calculate voxcel state...")
	
    N_trial = int(param['N_trial'])
    threshold_overlap = float(param['threshold_overlap'])

    print("calculate voxcel state...")
    for ci in cells.coords:

        vcoords = cells.voxcel_list[ci]
        pcoords = []
        neighbours = cells.get_neighbour_coords(ci)
        for nci in neighbours:
            pcoords.extend(cells.particle_list[nci])

        if pcoords:
            pos_v = np.array([voxcels.pos[i] for i in vcoords])
            pos_c = particles.pos[pcoords]
            dist, ndist = pt.get_vdistance(pos_v, pos_c, box.lx)

            for i, j in np.argwhere(ndist < (particles.outer_radius+voxcels.outer_radius)):
                overlap_volume_i = particles.estimate_volume(
                    voxcels, vcoords[i], pcoords[j], box.lx, N_trial)
                if overlap_volume_i > threshold_overlap:
                    voxcels.set_to_filled(vcoords[i])

    print("generate links between empty voxcels...")
    voxcels.get_links()

    print("get all pore volumes and domain lengths ")
    pore_volumes, domain_lengths = pt.get_pore_volume(voxcels)

    outfile = param['outfile']
    with open(outfile, 'w') as f:
        for v, d in zip(pore_volumes, domain_lengths):
            f.write("{},{}\n".format(v, d))
	