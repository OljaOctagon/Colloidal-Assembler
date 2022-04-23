import numpy as np 
import pandas as pd 
import networkx as nx 
import matplotlib as pyplot 
from collections import defaultdict 

import pore_tool as pt 

if __name__ == '__main__':
	print("read data...")
	fdir = "/Users/ada/Documents/Code_Development_2022/biogel/std_conditions_xyz/00115/poly/xyz-files"
	fdata = "biogel_9900.0.xyz"

	df_pos = pd.read_csv("{}/{}".format(fdir,fdata), delim_whitespace=True, names=["type,",'x','y','z'])
	pos = df_pos[['x','y','z']].values 

	print("initialize variables...")

	ndim=3

	# generate particles object 
	particles = pt.Spheres(pos,ndim)

	# generate box object 
	lx = 30 
	origin = (-15)*np.ones(ndim)
	box = pt.Box(origin,lx,ndim)

	# initalize cells object
	# absolute lmin = particles.sigma 
	lmin = particles.sigma 
	cells = pt.Cells(lmin,box.lx,ndim)

	# initialize voxels object
	voxcels_per_cell = 2
	vxl = cells.lx/voxcels_per_cell
	vxn = cells.nx*voxcels_per_cell 
	voxcels = pt.Voxcels(vxl,vxn, box.origin,ndim)
	
	print("generate cell lists...")
	# initalize cells object 
	cells.list_generate(particles.pos,voxcels.coords, voxcels.pos, box.origin)

	print("calculate voxcel state...")
	print(len(cells.coords))
	for ci in cells.coords:
		
		# Overlap lists 
		vcoords = cells.voxcel_list[ci]
		pcoords = []
		neighbours = cells.get_neighbour_coords(ci)
		for nci in neighbours:
		    pcoords.extend(cells.particle_list[nci])
		
		# get distances between partices and voxcels   
		pos_v = np.array([ voxcels.pos[i] for i in vcoords])
		pos_c = particles.pos[pcoords]
		dist, ndist = get_vdistance(pos_v,pos_c,box.lx)
		
		# overlap for sure 
		for i, j in np.argwhere(ndist<(particles.radius+voxcels.inner_radius)):
		    overlap_volume_i = estimate_volume(voxcels, vcoords[j], particles, pcoords[i], box.lx)
		    if overlap_volume_i>0.4:
		        voxcels.set_to_filled(vcoords[j])
	
	print("generate links between empty voxcels...")
	voxcels.get_links()
	
	pirnt("get all pore volumes and domain lengths ")
	pore_volumes, domain_lengths = get_pore_volume(voxcels.links)
	print(domain_lengths)

	