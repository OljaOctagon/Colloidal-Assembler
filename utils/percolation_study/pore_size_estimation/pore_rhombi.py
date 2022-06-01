import numpy as np
import pandas as pd
import networkx as nx
import matplotlib as pyplot
from collections import defaultdict
import cv2 as cv
import imutils
import pore_tool as pt
from math import ceil
import yaml
from yaml.loader import SafeLoader
import struct 

# TODO find better output format
# TODO Think of how to run for all systems and output

def draw(particles, voxcels, box, cells, frame_name):
    scale = 100
    img = np.full(
        (ceil(box.lx * scale), ceil(box.ly * scale), 3), 255, np.uint8)
    for vert_i in particles.vertices:
        cv.fillPoly(img, np.int32([vert_i * scale]), (0, 0, 0))

    for coord_vi in voxcels.coords:
        if voxcels.fill_state[coord_vi] == 0:
            vert_i = voxcels.get_vertices(coord_vi)
            cv.rectangle(img, np.int32(
                vert_i[2] * scale), np.int32(vert_i[0] * scale), (255, 0, 0), 2)

    outsize = (10000, 10000)
    out = cv.resize(img, outsize)
    cv.imwrite(frame_name, out)

def write_voxcel_edge_points(domains,voxcels):
    arr=[]
    for di, domain in enumerate(domains):
        for coord_vi in domain:
            vert_i = voxcels.get_vertices(coord_vi)
            for vi in vert_i:
                arr.append([di, vi[0],vi[1]])
               
    arr=np.array(arr)
    arr.tofile("edge_points_voxcels.bin")

def write_voxcel_pos(domains, voxcels):
    arr=[]
    for di, domain in enumerate(domains):
        for coord_vi in domain:
            pos_i = voxcels.pos[coord_vi]
            arr.append([di, pos_i[0],pos_i[1]])

    arr=np.array(arr)
    arr.tofile("pos_voxcels.bin")

def write_voxcel_coords(domains,voxcels):
    arr=[]
    for di, domain in enumerate(domains):
        for coord_vi in domain:
            arr.append([di, coord_vi[0],coord_vi[1]])

    arr=np.array(arr)
    arr.tofile("coord_voxcels.bin")


if __name__ == '__main__':

    with open('param_pore_rhombi.yaml') as f:
        param = yaml.load(f, Loader=SafeLoader)


    print("read data...")

    ptype=param['ptype']
    fdir = "example_data/rhombi/{}".format(ptype)
    pos_file = "positions.bin"
    orient_file = "orientations.bin"
    box_file = "Box.bin"

    box_all = np.fromfile("{}/Box.bin".format(fdir))
    box_x_center = box_all[0]
    box_y_center = box_all[1]
    blx = box_all[3]
    box_ly = box_all[4]

    pos = np.fromfile("{}/positions.bin".format(fdir))
    pos = np.reshape(pos, (-1, 3))
    pos = pos[:, :2]

    orient = np.fromfile("{}/orientations.bin".format(fdir))
    orient = np.reshape(orient, (-1, 5))[:, 4]

    print("initialize variables...")
    

    ndim = int(param['dimension'])

    # generate particles object
    rlx = float(param['side_length'])
    Np = len(pos)
    particles = pt.Rhombi(pos, orient, rlx, ndim)

    # generate box object
    origin = box_all[:2] - blx/2.
    box = pt.Box(origin, blx, ndim)

    # initalize cells object
    # absolute lmin = particles.sigma
    lmin = particles.sigma
    cells = pt.Cells(lmin, box.lx, ndim)

    # initialize voxels object
    voxcels_per_cell = int(param['voxcels_per_cell'])
    vxl = cells.lx/voxcels_per_cell
    vxn = cells.nx*voxcels_per_cell
    voxcels = pt.Voxcels(vxl, vxn, box.origin, ndim)

    # initalize cells object
    print("generate cell lists...")
    cells.list_generate(particles.pos, voxcels.coords, voxcels.pos, box.origin)

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
    pore_volumes, domain_lengths, domains = pt.get_pore_volume(voxcels)

    output_frame = param['output_frame']
    if output_frame:
        frame_name = param['frame_name']
        draw(particles, voxcels, box, cells, frame_name)

    outfile = param['outfile']
    with open(outfile, 'w') as f:
        for v, d in zip(pore_volumes, domain_lengths):
            f.write("{},{}\n".format(v, d))

    param_out = "parameters.yaml"
    with open(param_out, 'w') as f:
        f.write("ptype: {}\n".format(ptype))
        f.write("cell_lx: {}\n".format(cells.lx))
        f.write("cell_nx: {}\n".format(cells.nx))
        f.write("voxcel_lx: {}\n".format(voxcels.lx)) 
        f.write("voxcel_nx: {}\n".format(voxcels.nx))
        f.write("N_particles: {}\n".format(Np))

    write_voxcel_edge_points(domains,voxcels)
    write_voxcel_pos(domains,voxcels)
    write_voxcel_coords(domains,voxcels)

