import numpy as np
import pandas as pd
import networkx as nx
import matplotlib as pyplot
from collections import defaultdict
import cv2 as cv
import pore_tool as pt
from math import ceil
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
import configparser
from os.path import exists
import networkx as nx
from collections import defaultdict
import multiprocessing
import glob
import re
import argparse


def generator_from_fsys(fsys_iterator):

    for dir_i in fsys_iterator:
        config = configparser.ConfigParser()
        config.read('{}para.ini'.format(dir_i))

        N = int(config['System']['Number_of_Particles'])
        phi = float(config['System']['Packing_Fraction'])
        temperature = float(config['System']['Temperature'])
        ptype = config['Rhombus']['rhombus_type']
        delta = config['Rhombus']['patch_delta']
        patch_size = config['Rhombus']['patch_size']

        pos_files = glob.glob('{}positions_*.bin'.format(dir_i))

        # get the last value from the string
        def g(x): return int(re.findall(r'\d+', x)[-1])

        mc_times = list(map(g, pos_files))

        last_time = np.max(mc_times)

        pos_file = "{}positions_{}.bin".format(dir_i, last_time)
        pos = np.fromfile(pos_file)
        pos = np.reshape(pos, (-1, 3))
        pos = pos[:, :2]
        orient_file = "{}orientations_{}.bin".format(dir_i, last_time)
        orient = np.fromfile(orient_file)
        orient = np.reshape(orient, (-1, 5))[:, 4]

        box_file = "{}Box_{}.bin".format(dir_i, last_time)
        box = np.fromfile(box_file)

        yield (ptype, phi, temperature, delta, last_time, pos, orient, box)


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


def calculate(vals):
    ptype, phi, temperature, delta, last_time, pos, orient, box = vals

    print("read data...")

    box_x_center = box[0]
    box_y_center = box[1]
    blx = box[3]
    box_ly = box[4]

    print("initialize variables...")
    ndim = 2
    side_length = 1
    voxcels_per_cell = 7
    threshold_overlap = 0.25
    outfile = "pore_sizes_vs_7.dat"
    frame_name = "voxcels.png"
    N_trial = 100

    # generate particles object
    rlx = side_length
    Np = len(pos)
    particles = pt.Rhombi(pos, orient, rlx, ndim)

    # generate box object
    origin = box[:2] - blx/2.
    box = pt.Box(origin, blx, ndim)

    # initalize cells object
    # absolute lmin = particles.sigma
    lmin = particles.sigma
    cells = pt.Cells(lmin, box.lx, ndim)

    # initialize voxels object
    vxl = cells.lx/voxcels_per_cell
    vxn = cells.nx*voxcels_per_cell
    voxcels = pt.Voxcels(vxl, vxn, box.origin, ndim)

    # initalize cells object
    print("generate cell lists...")
    cells.list_generate(particles.pos, voxcels.coords, voxcels.pos, box.origin)

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
    pore_areas, domain_lengths, domains, G = pt.get_pore_volume(voxcels)

    print("write output")


if __name__ == '__main__':

    # read data either through files system via glob or via db
    parser = argparse.ArgumentParser()
    parser.add_argument('-run_id', type=str)
    parser.add_argument('-ncores', type=int)

    args = parser.parse_args()

    gen_fsys = generator_from_fsys(glob.glob("double*/double*/"))

    N_CORES = int(args.ncores)
    N_CORES_MAX = 12

    if N_CORES > 1 and N_CORES <= N_CORES_MAX:
        print("Multiprocessing with {} cores".format(N_CORES))
        pool = multiprocessing.Pool(N_CORES)
        new_results = pool.map(calculate, gen_fsys)
        pool.close()
        pool.join()
        df = pd.concat(new_results)

    if N_CORES == 1:
        print("single core job")
        for vals in gen_fsys:
            new_results = calculate(vals)
            df = df.append(new_results, ignore_index=True)

    if N_CORES > N_CORES_MAX:
        print("Too many cores allocated, please do not use more than {} cores".format(
            N_CORES_MAX))
