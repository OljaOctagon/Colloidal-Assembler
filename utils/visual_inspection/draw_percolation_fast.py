import numpy as np
import cv2 as cv
import glob
from math import ceil
import pore_tool as pt
import multiprocessing
import argparse
import configparser
import re


def generator_from_fsys():

    fsys_iterator = glob.glob("double*/double*/")
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
        fid = dir_i

        yield (fid, ptype, phi, temperature, delta, last_time, pos, orient, box)


def generator_from_single_dir():

    pos_iterator = glob.glob('positions_*.bin')
    # get the last value from the string
    def g(x): return int(re.findall(r'\d+', x)[-1])

    for posf_i in pos_iterator:
        stime = g(posf_i)
        pos_file = "positions_{}.bin".format(stime)
        pos = np.fromfile(pos_file)
        pos = np.reshape(pos, (-1, 3))
        pos = pos[:, :2]
        orient_file = "{}orientations_{}.bin".format(stime)
        orient = np.fromfile(orient_file)
        orient = np.reshape(orient, (-1, 5))[:, 4]

        box_file = "{}Box_{}.bin".format(stime)
        box = np.fromfile(box_file)

        yield (stime, pos, orient, box)


def draw(particles, box, frame_name):
    scale = 100
    img = np.full(
        (ceil(box.lx * scale), ceil(box.ly * scale), 3), 255, np.uint8)
    for vert_i in particles.vertices:
        cv.fillPoly(img, np.int32([vert_i * scale]), (255, 153, 153))

    outsize = (10000, 10000)
    out = cv.resize(img, outsize)
    cv.imwrite(frame_name, out)


def get_frames(vals):
    fid, ptype, phi, temperature, delta, last_time, pos, orient, box = vals
    blx = box[3]
    ndim = 2
    side_length = 1

    # Generate particles object
    Np = len(pos)
    particles = pt.Rhombi(pos, orient, side_length, ndim)

    # Generate box object
    origin = box[:2] - blx/2.
    box = pt.Box(origin, blx, ndim)

    frame_name = "{}_phi_{}_delta_{}_temp_{}_rhombi.png".format(
        ptype, phi, delta, temperature)

    draw(particles, box, frame_name)


def get_all_frames(vals):
    stime, pos, orient, box = vals
    blx = box[3]
    ndim = 2
    side_length = 1

    # Generate particles object
    Np = len(pos)
    particles = pt.Rhombi(pos, orient, side_length, ndim)

    # Generate box object
    origin = box[:2] - blx/2.
    box = pt.Box(origin, blx, ndim)

    frame_name = "rhombi_{}.png".format(stime)

    draw(particles, box, frame_name)


if __name__ == '__main__':

    # read data either through files system via glob or via db
    parser = argparse.ArgumentParser()
    parser.add_argument('-ncores', type=int)
    parser.add_argument('-it', type=str, choices=('all_dirs', 'all_pos'))

    args = parser.parse_args()

    fdict = {'all_dirs': generator_from_fsys,
             'all_pos': generator_from_single_dir}

    gdict = {'all_dirs': get_frames,
             'all_pos': get_all_frames}

    gen_fsys = fdict[args.it]()
    gframes = gdict[args.it]

    N_CORES = int(args.ncores)
    N_CORES_MAX = 8

    if N_CORES > 1 and N_CORES <= N_CORES_MAX:
        print("Multiprocessing with {} cores".format(N_CORES))
        pool = multiprocessing.Pool(N_CORES)
        #pool.map(get_frames, gen_fsys)
        pool.map(gframes, gen_fsys)
        pool.close()
        pool.join()

    if N_CORES == 1:
        print("single core job")
        for vals in gen_fsys:
            results = gframes(vals)

    if N_CORES > N_CORES_MAX:
        print("Too many cores allocated, please do not use more than {} cores".format(
            N_CORES_MAX))
        exit()
