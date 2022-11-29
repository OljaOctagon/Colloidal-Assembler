import numpy as np
import cv2 as cv
import glob
from math import ceil
import pore_tool as pt
import multiprocessing
import argparse
import configparser


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
        fid = dir_i

        yield (fid, ptype, phi, temperature, delta, last_time, pos, orient, box)


def draw(particles, box, frame_name):
    scale = 100
    img = np.full(
        (ceil(box.lx * scale), ceil(box.ly * scale), 3), 255, np.uint8)
    for vert_i in particles.vertices:
        cv.fillPoly(img, np.int32([vert_i * scale]), (153, 153, 255))

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

    frame_name = frame_name = "{}_phi_{}_delta_{}_temp_{}_rhombi".format(
        ptype, phi, delta, temperature)

    draw(particles, box, frame_name)


if __name__ == '__main__':

    # read data either through files system via glob or via db
    parser = argparse.ArgumentParser()
    parser.add_argument('-ncores', type=int)

    args = parser.parse_args()

    gen_fsys = generator_from_fsys(glob.glob("double*/double*/"))

    N_CORES = int(args.ncores)
    N_CORES_MAX = 8

    if N_CORES > 1 and N_CORES <= N_CORES_MAX:
        print("Multiprocessing with {} cores".format(N_CORES))
        pool = multiprocessing.Pool(N_CORES)
        results = pool.map(get_frames, gen_fsys)
        pool.close()
        pool.join()

    if N_CORES == 1:
        print("single core job")
        for vals in gen_fsys:
            results = get_frame(vals)

    if N_CORES > N_CORES_MAX:
        print("Too many cores allocated, please do not use more than {} cores".format(
            N_CORES_MAX))
        exit()
