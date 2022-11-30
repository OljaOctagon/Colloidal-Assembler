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


def draw(particles, box):
    scale = 100
    img = np.full(
        (ceil(box.lx * scale), ceil(box.ly * scale), 3), 255, np.uint8)
    for vert_i in particles.vertices:
        cv.fillPoly(img, np.int32([vert_i * scale]), (255, 153, 153))

    return img


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

    img = draw(particles, box)

    meta = {}
    meta["fid"] = "{}_{}".format(fid, last_time)
    meta["ptype"] = ptype
    meta["phi"] = ptype
    meta["temperature"] = temperature
    meta["delta"] = delta
    meta["last_time"] = last_time
    meta["run_id"] = 0

    return meta, img


if __name__ == '__main__':

    # read data either through files system via glob or via db
    parser = argparse.ArgumentParser()
    parser.add_argument('-ncores', type=int)

    args = parser.parse_args()

    gen_fsys = generator_from_fsys()
    gframes = get_frames

    N_CORES = int(args.ncores)
    N_CORES_MAX = 8

    if N_CORES > 1 and N_CORES <= N_CORES_MAX:
        print("Multiprocessing with {} cores".format(N_CORES))
        pool = multiprocessing.Pool(N_CORES)
        img_results = pool.map(gframes, gen_fsys)
        pool.close()
        pool.join()

    if N_CORES == 1:
        print("single core job")
        for vals in gen_fsys:
            img_results = gframes(vals)

    if N_CORES > N_CORES_MAX:
        print("Too many cores allocated, please do not use more than {} cores".format(
            N_CORES_MAX))
        exit()

    N_delta = 3
    delta = [0.2, 0.3, 0.4]

    #temperature = [0.16,0.14,0.12,0.10,0.09,0.07,0.05,0.04,0.01]
    #phi = [0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]

    temperature = [0.10, 0.05]
    phi = [0.1, 0.2]

    N_temperature = len(temperature)
    N_phi = len(phi)

    from collections import defaultdict
    img_dict = defaultdict(np.empty((N_phi, N_temperature)))
    for res in img_results:
        tag = "{}_{}".format(res[0]['ptype'], res[0]['delta'])
        phi_i = res[0]['phi']
        temp_i = res[0][temperature]

        if phi_i in phi and temp_i in temperature:
            index_i = phi.index(phi_i)
            index_j = temperature.index(temp_i)

            img_dict[tag][index_i][index_j] = res[1]

    for key in img_dict:
        frame_name = "panel_{}.png".format(tag)
        outsize = (10000, 10000)

        out = cv.resize(img_dict[key], outsize)
        cv.imwrite(frame_name, out)
