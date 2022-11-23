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
import random
import h5py


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


def draw_pos(voxcel_pos, blx, frame_name, max_id, min_id, max_domain_id, voxcels):

    scale = 100
    img = np.full(
        (ceil(blx * scale), ceil(blx * scale), 3), 255, np.uint8)

    for di in range(int(min_id), int(max_id)+1):
        if di != max_domain_id:
            points = voxcel_pos[voxcel_pos[:, 0] == di][:, 1:3]

            color = (random.randint(0, 255),
                     random.randint(0, 255),
                     random.randint(0, 255))

            for p in points:
                vert_i = voxcels.get_vertices(p)
                cv.rectangle(img, np.int32(
                    vert_i[2] * scale), np.int32(vert_i[0] * scale), color, 2)

    outsize = (10000, 10000)
    out = cv.resize(img, outsize)
    cv.imwrite(frame_name, out)


def get_voxel_array(domains, voxcels):
    arr = []
    for di, domain in enumerate(domains):
        for coord_vi in domain:
            arr.append([di, coord_vi[0], coord_vi[1]])

    arr = np.array(arr)

    arr = np.array(arr)
    arr = np.reshape(arr, (-1, 3))
    return arr


def stitch_cluster(G, next_i, old_coords, new_coords):
    # print(list(G.nodes))

    while len(old_coords) > 0:
        #print("next_i", next_i)
        neigh = [n for n in G.neighbors(next_i)]
        leftover_neigh = [
            elem_1 for elem_1 in neigh for elem_2 in old_coords if elem_1 == elem_2]
        for ni in leftover_neigh:
            nxi = ni[0]
            nyi = ni[1]

            dxi = next_i[0] - ni[0]
            dyi = next_i[1] - ni[1]

            if np.abs(dxi) > 1:
                nxi = next_i[0] - np.sign(nxi)
            if np.abs(dyi) > 1:
                nyi = next_i[1] - np.sign(nyi)

            new_ni = (nxi, nyi)
            new_coords.append(new_ni)
            old_coords.remove(ni)
            #print("ni, new_ni", ni, new_ni)
            nx.relabel_nodes(G, {ni: new_ni}, copy=False)

        next_i = new_ni

        if len(leftover_neigh) == 0:
            #print("empty step")
            k = 1
            # print(old_coords)
            while len(leftover_neigh) == 0:
                last_i = new_coords[new_coords.index(next_i)-k]
                neigh = [n for n in G.neighbors(last_i)]
                i = 0

                while len(leftover_neigh) == 0 and i < len(neigh):
                    next_test = neigh[i]
                    neighk = [n for n in G.neighbors(next_test)]
                    leftover_neigh = [
                        elem_1 for elem_1 in neighk for elem_2 in old_coords if elem_1 == elem_2]
                    i += 1

                k += 1
            next_i = next_test
            #print("reduce, next_i", next_i)

    return old_coords, new_coords


def get_stitched_pos(voxcels, box, domain_obs, G, arr):

    min_id, max_id, max_domain_id = domain_obs
    voxcel_pos = []
    # random.seed(10)
    for di in range(int(min_id), int(max_id)+1):
        if di != max_domain_id:
            coords = arr[arr[:, 0] == di][:, 1:3]
            old_coords = [(i, j) for i, j in coords]
            rand_c = coords[random.randint(0, len(coords)-1)]
            ci = (rand_c[0], rand_c[1])
            new_coords = [ci]
            old_coords.remove(ci)

            old_coords, new_coords = stitch_cluster(
                G, ci, old_coords, new_coords)

            for entry in new_coords:
                nni = entry
                vx = box.origin[0] + voxcels.lx*nni[0] + voxcels.lx/2
                vy = box.origin[1] + voxcels.lx*nni[1] + voxcels.lx/2
                voxcel_pos.append([di, vx, vy])

    voxcel_pos = np.array(voxcel_pos)
    draw_pos(voxcel_pos, box.lx, "voxcels_shifted.png",
             max_id, min_id, max_domain_id, voxcels)

    edge_pos = []
    for di in range(int(min_id), int(max_id)+1):
        if di != max_domain_id:
            coords = voxcel_pos[voxcel_pos[:, 0] == di][:, 1:3]
            for coord_vi in coords:
                vert_i = voxcels.get_vertices(coord_vi)
                for vi in vert_i:
                    edge_pos.append([di, vi[0], vi[1]])

    edge_pos = np.array(edge_pos)

    return voxcel_pos, edge_pos


def get_ALL_circumferences(G, domains, n_edges, domain_obs, vlx):

    min_id, max_id, max_domain_id = domain_obs

    def ete_degree(G):
        degree = defaultdict(int)
        for u, v in ((u, v) for u, v, d in G.edges(data=True) if d['ete'] == 1):
            degree[u] += 1
            degree[v] += 1
        return degree

    degree_dict = ete_degree(G)

    def get_ete_degree(u):
        return degree_dict[u]

    domain_dg = []
    for di in range(int(min_id), int(max_id)+1):
        if di != max_domain_id:
            domain_dg.append(list(map(get_ete_degree, domains[di])))

    def get_circumference(domain_dgi):
        arr = (n_edges - np.array(domain_dgi))*vlx
        cf = np.sum(arr)
        return cf

    circumferences = [get_circumference(dd) for dd in domain_dg]
    return circumferences


def calculate(vals):
    ptype, phi, temperature, delta, last_time, pos, orient, box = vals

    print("read data...")
    blx = box[3]

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

    # RESULT: pore area/ domain sizes
    print("get all pore volumes and domain lengths ")
    pore_areas, domain_lengths, domains, G = pt.get_pore_volume(voxcels)

    arr = get_voxel_array(domains, voxcels)

    min_id = np.min(arr[:, 0])
    max_id = np.max(arr[:, 0])
    uniques, counts = np.unique(arr[:, 0], return_counts=True)
    max_domain_id = uniques[np.argmax(counts)]
    max_domain = np.max(pore_areas)

    domain_obs = (min_id, max_id, max_domain_id)

    # RESULT: pore circumference
    circumferences = get_ALL_circumferences(G, domains,
                                            particles.N_edges,
                                            domain_obs,
                                            voxcels.lx)

    # get PBC stitched clusters for convex hull to get asymmetry measures:
    # stitch together cluster over pbcs by adapting voxcel pos and edge pos
    voxcel_pos, edge_pos = get_stitched_pos(voxcels, box, domain_obs, G, arr)

    # RESULT: ratio pore volume and convex hull: general asymmetry measure
    hull = []
    pore_area_no_max = []
    for di in range(int(min_id), int(max_id)+1):
        if di != max_domain_id:
            points = edge_pos[edge_pos[:, 0] == di][:, 1:3]
            pore_area = (len(points)/particles.N_edges)*np.power(voxcels.lx, 2)
            pore_area_no_max.append(pore_area)
            chull = ConvexHull(points)
            hull.append([pore_area, chull.volume])

    hull = np.array(hull)
    hull_ratio = hull[:, 0]/hull[:, 1]

    # RESULT: explained variance by largest principal component:
    # how elongated are pores
    pca = PCA(n_components=2)
    xlambda = []
    for di in range(int(min_id), int(max_id)+1):
        if di != max_domain_id:
            points = edge_pos[edge_pos[:, 0] == di][:, 1:3]
            pca.fit(points)
            xlambda.append(pca.explained_variance_ratio_[0])

    df = pd.DataFrame()
    df['pore area'] = pore_area_no_max
    df['percent_explained_variance'] = xlambda
    df['convex_hull_ratio'] = hull_ratio
    df['circumference'] = circumferences

    meta = (ptype, phi, temperature, delta, last_time, max_domain)

    return meta, df


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
        meta_data, results = pool.map(calculate, gen_fsys)
        pool.close()
        pool.join()
        #df = pd.concat(new_results)

    if N_CORES == 1:
        print("single core job")
        for vals in gen_fsys:
            meta_data, results = calculate(vals)
            print(meta_data, results)
            #df = df.append(new_results, ignore_index=True)

    if N_CORES > N_CORES_MAX:
        print("Too many cores allocated, please do not use more than {} cores".format(
            N_CORES_MAX))
