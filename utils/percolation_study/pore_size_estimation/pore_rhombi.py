import numpy as np
import pandas as pd
import networkx as nx
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
        fid = dir_i

        yield (fid, ptype, phi, temperature, delta, last_time, pos, orient, box)


def get_edge_points(p, axes, sign_p):
    vertex_n = np.zeros(2)
    vertex_n = p + sign_p[0]*axes[:, 0]/2. + sign_p[1]*axes[:, 1]/2.
    return vertex_n


def get_vertices(p, axes):
    vertices = np.zeros((4, 2))
    vertices[0] = get_edge_points(p, axes, np.array([-1, -1]))
    vertices[1] = get_edge_points(p, axes, np.array([+1, -1]))
    vertices[2] = get_edge_points(p, axes, np.array([+1, +1]))
    vertices[3] = get_edge_points(p, axes, np.array([-1, +1]))

    return vertices


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


def draw_pos(voxcel_pos, blx, frame_name, max_id, min_id, max_domain_id, voxcels):

    scale = 100
    img = np.full(
        (ceil(blx * scale), ceil(blx * scale), 3), 255, np.uint8)

    axes = np.array([[1, 0], [0, 1]])*voxcels.lx

    for di in range(int(min_id), int(max_id)+1):
        if di != max_domain_id:
            points = voxcel_pos[voxcel_pos[:, 0] == di][:, 1:3]

            color = (random.randint(0, 255),
                     random.randint(0, 255),
                     random.randint(0, 255))

            for p in points:
                vert_i = get_vertices(p, axes)
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

    while len(old_coords) > 0:
        neigh = [n for n in G.neighbors(next_i)]
        leftover_neigh = [
            elem_1 for elem_1 in neigh for elem_2 in old_coords if elem_1 == elem_2]
        for ni in leftover_neigh:
            nxi = ni[0]
            nyi = ni[1]

            dxi = next_i[0] - ni[0]
            dyi = next_i[1] - ni[1]

            if abs(dxi) > 1:
                nxi = next_i[0] + np.sign(dxi)
            if abs(dyi) > 1:
                nyi = next_i[1] + np.sign(dyi)

            new_ni = (nxi, nyi)
            new_coords.append(new_ni)
            old_coords.remove(ni)
            nx.relabel_nodes(G, {ni: new_ni}, copy=False)

        next_i = new_ni

        if len(leftover_neigh) == 0:
            k = 1
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

    return old_coords, new_coords


def get_stitched_pos(voxcels, box, domain_obs, G, arr):

    axes = np.array([[1, 0], [0, 1]])*voxcels.lx
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
    edge_pos = []
    for di in range(int(min_id), int(max_id)+1):
        if di != max_domain_id:
            coords = voxcel_pos[voxcel_pos[:, 0] == di][:, 1:3]
            for coord_vi in coords:
                vert_i = get_vertices(coord_vi, axes)
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
    fid, ptype, phi, temperature, delta, last_time, pos, orient, box = vals
    blx = box[3]
    ndim = 2
    side_length = 1
    voxcels_per_cell = 7
    threshold_overlap = 0.25
    outfile = "pore_sizes_vs_7.dat"
    frame_name = "voxcels.png"
    N_trial = 100

    # Generate particles object
    rlx = side_length
    Np = len(pos)
    particles = pt.Rhombi(pos, orient, rlx, ndim)

    # Generate box object
    origin = box[:2] - blx/2.
    box = pt.Box(origin, blx, ndim)

    # Initalize cells object
    # Absolute lmin = particles.sigma
    lmin = particles.sigma
    cells = pt.Cells(lmin, box.lx, ndim)

    # Initialize voxels object
    vxl = cells.lx/voxcels_per_cell
    vxn = cells.nx*voxcels_per_cell
    voxcels = pt.Voxcels(vxl, vxn, box.origin, ndim)

    # Initalize cells object
    cells.list_generate(particles.pos, voxcels.coords, voxcels.pos, box.origin)

    # Calculate voxcel state
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

    # Generate links between empty voxcels
    voxcels.get_links()
    frame_name = "{}_{}".format(fid, frame_name)
    draw(particles, voxcels, box, cells, frame_name)
    # RESULT: pore area/ domain sizes
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

    # Get PBC stitched clusters for convex hull to get asymmetry measures:
    # Stitch together cluster over pbcs by adapting voxcel pos and edge pos
    voxcel_pos, edge_pos = get_stitched_pos(voxcels, box, domain_obs, G, arr)

    shifted_frame_name = "{}_voxcels_shifted.png".format(fid)
    draw_pos(voxcel_pos, box.lx, shifted_frame_name,
             max_id, min_id, max_domain_id, voxcels)

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
    df['pore_area'] = pore_area_no_max
    df['percent_explained_variance'] = xlambda
    df['convex_hull_ratio'] = hull_ratio
    df['circumference'] = circumferences

    meta = {}

    fid = fid.split("/")[1]

    meta["fid"] = "{}_{}".format(fid, last_time)
    meta["ptype"] = ptype
    meta["phi"] = ptype
    meta["temperature"] = temperature
    meta["delta"] = delta
    meta["last_time"] = last_time
    meta["run_id"] = 0
    meta["max_domain"] = max_domain
    return meta, df


if __name__ == '__main__':

    # read data either through files system via glob or via db
    parser = argparse.ArgumentParser()
    parser.add_argument('-run_id', type=str)
    parser.add_argument('-ncores', type=int)

    args = parser.parse_args()

    gen_fsys = generator_from_fsys(glob.glob("double*/double*/"))

    N_CORES = int(args.ncores)
    N_CORES_MAX = 8

    if N_CORES > 1 and N_CORES <= N_CORES_MAX:
        print("Multiprocessing with {} cores".format(N_CORES))
        pool = multiprocessing.Pool(N_CORES)
        results = pool.map(calculate, gen_fsys)
        pool.close()
        pool.join()

    if N_CORES == 1:
        print("single core job")
        for vals in gen_fsys:
            results = calculate(vals)

    if N_CORES > N_CORES_MAX:
        print("Too many cores allocated, please do not use more than {} cores".format(
            N_CORES_MAX))
        exit()

    # create hd5f file
    f = h5py.File("pore_measures_{}.h5".format(args.run_id), 'w')
    for i, res in enumerate(results):
        grp = f.create_group(res[0]['fid'])

        for key in res[0]:
            grp.attrs[key] = res[0][key]

        dfi = res[1]
        for col in dfi.columns:
            dset = grp.create_dataset(col, len(dfi), dtype='f')
            dset[...] = dfi[col]

    f.close()
