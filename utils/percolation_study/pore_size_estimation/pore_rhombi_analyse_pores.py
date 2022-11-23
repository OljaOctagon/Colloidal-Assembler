import networkx as nx
import numpy as np
import yaml
from yaml.loader import SafeLoader
from collections import defaultdict
import matplotlib.pyplot as plt
import random
from math import ceil
import cv2 as cv
import imutils
import struct
import sys


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


def draw_pos(voxcel_pos, blx, frame_name, max_id, min_id, max_domain_id, axes):

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
                vert_i = get_vertices(p, axes)
                cv.rectangle(img, np.int32(
                    vert_i[2] * scale), np.int32(vert_i[0] * scale), color, 2)

    outsize = (10000, 10000)
    out = cv.resize(img, outsize)
    cv.imwrite(frame_name, out)


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


if __name__ == '__main__':

    # load/set params
    with open('parameters.yaml') as f:
        param = yaml.load(f, Loader=SafeLoader)

    ndim = 2
    voxcel_lx = float(param['voxcel_lx'])
    voxcel_area = np.power(voxcel_lx, ndim)
    n_edges = 4

    lx = 1
    alpha = np.pi/3.
    sigma = lx*np.sqrt(2+2*np.cos(alpha))
    vs = 7
    ptype = param['ptype']

    box_all = np.fromfile("example_data/rhombi/{}/Box.bin".format(ptype))
    blx = box_all[3]
    total_area = blx*blx
    origin = box_all[:2] - blx/2.

    phi_particles = 0.125
    phi_all_pores = 1 - phi_particles
    area_rhombus = lx*lx * np.sin(alpha)

    axes = np.array([[1, 0], [0, 1]])*voxcel_lx

    # load files
    G = nx.read_gpickle("pore_network.gpickle")

    arr = np.reshape(np.fromfile("coord_voxcels.bin", dtype=int), (-1, 3))
    min_id = np.min(arr[:, 0])
    max_id = np.max(arr[:, 0])
    uniques, counts = np.unique(arr[:, 0], return_counts=True)
    max_domain_id = uniques[np.argmax(counts)]

    # stitch together cluster over pbcs by adapting voxcel pos and edge pos
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
                vx = origin[0] + voxcel_lx*nni[0] + voxcel_lx/2
                vy = origin[1] + voxcel_lx*nni[1] + voxcel_lx/2
                voxcel_pos.append([di, vx, vy])

    voxcel_pos = np.array(voxcel_pos)
    draw_pos(voxcel_pos, blx, "voxcels_shifted.png",
             max_id, min_id, max_domain_id, axes)

    edge_pos = []
    for di in range(int(min_id), int(max_id)+1):
        if di != max_domain_id:
            coords = voxcel_pos[voxcel_pos[:, 0] == di][:, 1:3]
            for coord_vi in coords:
                vert_i = get_vertices(coord_vi, axes)
                for vi in vert_i:
                    edge_pos.append([di, vi[0], vi[1]])

    edge_pos = np.array(edge_pos)

    # 1. get domain lengths and pore area
    domains = list(nx.connected_components(G))
    domain_lengths = np.array([len(domain) for domain in domains])
    pore_area = voxcel_area*domain_lengths

    fig, ax = plt.subplots()
    pore_area_excluded_void = np.sort(pore_area)[:-1]
    hist, bin_edges = np.histogram(pore_area_excluded_void, density=True)
    start = bin_edges[0]
    x = (bin_edges[:-1] + (bin_edges[1]-bin_edges[2])/2)/area_rhombus
    plt.plot(x, hist, lw=1, label='$l_p/l_r = {}$'.format(np.round(sigma/vs, 3)))

    plt.legend()
    plt.xlabel("$A_{p}  / A_{r}$", size=15)
    plt.ylabel("P", size=15)
    plt.xlim([0, 110])
    plt.ylim([0, 0.06])
    plt.savefig("results_pore_area.pdf")

    # 2. pore circumference
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
    for domain in domains:
        domain_dg.append(list(map(get_ete_degree, domain)))

    def get_circumference(domain_dgi):
        arr = (n_edges - np.array(domain_dgi))*voxcel_lx
        cf = np.sum(arr)
        return cf

    circumferences = [get_circumference(dd) for dd in domain_dg]

    fig, ax = plt.subplots()
    pore_cf_excluded_void = np.sort(circumferences)[:-1]
    hist, bin_edges = np.histogram(pore_cf_excluded_void, density=True)
    start = bin_edges[0]
    x = (bin_edges[:-1] + (bin_edges[1]-bin_edges[2])/2)/lx
    plt.plot(x, hist, lw=1, label='$l_p/l_r = {}$'.format(np.round(sigma/vs, 3)))

    plt.legend()
    plt.xlabel("$Cf_{p} / l_{r}$", size=15)
    plt.xlim([0, 90])
    plt.ylim([0, 0.06])

    plt.ylabel("P", size=15)
    plt.savefig("results_pore_cf.pdf")

    # 3. pore assymmetry
    # 3.1 compare convex hull with actuall hull: measure of concavity
    # if A/Ahull = 1 -- convex cluster, if A/Ahull <1 concave
    from scipy.spatial import ConvexHull

    hull = []
    for di in range(int(min_id), int(max_id)+1):
        if di != max_domain_id:
            points = edge_pos[edge_pos[:, 0] == di][:, 1:3]
            pore_area = (len(points)/n_edges)*np.power(voxcel_lx, 2)
            chull = ConvexHull(points)
            hull.append([pore_area, chull.volume])

    hull = np.array(hull)
    hull_ratio = hull[:, 0]/hull[:, 1]

    fig, ax = plt.subplots()

    hist, bin_edges = np.histogram(hull_ratio, density=True)
    start = bin_edges[0]
    x = (bin_edges[:-1] + (bin_edges[1]-bin_edges[2])/2)/lx
    plt.plot(x, hist, lw=1, label='$l_p/l_r = {}$'.format(np.round(sigma/vs, 3)))
    plt.xlabel("$A_{p}/ A_{ch}$", size=15)
    plt.ylabel("P", size=15)
    plt.xlim([0, 1])
    plt.ylim([0, 10])
    plt.savefig("compare_convex_hull.pdf")

    # 3.2 PCA
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    xlambda = []
    for di in range(int(min_id), int(max_id)+1):
        if di != max_domain_id:
            points = edge_pos[edge_pos[:, 0] == di][:, 1:3]
            pca.fit(points)
            xlambda.append(pca.explained_variance_ratio_[0])

    fig, ax = plt.subplots()

    hist, bin_edges = np.histogram(xlambda, density=True)
    start = bin_edges[0]
    x = (bin_edges[:-1] + (bin_edges[1]-bin_edges[2])/2)/lx
    plt.plot(x, hist, lw=1, label='$l_p/l_r = {}$'.format(np.round(sigma/vs, 3)))
    plt.xlabel("explained variance ratio largest component", size=15)
    plt.ylabel("P", size=15)
    plt.xlim([0, 1])
    plt.ylim([0, 10])
    plt.savefig("pca_lambda.pdf")

    # 3.3 mean curvature

    # find out surface voxcels edgepoints and their sequence
    # smooth curve
    # calculate tangent vectors by numerical differentiation
    # calculate derivative of the tangent vectors to obtain curvature

    # list of surface voxels, and their positions
    # for domain in domains:
    #    domain_dg.append(list(map(get_ete_degree, domain)))
