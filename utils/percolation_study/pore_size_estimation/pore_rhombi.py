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
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA

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


def ana_run(fdir, fpos, forient, fbox, state):

    with open('param_pore_rhombi.yaml') as f:
        param = yaml.load(f, Loader=SafeLoader)

    print("read data...")

    box_all = np.fromfile("{}/{}".format(fdir, fbox))
    box_x_center = box_all[0]
    box_y_center = box_all[1]
    blx = box_all[3]
    box_ly = box_all[4]

    pos = np.fromfile("{}/{}".format(fdir, fpos))
    pos = np.reshape(pos, (-1, 3))
    pos = pos[:, :2]

    orient = np.fromfile("{}/{}".format(fdir, forient))
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
    pore_areas, domain_lengths, domains, G  = pt.get_pore_volume(voxcels)

##############################

    
    #phi_particles =0.125    
    #phi_all_pores = 1 - phi_particles 
   
    arr=[]
    for di, domain in enumerate(domains):
        for coord_vi in domain:
            arr.append([di, coord_vi[0],coord_vi[1]])

    arr=np.array(arr)

    min_id = np.min(arr[:,0])
    max_id = np.max(arr[:,0])
    uniques, counts = np.unique(arr[:,0], return_counts=True)
    max_domain_id = uniques[np.argmax(counts)]

    def unravel_pbc(arr,G, box, max_domain_id,min_id,max_id):
        new_voxcel_pos = []
        #random.seed(10)
        for di in range(int(min_id),int(max_id)+1): 
            if di!=max_domain_id:
                coords = arr[arr[:,0]==di][:,1:3]
                old_coords = [ (i,j) for i,j in coords]
                rand_c = coords[random.randint(0,len(coords)-1)]
                ci = (rand_c[0],rand_c[1])
                new_coords = [ci]   
                old_coords.remove(ci)
                
                old_coords, new_coords = stitch_cluster(G,ci, old_coords,new_coords)
                
                for entry in new_coords:
                    nni = entry 
                    vx = box.origin[0] + voxcel_lx*nni[0] + voxcel_lx/2 
                    vy = box.origin[1] + voxcel_lx*nni[1] + voxcel_lx/2 
                    new_voxcel_pos.append([di,vx,vy])

        new_voxcel_pos = np.array(new_voxcel_pos)

        new_edge_pos = []
        for di in range(int(min_id),int(max_id)+1): 
            if di!=max_domain_id:
                coords = voxcel_pos[new_voxcel_pos[:,0]==di][:,1:3]
                for coord_vi in coords:
                    vert_i = get_vertices(coord_vi,axes)
                    for vi in vert_i:
                        new_edge_pos.append([di, vi[0],vi[1]])
               
        new_edge_pos=np.array(new_edge_pos)

        return new_voxcel_pos, new_edge_pos 

    # stitch together cluster over pbcs by adapting voxcel pos and edge pos 
    new_voxcel_pos, new_edge_pos  = unravel_pbc(arr,G, box, max_domain_id,min_id,max_id)
   
    # draw 
    #draw_pos(voxcel_pos, blx, "voxcels_shifted.png",max_id,min_id,max_domain_id,axes)

    # 1. pore area histogram 
    def get_pore_area(pore_areas,particle_area):
        pore_area_excluded_void = np.sort(pore_areas)[:-1]
        hist, bin_edges = np.histogram(pore_area_excluded_void,density=True)
        start=bin_edges[0]
        x=(bin_edges[:-1] + (bin_edges[1]-bin_edges[2])/2)/particles_area
        return x,hist 
    # results: x, hist 

    x_area,hist_area = get_pore_area(pore_areas,particles.area)

    # 2. pore circumference histogram 
    def get_pore_circumference(G, domains, voxcels_lx):

        def ete_degree(G):  
            degree = defaultdict(int)
            for u,v in ((u,v) for u,v,d in G.edges(data=True) if d['ete']==1): 
                degree[u]+=1
                degree[v]+=1
            return degree

        degree_dict=ete_degree(G)

        def get_ete_degree(u):
            return degree_dict[u]

        domain_dg = []
        for domain in domains:
            domain_dg.append(list(map(get_ete_degree, domain)))

        def get_circumference(domain_dgi):
            arr=(n_edges - np.array(domain_dgi))*voxcel_lx
            cf = np.sum(arr)
            return cf

        circumferences = [ get_circumference(dd) for dd in domain_dg]

        pore_cf_excluded_void = np.sort(circumferences)[:-1]
        hist, bin_edges = np.histogram(pore_cf_excluded_void,density=True)
        start=bin_edges[0]
        x=(bin_edges[:-1] + (bin_edges[1]-bin_edges[2])/2)/lx

        return x,hist 

    x_cf, hist_cf = get_pore_circumference(G, domains, voxcels.lx)
    
    # 3. pore assymmetry 

    # 3.1 compare convex hull with actuall hull: measure of concavity
    # if A/Ahull = 1 -- convex cluster, if A/Ahull <1 concave 
    def estimate_convexity_wchull(new_edge_points, min_id, max_id, max_domain_id, n_edges, voxcel_lx):
        hull = []
        for di in range(int(min_id),int(max_id)+1): 
            if di!= max_domain_id:
                points = new_edge_pos[new_edge_pos[:,0]==di][:,1:3]
                pore_area = (len(points)/n_edges)*np.power(voxcel_lx,2)
                chull=ConvexHull(points)
                hull.append([pore_area, chull.volume])

        hull = np.array(hull)
        hull_ratio = hull[:,0]/hull[:,1]

        hist, bin_edges = np.histogram(hull_ratio,density=True)
        start = bin_edges[0]
        x = bin_edges[:-1] + (bin_edges[1]-bin_edges[2])/2

        return x, hist 
    
    x_ratio, hist_ratio = estimate_convexity_wchull(new_edge_points, min_id, max_id, max_domain_id, particles.N_edges, voxcels.lx):

    # 3.2 PCA 
    def estimate_asymmery_wpca(new_edge_points, min_id,max_id,max_domain_id):
        pca = PCA(n_components=2)
        xlambda = []
        for di in range(int(min_id),int(max_id)+1): 
            if di!=max_domain_id:
                points = new_edge_pos[new_edge_pos[:,0]==di][:,1:3]
                pca.fit(points)
                xlambda.append(pca.explained_variance_ratio_[0])


        hist, bin_edges = np.histogram(xlambda,density=True)
        start=bin_edges[0]
        x = bin_edges[:-1] + (bin_edges[1]-bin_edges[2])/2

        return x, hist 

    x_explained_variance, hist_explained_variance = estimate_asymmery_wpca(new_edge_points, min_id,max_id,max_domain_id)



##############################

    output_frame = param['output_frame']
    if output_frame:
        frame_name = param['frame_name']
        draw(particles, voxcels, box, cells, frame_name)







    '''
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
    ''' 

