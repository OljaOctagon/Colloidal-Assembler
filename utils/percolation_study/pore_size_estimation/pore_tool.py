import numpy as np 
import pandas as pd 
import networkx as nx 
import matplotlib as pyplot 
from collections import defaultdict 
import particle_tools as rt 

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import numba 

@numba.jit(nopython=True)
def numba_ray_tracing(x,y,poly):
    n = len(poly)
    inside = False
    p2x = 0.0
    p2y = 0.0
    xints = 0.0
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

class Cells:
    def __init__(self,lmin,box_lx,ndim):
        self.nx = int(np.floor(box_lx/lmin))
        self.lx = box_lx/self.nx 
        
        # generate neigbhour basis 
        en = np.array([-1,0,1])
        enl = len(en)

        self.ndim=ndim

        if self.ndim==2:
            cx,cy = np.meshgrid(en,en)
            ccx = np.reshape(cx,(1,enl*enl))
            ccy = np.reshape(cy,(1,enl*enl))
            arr = np.transpose(np.concatenate((ccx,ccy), axis=0))
            self.neighbour_cell_basis = [ (i,j) for i,j in arr]
                      
            xx,yy= np.meshgrid(range(self.nx),range(self.nx))
            x = np.reshape(xx,(1,self.nx*self.nx))
            y = np.reshape(yy,(1,self.nx*self.nx))
            arr = np.transpose(np.concatenate((x,y), axis=0))
            self.coords = [ (i,j) for i,j in arr]
            
        if self.ndim==3:

            cx,cy,cz = np.meshgrid(en,en,en)
            ccx = np.reshape(cx,(1,enl*enl*enl))
            ccy = np.reshape(cy,(1,enl*enl*enl))
            ccz = np.reshape(cz,(1,enl*enl*enl))
            arr = np.transpose(np.concatenate((ccx,ccy,ccz), axis=0))
            self.neighbour_cell_basis = [ (i,j,k) for i,j,k in arr]
                      
            xx,yy,zz = np.meshgrid(range(self.nx),range(self.nx),range(self.nx))
            x = np.reshape(xx,(1,self.nx*self.nx*self.nx))
            y = np.reshape(yy,(1,self.nx*self.nx*self.nx))
            z = np.reshape(zz,(1,self.nx*self.nx*self.nx))
            arr = np.transpose(np.concatenate((x,y,z), axis=0))
            self.coords = [ (i,j,k) for i,j,k in arr]     
        
    def list_generate(self, ppos, vcoords, vpos, origin):
        
        self.particle_list = defaultdict(list)
        for i,pos_i in enumerate(ppos):
            idx = np.zeros(self.ndim)
            for dim_i in range(self.ndim):
                idx[dim_i] = np.floor((pos_i[dim_i] - origin[dim_i])/self.lx).astype(int)
           
            self.particle_list[tuple(idx)].append(i)
        
        self.voxcel_list = defaultdict(list)
        for coord_i in vcoords:
            pos_i = vpos[coord_i]
            idx=np.zeros(self.ndim)
            for dim_i in range(self.ndim):
                idx[dim_i] = np.floor((pos_i[dim_i] - origin[dim_i])/self.lx).astype(int)

            self.voxcel_list[tuple(idx)].append(coord_i)
        
          
    def get_neighbour_coords(self, coord_i):
         # get particles from neighbour cells   
        nncell = []
        for ncell in self.neighbour_cell_basis: 
            celli = np.array(coord_i) + np.array(ncell)
            for dim_i in range(self.ndim):
                celli[dim_i] = celli[dim_i]%self.nx
            nncell.append(tuple(celli))      
        return nncell 

class Box:
    def __init__(self,origin,lx,ndim):
        self.ndim = ndim 
        self.origin = origin
        self.lx = lx 
        self.ly = lx 
        self.center = origin + self.lx/2*np.ones(self.ndim)
    
class Spheres:
    sigma=1
    radius=sigma/2 
    inner_radius = radius 
    outer_radius = radius 
    def __init__(self, pos,ndim): 
        self.pos = pos
        self.ndim = ndim 
       
    def get_pos(self,id_i):
        return self.pos[id_i]

    
class Voxcels:
    def __init__(self,lx,vxn,origin,ndim):
        self.lx=lx
        self.ndim=ndim
        self.outer_radius = lx*np.sqrt(self.ndim)/2
        self.inner_radius = lx/2   
        self.nx = vxn
        self.origin = origin
        self.volume = np.power(self.lx, self.ndim)
        self.N_edges = np.power(2,ndim)
        
        en = np.array([-1,0,1])
        enl = len(en)

        #xyz = np.meshgrid(en,en,en)
        #cxyz = (np.empty(),np.empty(),np.empty())

        #for dim_i in range(self.ndim):
        #    cxyz[dim_i] = np.reshape(cxyz[dim_i],(1,np.power(enl,3)))

        #arr=np.transpose(np.concatenate(cxyz, axis=0))

        if self.ndim == 3: 
            cx,cy,cz = np.meshgrid(en,en,en)
            ccx = np.reshape(cx,(1,enl*enl*enl))
            ccy = np.reshape(cy,(1,enl*enl*enl))
            ccz = np.reshape(cz,(1,enl*enl*enl))
            arr = np.transpose(np.concatenate((ccx,ccy,ccz), axis=0))
            self.neighbour_cell_basis = [ (i,j,k) for i,j,k in arr]
             
             # initalize voxcel positions 
            xx,yy,zz = np.meshgrid(range(vxn),range(vxn), range(vxn))
            vx = np.reshape(xx,(1,self.nx*self.nx*self.nx))
            vy = np.reshape(yy,(1,self.nx*self.nx*self.nx))
            vz = np.reshape(zz,(1,self.nx*self.nx*self.nx))
            vxyz=np.transpose(np.concatenate((vx,vy,vz), axis=0))
            self.coords = [(i,j,k) for i,j,k in vxyz]

            en=np.array([-1,1])
            enl = len(en)
            x,y,z = np.meshgrid(en,en,en)
            dx = np.reshape(x,(1,np.power(enl,self.ndim)))
            dy = np.reshape(x,(1,np.power(enl,self.ndim)))
            dz = np.reshape(x,(1,np.power(enl,self.ndim)))
            self.vertex_basis = np.transpose(np.concatenate((dx,dy), axis=0))

            self.axes = np.array([[1,0,0],[0,1,0],[0,0,1]])*self.lx
        
        
        if self.ndim == 2:
            cx,cy = np.meshgrid(en,en)
            ccx = np.reshape(cx,(1,enl*enl))
            ccy = np.reshape(cy,(1,enl*enl))
            arr = np.transpose(np.concatenate((ccx,ccy), axis=0))
            self.neighbour_cell_basis = [ (i,j) for i,j in arr]
             
             # initalize voxcel positions 
            xx,yy = np.meshgrid(range(vxn),range(vxn))
            vx = np.reshape(xx,(1,self.nx*self.nx))
            vy = np.reshape(yy,(1,self.nx*self.nx))
            vxy=np.transpose(np.concatenate((vx,vy), axis=0))
            self.coords = [ (i,j) for i,j in vxy ]

            en=np.array([-1,1])
            enl=len(en)
            x,y = np.meshgrid(en,en)
            dx = np.reshape(x,(1,np.power(enl,self.ndim)))
            dy = np.reshape(x,(1,np.power(enl,self.ndim)))
            self.vertex_basis = np.transpose(np.concatenate((dx,dy), axis=0))
            
            self.axes = np.array([[1,0],[0,1]])*self.lx
        

        self.pos = defaultdict(list)
        for ci in self.coords:
            self.pos[ci] = self.origin + self.lx*np.array(ci) + self.lx/2

        self.fill_state = defaultdict(list)
        for ci in self.coords: 
            self.fill_state[ci] = 0

        self.links = []
    
    # TODO only in 2D NOW!!! implement generally after debugging.
    def get_vertices(self, coord_i):
        p = self.pos[coord_i]
        
        def get_edge_points(p,axes,sign_p):
            vertex_n = np.zeros(2)
            vertex_n = p + sign_p[0]*axes[:,0]/2. + sign_p[1]*axes[:,1]/2.
            return vertex_n

        vertices = np.zeros((4,2))
        
        vertices[0] = get_edge_points(p,self.axes,np.array([-1,-1]))
        vertices[1] = get_edge_points(p,self.axes,np.array([+1,-1]))
        vertices[2] = get_edge_points(p,self.axes,np.array([+1,+1]))
        vertices[3] = get_edge_points(p,self.axes,np.array([-1,+1]))
        return vertices
    
    def set_to_filled(self, coord_i):
        self.fill_state[coord_i] = 1 
        
    def set_to_empty(self, coord_i):
         self.voxcel_state[coord_i] = 0 
         
    def get_neighbour_coords(self, coord_i):
         # get particles from neighbour cells   
        nncell = []
        for ncell in self.neighbour_cell_basis: 
            celli = np.array(coord_i) + np.array(ncell)
            for dim_i in range(self.ndim):
                celli[dim_i] = celli[dim_i]%self.nx
            nncell.append(tuple(celli))
        return nncell 

    def get_links(self):
        for coord_vi in self.coords:
            if self.fill_state[coord_vi] == 0:
                for coord_ni in self.get_neighbour_coords(coord_vi):
                    if (self.fill_state[coord_ni] == 0) and (coord_vi != coord_ni):
                        self.links.append((coord_vi,coord_ni))


class Rhombi:
    def __init__(self,pos,orient,lx,ndim):

        self.lx=lx
        self.ly = lx  
        self.ndim = 2
        self.N_axes=2
        self.Np = len(pos)
        self.pos = pos 
        self.orient = orient    
        self.alpha=np.pi/3
        self.h = self.lx * np.sin(self.alpha)
        self.a_x = self.ly * np.cos(self.alpha)

        #self.outer_radius = np.sqrt((self.lx + self.a_x) * (self.lx + self.a_x) + self.h * self.h)/2.
        #self.inner_radius = np.sqrt((self.lx - self.a_x) * (self.lx - self.a_x) + self.h * self.h)/2. + 0.15
    

        self.long_diagonal = self.lx*np.sqrt(2+2*np.cos(self.alpha))
        self.short_diagonal = self.lx*np.sqrt(2-2*np.cos(self.alpha))

        self.outer_radius = self.long_diagonal/2.
        self.inner_radius = self.lx*np.sin(self.alpha)/2.

        self.sigma = self.long_diagonal
        self.volume = self.lx*self.h
        self.N_edges = 4
        self.N_patches = 4 
       
        #en=np.array([-1,1])
        #enl=len(en)
        #x,y = np.meshgrid(en,en)
        #dx = np.reshape(x,(1,np.power(enl,self.ndim)))
        #dy = np.reshape(x,(1,np.power(enl,self.ndim)))
        #self.vertex_basis = np.transpose(np.concatenate((dx,dy), axis=0))

        def get_edge_points(p,ax_n,sign_p):
            vertex_n = np.zeros(2)
            vertex_n = p + sign_p[0]*ax_n[:,0]/2. + sign_p[1]*ax_n[:,1]/2.
            return vertex_n

        sin60 = np.sin(np.pi/3.)
        cos60 = np.cos(np.pi/3.)

        self.ax0 = np.array([[1,cos60],[0,sin60]])
        self.vertices = np.zeros((self.Np,4,2))
        self.ax_n = np.zeros((self.Np,2,2))

        def rotation_matrix(theta):
            rot_mat = np.zeros((2,2))
            rot_mat[0,0] = np.cos(theta) 
            rot_mat[0,1] = -np.sin(theta)
            rot_mat[1,0] = np.sin(theta)
            rot_mat[1,1] = np.cos(theta)

            return rot_mat

        def get_orient(v, rot_mat):
            return rot_mat.dot(v)

        for i in range(self.Np):
            p = self.pos[i]
            o = self.orient[i]

            rotmat_i = rotation_matrix(o)
            self.ax_n[i] = get_orient(self.ax0, rotmat_i)
            self.vertices[i,0] = get_edge_points(p,self.ax_n[i],np.array([-1,-1]))
            self.vertices[i,1] = get_edge_points(p,self.ax_n[i],np.array([+1,-1]))
            self.vertices[i,2] = get_edge_points(p,self.ax_n[i],np.array([+1,+1]))
            self.vertices[i,3] = get_edge_points(p,self.ax_n[i],np.array([-1,+1]))

    def point_intersect(self,rel_point,pi):
        intersect = 0
        point_rhs = np.zeros((self.ndim,1))

        #e01 = self.ax_n[coord_pi, :, 0]/self.lx 
        #e02 = self.ax_n[coord_pi, :, 1]/self.lx 

        e01 = (self.vertices[pi, 1] - self.vertices[pi, 0])/self.lx 
        e02 = (self.vertices[pi, 3] - self.vertices[pi, 0])/self.lx

        point_rhs = np.array([
        np.dot(rel_point,e01),
        np.dot(rel_point,e02)])

        length_point_rhs = np.linalg.norm(point_rhs)

        if length_point_rhs < 1:
            intersect=1

        return intersect 

def get_distance(v_pos,p_pos, blx):
    dist = p_pos - v_pos 
    dist = dist - blx*np.rint(dist/blx)
    return dist   


def get_vdistance(pos_v,pos_c,blx):
    
    m = len(pos_c)
    n = len(pos_v)
    ndim = pos_v.shape[1]
   
    pos_c_tiled = np.tile(pos_c,(n,1))
    pos_v_repeat = np.repeat(pos_v,m,axis=0)
    
    #print("tiled", pos_c_tiled)
    #print("repeat", pos_v_repeat)

    dist = pos_c_tiled - pos_v_repeat 
    dist = dist - blx*np.rint(dist/blx)
    ndist = np.linalg.norm(dist,axis=1)
    dist = np.reshape(dist,(n,m,ndim))
    ndist = np.reshape(ndist,(n,m))
    return dist, ndist 


def point_intersect(pdist, voxcels, coord_vi):
    
    intersect = 0
    vertices = voxcels.get_vertices(coord_vi)

    e01 = (vertices[1] - vertices[0])/voxcels.lx 
    e02 = (vertices[3] - vertices[0])/voxcels.lx 
    e03 = (vertices[4] - vertices[0])/voxcels.lx

    point_rhs = np.array([
        np.dot(pdist,e01),
        np.dot(pdist,e02),
        np.dot(pdist,e03)])
    
    length_point_rhs = np.linalg.norm(point_rhs)
    if length_point_rhs < 1:
        intersect=1

    return intersect 

def estimate_volume(voxcels, coord_vi, particles, coord_pi,blx):
    Ntrial = 100
    p_pos = np.reshape(particles.pos[coord_pi],(particles.ndim,1))
    total_intersect=0
    points = voxcels.lx*(0.5 - np.random.sample((Ntrial,voxcels.ndim)))

    for dim_i in range(voxcels.ndim):
        points[:,dim_i] =  points[:,dim_i] + voxcels.pos[coord_vi][dim_i]  


    pos_c = particles.pos[coord_pi]
    pos_c = np.reshape(pos_c, (1,-1))
    dist, ndist = get_vdistance(points,pos_c,blx) 

    total_intersect = len(ndist[ndist<particles.inner_radius].flatten())

    #for i,j in np.argwhere(ndist<particles.inner_radius)&(ndist<particles.outer_radius)):
    overlap_cand = np.argwhere((ndist>particles.inner_radius)&(ndist<particles.outer_radius))
    #print("oc", overlap_cand)
    polygon=particles.vertices[coord_pi]
    inside1 = [numba_ray_tracing(points[i,0], points[i,1], polygon) for i in overlap_cand[:,0] ]
    #print("inside", inside1)
    total_intersect += np.count_nonzero(inside1)
    #for i,j in np.argwhere(ndist<particles.inner_radius)&(ndist<particles.outer_radius)):
    #    point = Point(tuple(points[i]))
    #    p1 = tuple(particles.vertices[coord_pi,0])
    #    p2 = tuple(particles.vertices[coord_pi,1])
    #    p3 = tuple(particles.vertices[coord_pi,2])
    #    p4 = tuple(particles.vertices[coord_pi,3])

    #    polygon = Polygon([p1,p2,p3,p4])
    #    total_intersect += polygon.contains(point)

        #total_intersect+= particles.point_intersect(dist[i,j],coord_pi)

    overlap_volume = total_intersect/Ntrial 
    #print("overlap_volume ", overlap_volume)
    return overlap_volume


def get_overlap_volume(voxcels, coord_vi, particles, coord_pi, blx):
    v_pos = voxcels.pos[coord_vi]
    p_pos = particles.get_pos(coord_pi)
    overlap_volume_i = 0 
    # get distance 
    dist_vp = get_distance(v_pos,p_pos,blx)
    ldist_vp = np.linalg.norm(dist_vp)
                         
    # no overlap if distance larger than outer radius of voxcels + radius of sphere
    if ldist_vp > (voxcels.outer_radius + particles.radius):
        return overlap_volume_i 
                         
    # else, check inner radius of voxcel 
    else:
        # overlap, if distance smaller than inner radius:
        if ldist_vp < (voxcels.inner_radius + particles.radius):
            #print("overlap because dist smaller than inner radius")
            overlap_volume_i = estimate_volume(voxcels, coord_vi, particles, coord_pi, blx)
            #print("overlap volume MUST BE BIGGER THAN ZERO", overlap_volume_i)
            return overlap_volume_i 
                         
        # otherwise, detailed overlap check 
        else: 
            nearest_point_distance = (dist_vp/ldist_vp)*(1 - particles.radius)                    
            intersect = point_intersect(nearest_point_distance,voxcels, coord_vi)

            if intersect==0:
                return overlap_volume_i

            else:
                overlap_volume_i = estimate_volume(voxcels, coord_vi, particles, coord_pi, blx)
                return overlap_volume_i

def get_pore_volume(voxcels):
        G = nx.Graph()
        G.add_edges_from(voxcels.links)
        domains = list(nx.connected_components(G))
        domain_lengths = np.array([ len(domain) for domain in domains ])
        pore_volumes = voxcels.volume*domain_lengths
        return pore_volumes, domain_lengths




