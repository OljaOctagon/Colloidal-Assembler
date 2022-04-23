import numpy as np
import cv2 as cv
import imutils
from math import ceil

scale = 100

def get_orient(v, rot_mat):
    return rot_mat.dot(v)

def rotation_matrix(theta):
    rot_mat = np.zeros((2,2))

    rot_mat[0,0] = np.cos(theta) 
    rot_mat[0,1] = -np.sin(theta)
    rot_mat[1,0] = np.sin(theta)
    rot_mat[1,1] = np.cos(theta)

    return rot_mat

def get_rhombi_vertices(p, o):

    def get_edge_points(p,ax_n,sign_p):
        vertex_n = np.zeros(2)
        vertex_n = p + sign_p[0]*ax_n[:,0]/2. + sign_p[1]*ax_n[:,1]/2.
        return vertex_n

    sin60 = np.sin(np.pi/3.)
    cos60 = np.cos(np.pi/3.)

    ax0 = np.array([[1,cos60],[0,sin60]])
    vertices = np.zeros((4,2))
    ax_n = np.zeros((2,2))

    rotmat_i = rotation_matrix(o)
    ax_n = get_orient(ax0, rotmat_i)

    vertices[0] = get_edge_points(p,ax_n,np.array([-1,-1]))
    vertices[1] = get_edge_points(p,ax_n,np.array([+1,-1]))
    vertices[2] = get_edge_points(p,ax_n,np.array([+1,+1]))
    vertices[3] = get_edge_points(p,ax_n,np.array([-1,+1]))

    return vertices

def get_patch_positions(vertices):
    res = np.zeros((4, 2))

    res[0] =  vertices[0] + ((vertices[1] - vertices[0]) * .8)
    res[1] =  vertices[2] + ((vertices[1] - vertices[2]) * .8)
    res[2] =  vertices[3] + ((vertices[2] - vertices[3]) * .8)
    res[3] =  vertices[3] + ((vertices[0] - vertices[3]) * .8)
    return res 


def read_binaries():
    box_all = np.fromfile("box.bin")
    box_x_center = box_all[0]  
    box_y_center = box_all[1]    

    box_lx = box_all[3]
    box_ly = box_all[4]

    pos_i = np.fromfile("positions.bin")
    pos_i = np.reshape(pos_i, (-1,3))
    pos_i = pos_i[:,:2]

    orient_i = np.fromfile("orientations.bin")
    orient_i = np.reshape(orient_i, (-1,5))[:,4]

    return box_x_center, box_y_center, box_lx, box_ly, pos_i, orient_i 

if __name__ == '__main__':
    box_x_center, box_y_center, box_lx, box_ly, pos_i, orient_i  = read_binaries()
    img = np.full((ceil(box_lx * scale), ceil(box_ly * scale), 3), 255, np.uint8)

    for p, o in zip(pos_i, orient_i):
        vert = get_rhombi_vertices(p, o)
        patches = get_patch_positions(vert)
        cv.fillPoly(img, np.int32([vert * scale]), (0,0,0))
        for patch in patches:
            cv.circle(img,(ceil(patch[0] * scale), ceil(patch[1] * scale)), ceil(.1 * scale), (0,0,0), -1)

    #cimg = cv.cvtColor(img, cv.COLOR_BGR2GRAY)

    outsize = (5000,5000)
    out = cv.resize(img, outsize)
    cv.imwrite('out.png', out)










