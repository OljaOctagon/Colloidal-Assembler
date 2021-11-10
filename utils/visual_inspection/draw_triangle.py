import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib as mpl
import os
import glob
import argparse
import matplotlib.style as style
mpl.rcParams['font.family'] = "sans-serif"

class Triangle:
    # needs self.vertices and self.patch_pos 
    def __init__(self, id, coord, theta, type, delta, radius):
        self.id = id
        self.center = coord
        self.vertices = np.zeros((3, 2))
        
        self.calc_vertices(theta)
        self.calc_patch_pos(type)

    def calc_vertices(self, theta):
        # how do i even python
        L = 1
        alpha = (60 * np.pi) / 180.0
        sinus = L * np.sin(alpha) * 1 / 3

        self.vertices[0,0] = - L/2
        self.vertices[0,1] = - sinus

        self.vertices[1,0] = L/2
        self.vertices[1,1] = - sinus

        self.vertices[2,0] = 0
        self.vertices[2,1] = 2 * sinus

        c, s = np.cos(theta), np.sin(theta) 
        R = np.array(((c, -s), (s, c)))

        self.vertices[0] = R.dot(self.vertices[0])
        self.vertices[1] = R.dot(self.vertices[1])
        self.vertices[2] = R.dot(self.vertices[2])

        self.vertices[:,0] += self.center[0]
        self.vertices[:,1] += self.center[1]

    def calc_patch_pos(self, type):
        if(type=="6patch"):
                #six_patch
                self.patch_pos = np.zeros((6,2))
                self.patch_pos[0] = self.vertices[0] + delta * (self.vertices[1]-self.vertices[0])
                self.patch_pos[1] = self.vertices[0] + (1 - delta) * (self.vertices[1]-self.vertices[0])
                self.patch_pos[2] = self.vertices[1] + delta * (self.vertices[2]-self.vertices[1])
                self.patch_pos[3] = self.vertices[1] + (1 - delta) * (self.vertices[2]-self.vertices[1])
                self.patch_pos[4] = self.vertices[2] + delta * (self.vertices[0]-self.vertices[2])
                self.patch_pos[5] = self.vertices[2] + (1 - delta) * (self.vertices[0]-self.vertices[2])

        elif(type=="3asym"):
                #three_asymm
                self.patch_pos = np.zeros((3,2))
                self.patch_pos[0] = self.vertices[0] + delta * (self.vertices[1]-self.vertices[0])
                self.patch_pos[1] = self.vertices[1] + delta * (self.vertices[2]-self.vertices[1])
                self.patch_pos[2] = self.vertices[2] + delta * (self.vertices[0]-self.vertices[2])

        elif(type=="2nfc"):
                #two_neighbour_fixedcorner
                self.patch_pos = np.zeros((2,2))
                self.patch_pos[0] = self.vertices[0]
                self.patch_pos[1] = self.vertices[0] + delta * (self.vertices[1]-self.vertices[0])

        elif(type=="2ofc"):
                #two_opposite_fixedcorner
                self.patch_pos = np.zeros((2,2))
                self.patch_pos[0] = self.vertices[0]
                self.patch_pos[1] = self.vertices[1] + delta * (self.vertices[2]-self.vertices[1])

if __name__ == '__main__':

    # command line
    parser = argparse.ArgumentParser(description='Particle drawing methods')
    parser.add_argument('-ptype', type=str, choices=['6patch', '3asym', '2nfc', '2ofc'])
    parser.add_argument('-delta', type=float)
    parser.add_argument('-radius', type=float)

    # import arguments
    args = parser.parse_args()
    delta = args.delta
    radius = args.radius
    ptype = args.ptype

    # get all check point values and sort them
    checkpoints= glob.glob("Box*.bin")
    check_point_values = np.sort(
    [ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])

    # make frame directory if it doesn't exist
    if not os.path.isdir("./frames"):
        os.mkdir("./frames")

    for j,val in enumerate(check_point_values[-1:]):
        position_list = np.fromfile("positions_{}.bin".format(val))
        position_list = np.reshape(position_list, (-1,3))
        position_list = position_list[:,:2]
        orientation_list = np.fromfile("orientations_{}.bin".format(val))
        orientation_list = np.reshape(orientation_list, (-1,5))[:,4]

        N=len(position_list)
        triangles = np.ndarray((N),dtype=Triangle)

        fig,ax = plt.subplots()
        ax.set_aspect('equal', 'box')
        for j in range(N):
            triangles[j] = Triangle(j, position_list[j], orientation_list[j], ptype, delta, radius)
            polygon = patches.Polygon(triangles[j].vertices, linewidth=0.5, edgecolor='k',facecolor='g', alpha=0.7)
            ax.add_patch(polygon)

            for patch_id in range(len(triangles[j].patch_pos)):
                patch = patches.Circle((triangles[j].patch_pos[patch_id,0],
                                        triangles[j].patch_pos[patch_id,1]),
                                        radius=radius,facecolor='r')
                ax.add_patch(patch)

        ax.set_title("Frame {}".format(j))
        plt.xlim((-1,30))
        plt.ylim((-1,40))
        plt.axis("equal")
        plt.axis('off')
        plt.savefig("./frames/frame_{}.png".format(j), dpi=500)
        plt.close()