import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib as mpl
import os
import glob
checkpoints= glob.glob("Box*.bin")
check_point_values = np.sort([ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])
Lx=1.0
Ly=2.0
a=0.2
b=0.8
radius=0.1

def get_patches(Lx,Ly,a,b):

	l = np.array([Lx/2,Ly/2])

	edges = np.zeros((4,2))
	edges[0] = np.array([-1,-1])*l
	edges[1] = np.array([ 1,-1])*l
	edges[2] = np.array([ 1, 1])*l
	edges[3] = np.array([-1, 1])*l

	patches = np.zeros((2,2))
	patches[0] = edges[0] + a*(edges[3]-edges[0])
	patches[1] = edges[2] + b*(edges[3]-edges[2])

	return patches

particle_patches = get_patches(Lx,Ly,a,b)
os.mkdir("./frames")
for j,val in enumerate(check_point_values):
	pos_i = np.fromfile("positions_{}.bin".format(val))
	pos_i = np.reshape(pos_i, (-1,3))
	pos_i = pos_i[:,:2]
	orient_i = np.fromfile("orientations_{}.bin".format(val))
	orient_i = np.reshape(orient_i, (-1,5))[:,4]

	N=len(pos_i)
	fig,ax = plt.subplots()
	ax.set_aspect('equal', 'box')

	for i in range(N):
		rect = patches.Rectangle((-Lx/2,-Ly/2),
		Lx,Ly, linewidth=0.5, edgecolor='k',facecolor='#9AE0D1', alpha=0.7)
		theta = (orient_i[i]*180)/(np.pi)
		r = mpl.transforms.Affine2D().rotate_deg_around(0,0,theta)
		t = mpl.transforms.Affine2D().translate(pos_i[i,0],pos_i[i,1])

		patch_1 = patches.Circle((particle_patches[0,0],particle_patches[0,1]),radius=radius,
		facecolor='r')
		patch_2 = patches.Circle((particle_patches[1,0],particle_patches[1,1]),radius=radius,
		facecolor='r')

		tra = r + t + ax.transData
		rect.set_transform(tra)
		patch_1.set_transform(tra)
		patch_2.set_transform(tra)

		ax.add_patch(rect)
		ax.add_patch(patch_1)
		ax.add_patch(patch_2)

	ax.scatter(pos_i[:,0], pos_i[:,1],s=1)
	plt.xlim((-1,55))
	plt.ylim((-1,50))
	plt.grid()
	plt.savefig("./frames/frame_{}.png".format(j), dpi=600)
	plt.close()
