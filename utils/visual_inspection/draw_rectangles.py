import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib as mpl
import os
import glob
import networkx as nx

def read_bonds(filen):
	first_line_pair = [0,0,0,0]
	cut=False
	with open(filen, 'r') as f:
		network_list = []
		for line in f:
			if "#" in line:
				network_list.append([])
				first_line_pair = [0,0,0,0]
				cut=False

			else:
				line_counter=len(network_list[-1])
				pairs = list(map(int, line.split(" ")))
				if pairs == first_line_pair or cut==True:
					cut=True
				else:
					network_list[-1].append(np.array(pairs))

				if line_counter == 0:
					first_line_pair = pairs
	network_list = [ np.array(item) for item in network_list]

	return network_list

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

def get_hexcolor(i, cmap):
	rgb = cmap(i)[:3]
	return mpl.colors.rgb2hex(rgb)

def get_domain_colors(N_particles, bond_arr, length_color_dict):
	G = nx.Graph()
	G.add_edges_from(bond_arr[:,:2])
	domains = list(nx.connected_components(G))
	domain_colors = [length_color_dict[0]]*N
	for cluster in domains:
		length = len(cluster)
		if length >=5:
			length=5
		for particle in cluster:
			domain_colors[particle] = length_color_dict[length]

	return domain_colors

if __name__ == '__main__':
	# get all check point values and sort them
	checkpoints= glob.glob("Box*.bin")
	check_point_values = np.sort(
	[ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])

	# import the bonds file:
	# file format:
	# ------------------------
	# #new time
	# particle_id1 particle_id2  patch_id1 patch_id2]
	#  ....
	# --------------------------

	pn_file = "patch_network.dat"

	# colormap for cluster size
	cmap = plt.cm.get_cmap('cividis', 6)
	hex_color = [ get_hexcolor(i, cmap) for i in range(6)]
	hex_color[4] = '#8A2BE2'
	length_color_dict = dict(zip(np.arange(6),hex_color))

	# network_arr format: network_arr.shape = ( frame_i, bond_rows_frame_i )
	network_arr = read_bonds("patch_network.dat")
	# patch position calculation
	Lx=1.0
	Ly=2.0
	a=0.25
	b=0.5
	radius=0.2
	particle_patches = get_patches(Lx,Ly,a,b)

	patch_color_dict = {0:'red', 2:'yellow'}

	# make frame directory if it doesn't exist
	if not os.path.isdir("./frames"):
		os.mkdir("./frames")

	for j,val in enumerate(check_point_values[1:]):
		pos_i = np.fromfile("positions_{}.bin".format(val))
		pos_i = np.reshape(pos_i, (-1,3))
		pos_i = pos_i[:,:2]
		orient_i = np.fromfile("orientations_{}.bin".format(val))
		orient_i = np.reshape(orient_i, (-1,5))[:,4]

		fb = open("patch_energy_{}.bin".format(val), "rb")
		patch_i = np.fromfile(fb, dtype=np.int32)
		patch_i = np.reshape(patch_i, (-1,4))[:,:2]
		print(patch_i)
		N=len(pos_i)
		patch_i = patch_i[:N]

		fig,ax = plt.subplots()
		ax.set_aspect('equal', 'box')

		# domain_colors.shape = (N,) ( colors per particle )
		domain_colors = get_domain_colors(N, network_arr[j], length_color_dict)

		for i in range(N):
			print(i, patch_i[i,0], patch_i[i,1])
			rect = patches.Rectangle((-Lx/2,-Ly/2),
			Lx,Ly, linewidth=0.5, edgecolor='k',facecolor=domain_colors[i], alpha=0.7)
			theta = (orient_i[i]*180)/(np.pi)
			r = mpl.transforms.Affine2D().rotate_deg_around(0,0,theta)
			t = mpl.transforms.Affine2D().translate(pos_i[i,0],pos_i[i,1])

			patch_1 = patches.Circle((particle_patches[0,0],particle_patches[0,1]),radius=radius,
			facecolor=patch_color_dict[patch_i[i,0]])
			patch_2 = patches.Circle((particle_patches[1,0],particle_patches[1,1]),radius=radius,
			facecolor=patch_color_dict[patch_i[i,1]])

			tra = r + t + ax.transData
			rect.set_transform(tra)
			patch_1.set_transform(tra)
			patch_2.set_transform(tra)

			ax.add_patch(rect)
			ax.add_patch(patch_1)
			ax.add_patch(patch_2)

		ax.scatter(pos_i[:,0], pos_i[:,1],s=1)
		ax.set_title("Frame {}".format(j))
		plt.xlim((-1,80))
		plt.ylim((-1,70))
		plt.grid()
		plt.savefig("./frames/frame_{}.png".format(j), dpi=300)
		plt.close()
