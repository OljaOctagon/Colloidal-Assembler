import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib as mpl
# Set lenghts

#checkpoints= glob.glob("Box*.bin")
#check_point_values = np.sort([ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])

check_point_values = [200000]

for val in check_point_values:
	pos_i = np.fromfile("positions_{}.bin".format(val))
	pos_i = np.reshape(pos_i, (-1,3))
	pos_i = pos_i[:,:2]
	orient_i = np.fromfile("orientations_{}.bin".format(val))
	orient_i = np.reshape(orient_i, (-1,5))[:,4]

	N=len(pos_i)
	fig,ax = plt.subplots()
	ax.set_aspect('equal', 'box')

	for i in range(N):
		rect = patches.Rectangle((-0.5,-1.0),
		1.0,2.0, linewidth=0.5, edgecolor='k',facecolor='#9AE0D1', alpha=0.7)
		theta = (orient_i[i]*180)/(np.pi)
		r = mpl.transforms.Affine2D().rotate_deg_around(0,0,theta)
		t = mpl.transforms.Affine2D().translate(pos_i[i,0],pos_i[i,1])
		tra = r + t + ax.transData
		rect.set_transform(tra)
		ax.add_patch(rect)
	ax.scatter(pos_i[:,0], pos_i[:,1],s=1)
	plt.xlim((-1,55))
	plt.ylim((-1,55))
	plt.grid()
	plt.savefig("rect.png", dpi=600)
