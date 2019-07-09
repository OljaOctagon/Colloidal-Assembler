
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import rc 
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


df = pd.read_csv("b.txt", header=None, delim_whitespace=True, 
    names=['time', 'nmoves', 'nbroken', 'phi_i', 'patch_ix', 'patch_iy', 'phi_j','patch_jx', 'patch_jy','wx', 'wy'])
df['dpatch_x'] = df.patch_ix - df.patch_jx
df['dpatch_y'] = df.patch_iy - df.patch_jy
df['dphi'] = np.arctan2(df.wx,df.wy)

points = df[['dpatch_x', 'dpatch_y', 'dphi']].values

from scipy.spatial import ConvexHull  
from matplotlib.patches import Polygon

hull = ConvexHull(points)

fig = plt.figure()
ax = fig.gca(projection='3d',proj_type = 'ortho')
ax.grid(False)
#ax.set_aspect(0.5,'datalim')
ax.set_facecolor('#1e0768')
fig.patch.set_facecolor('#1e0768')
# Hide axes ticks
#ax.set_xticks([])
#ax.set_yticks([])
#ax.set_zticks([])

ax.set_xlim([0,0.1])
ax.set_ylim([-0.1,0.1])
ax.set_zlim([-np.pi,np.pi])

ax.set_xlabel('dx patch', size=14)
ax.set_ylabel('dy patch', size=14)
ax.set_zlabel('$\displaystyle \omega$', size=25)

ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

for s in hull.simplices:
    s = np.append(s, s[0])  
    ax.plot_trisurf(points[s, 0], points[s, 1], points[s, 2], lw=0.2, edgecolor='white', color='#efe010', alpha=0.9)
    #ax.plot(points[s, 0], points[s, 1], points[s, 2], "r-")
##00ccff'
ax.bar3d(0, -0.1, -np.pi, 0.1, 0.2, 2*np.pi, edgecolor='#ed5642', color=(1,1,1,0), alpha=0.1)


ax.tick_params(axis='both', labelsize=16)
ax.set_xticks([0,0.1])
ax.set_yticks([-0.1,0.1])
ax.set_zticks([-np.pi,0,np.pi])
plt.tight_layout()

ax.set_zticklabels(['$\displaystyle -\pi$','0', '$\displaystyle \pi$'])

ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.xaxis.label.set_color('white')
ax.yaxis.label.set_color('white')
ax.zaxis.label.set_color('white')
ax.tick_params(axis='both', colors='white')


#xyw
ax.view_init(elev=30., azim=-60)
plt.savefig("one_patch_parallel_3axis_view.pdf", facecolor='#1e0768')

#xw
ax.view_init(elev=0., azim=-90)
ax.set_yticks([])
ax.set_ylabel('')
plt.savefig("one_patch_parallel_xw_view.pdf")

ax.set_yticks([-0.1,0.1])
ax.set_ylabel('dy patch', size=14)

#yw
ax.set_xticks([])
ax.set_xlabel('')
ax.view_init(elev=0, azim=0)
plt.savefig("one_patch_parallel_yw_view.pdf")

ax.set_xticks([0,0.1])
ax.set_xlabel('dx patch', size=14)

#xy
ax.set_aspect('0.5')
ax.set_zticks([])
ax.set_zlabel('')
ax.view_init(elev=90., azim=-90)
plt.savefig("one_patch_parallel_xy_view.pdf")

plt.show()


#ax.text(0.05,-0.12,-0.22, 'dx patch','x', size=16, color='r')
#ax.text(-0.02,0,-0.22, 'dy patch', 'y', size=16, color='r')
#ax.text(0.12,0,0,'$\omega$', 'z', size=16, color='r')

#ax.text(0,-0.12,-0.22,'0', 'x', size=16, color='r')
#ax.text(0.1,-0.12,-0.22,'0.1', 'x', size=16, color='r')

#ax.text(-0.02,-0.1,-0.22, '-0.1','y', size=16, color='r')
#ax.text(-0.02,0.1,-0.22,'0.1', 'y', size=16, color='r')

#ax.text(0.12,0,-0.2,'-0.2', 'z', size=16, color='r')
#ax.text(0.12,0,0.2,'0.2', 'z', size=16, color='r')


##efe010



'''
 print(pts)
hull = ConvexHull(pts)

#ax.plot(pts.T[0], pts.T[1], pts.T[2], "ko")

for s in hull.simplices:
    s = np.append(s, s[0])  # Here we cycle back to the first coordinate
    #ax.plot(pts[s, 0], pts[s, 1], pts[s, 2], "r-")
    ax.contourf(pts[s, 0], pts[s, 1], pts[s, 2], color='grey', alpha=0.3)

plt.show()

'''