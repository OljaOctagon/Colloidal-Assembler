import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib as mpl
import glob
from PyQt5 import QtGui, uic
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas


checkpoints= glob.glob("Box*.bin")
check_point_values = np.sort([ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])

check_point_values = check_point_values

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

Lx=1.0
Ly=2.0
a=0.2
b=0.8
radius=0.1

particle_patches = get_patches(Lx,Ly,a,b)

fig,ax = plt.subplots()
ax.set_aspect('equal', 'box')
plt.xlim((-1,55))
plt.ylim((-1,55))
plt.grid()

def start():
    def run(i):
        plt.cla()
        plt.xlim((-1,55))
        plt.ylim((-1,55))
        plt.grid()

        val = check_point_values[i]
        pos_i = np.fromfile("positions_{}.bin".format(val))
        pos_i = np.reshape(pos_i, (-1,3))
        pos_i = pos_i[:,:2]
        orient_i = np.fromfile("orientations_{}.bin".format(val))
        orient_i = np.reshape(orient_i, (-1,5))[:,4]
        #ax.scatter(pos_i[:,0], pos_i[:,1],s=1)

        N=len(pos_i)
        for j in range(N):
            rect = patches.Rectangle((-Lx/2,-Ly/2),
            Lx,Ly, linewidth=0.5, edgecolor='k',facecolor='#9AE0D1', alpha=0.7)
            theta = (orient_i[j]*180)/(np.pi)
            r = mpl.transforms.Affine2D().rotate_deg_around(0,0,theta)
            t = mpl.transforms.Affine2D().translate(pos_i[j,0],pos_i[j,1])

            #patch_1 = patches.Circle((particle_patches[0,0],particle_patches[0,1]),radius=radius,
            #facecolor='r')
            #patch_2 = patches.Circle((particle_patches[1,0],particle_patches[1,1]),radius=radius,
            #facecolor='r')
            tra = r + t + ax.transData
            rect.set_transform(tra)
            #patch_1.set_transform(tra)
            #patch_2.set_transform(tra)

            ax.add_patch(rect)
            #ax.add_patch(patch_1)
            #ax.add_patch(patch_2)

    def stop():
         ani.event_source.stop()

    def borr():
        plt.clf()
        canvas.draw()

    def anim():
        ani.event_source.start()

    window.resume.clicked.connect(anim)
    window.pause.clicked.connect(stop)
    window.clean.clicked.connect(borr)

    ani = animation.FuncAnimation(fig, run,repeat=False)
    canvas.draw()

layout=QtGui.QVBoxLayout()
fig=plt.figure()
canvas=FigureCanvas(fig)
layout.addWidget(canvas)

app = QtGui.QApplication(sys.argv)
window = uic.loadUi("animacion.ui")

window.start.clicked.connect(start)
window.widget.setLayout(layout)
window.show()
sys.exit(app.exec_())

#ani = FuncAnimation(fig, animate, frames=100)
#plt.show()
