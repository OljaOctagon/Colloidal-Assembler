import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib as mpl
import glob
from matplotlib.animation import FuncAnimation


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

def on_press(event):
    if event.key.isspace():
        if anim.running:
            anim.event_source.stop()
        else:
            anim.event_source.start()
        anim.running ^= True
    elif event.key == 'left':
        anim.direction = -1
    elif event.key == 'right':
        anim.direction = +1

    # Manually update the plot
    if event.key in ['left','right']:
        t = anim.frame_seq.next()
        update_plot(t)
        plt.draw()

def update_time():
    t = 0
    t_max = 10
    while t<t_max:
        t += anim.direction
        yield t

def update_plot(i):
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

    N=len(pos_i)
    for j in range(N):
        rect = patches.Rectangle((-Lx/2,-Ly/2),
        Lx,Ly, linewidth=0.5, edgecolor='k',facecolor='#9AE0D1', alpha=0.7)
        theta = (orient_i[j]*180)/(np.pi)
        r = mpl.transforms.Affine2D().rotate_deg_around(0,0,theta)
        t = mpl.transforms.Affine2D().translate(pos_i[j,0],pos_i[j,1])
        tra = r + t + ax.transData
        rect.set_transform(tra)
        ax.add_patch(rect)

#fig.canvas.mpl_connect('key_press_event', on_press)
anim = FuncAnimation(fig, update_plot, repeat=True)
#anim.running = True
#anim.direction = +1
plt.show()
