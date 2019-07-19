import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib as mpl
import seaborn as sns
from matplotlib.widgets import Slider
import matplotlib.patches as mpatches
from matplotlib import rc 

sns.set(style='white')
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

colormap = plt.get_cmap('hsv')
#norm = mpl.colors.Normalize(0,np.pi/3.)

def get_pdist(df):
    Lx = df.Lx.unique()[0]
    Ly = df.Ly.unique()[0]
    l = len(df.x)
    pdist = np.zeros((l,l,2))

    p = np.reshape(df.x, (l,1))
    p = p - p.transpose()
    pdist[:,:,0] = p - Lx*np.rint(p/Lx)

    p = np.reshape(df.y, (l,1))
    p = p - p.transpose()
    pdist[:,:,1] = p - Ly*np.rint(p/Ly)

    N_pdist = np.sqrt( np.power( pdist[:,:,0], 2) + np.power( pdist[:,:,1], 2))
    return pdist, N_pdist

def hexatic_i(i,theta):
    return np.mean(np.exp(6*1j*theta))

def last_frame(df):
    last_time = df.time.values[-1]
    df = df[df.time == last_time ]
    return df
    

#------------------------------------------------

'''
import os 
import re
dir_list = next(os.walk('.'))[1]
parameter_names = ['mu', 'energy', 'patch', 'pressure', 'nrun']
filen = "/cluster_center_of_mass.dat"

# read in data and fill DataFrame
init=False
for dir in dir_list:
    parameter = re.findall(r"[-+]?\d*\.\d+|\d+", dir)
    parameter = list(map(float, parameter))
    file_trace = dir+filen
    print(file_trace)
 
    df = pd.read_csv(file_trace, header=None, delim_whitespace=True,
        names=['time','N','Lx', 'Ly', 'id', 'x', 'y'])

    df = last_frame(df)

    df = df.assign(patch = lambda x: parameter[2])
    df = df.assign(pressure = lambda x: parameter[3])
    df = df.assign(nrun = lambda x: parameter[4])

    if init==True: 
        data = data.append(df)

    if init==False:
        data = df 
        init=True



data.to_pickle("./dummy.pickle")
'''

pressure=100
patch=0.7
udata = pd.read_pickle("dummy.pickle")
data_sub = udata[ udata.pressure == pressure]
data_sub  = data_sub[data_sub.patch == patch]

Angle = []
Nruns = data_sub.nrun.unique()

for nrun in Nruns:
    data_sub_n = data_sub[data_sub.nrun == nrun]
    pdist, N_pdist = get_pdist(data_sub_n)
    sorted_dist = np.sort(N_pdist)[:,1:12]
    N=int(data_sub_n.N.unique()[0])
    for id in range(N):
        counts, bin_edges = np.histogram(sorted_dist[id,:])
        nneigh = np.argmax(counts)+np.argmin(counts[np.argmax(counts):])
        r_cutoff = bin_edges[nneigh+1]
        nn_id = np.argwhere(((N_pdist[id,:]<r_cutoff) 
            & (N_pdist[id,:]>0))).flatten()    
        
        theta = np.arccos(pdist[id, nn_id,0]/N_pdist[id,nn_id])
        psi = hexatic_i(id, theta)
        angle = np.arctan2(psi.real,psi.imag) + np.pi
        if angle > np.pi/3.:
            angle = angle - (np.pi/3.)

        Angle.append(angle)

data_sub['angle'] = Angle 

from pylab import plot, show, figure, scatter, axes, draw, subplots_adjust

fig = figure()
ax = fig.add_subplot(111)
subplots_adjust(left=0.15, bottom=0.25)

val0 = 1
x0 = data_sub[data_sub.nrun == val0].x.values
y0 = data_sub[data_sub.nrun == val0].y.values
c0 = data_sub[data_sub.nrun == val0].angle.values

scat = scatter(x0, y0, c=c0, 
    cmap=colormap, 
    label="P = "+str(pressure))

#plt.legend(prop={'size':15}, fancybox = True, loc='upper right')
#plt.xlabel("x", size=22)
#plt.ylabel("y", size=22)

axcolor = 'lightgoldenrodyellow'
ax_slider = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
run_slider = Slider(ax_slider, 'run', 1,8, valinit=val0, valstep=1)

def update(val):
    runi = run_slider.val
    ds = data_sub[data_sub.nrun == runi]
    xyi = ds[['x','y']].values
    ci = ds['angle'].values

    print(xyi)
    scat.set_array(ci)
    scat.set_offsets(xyi)
    draw()

run_slider.on_changed(update)
show(scat)


'''
fig,ax = plt.subplots()
plt.scatter(df_sub.x.values, df_sub.y.values, c=Angle, cmap=colormap, label="\displaystyle P = "+str(p))
plt.legend(prop={'size':15}, fancybox = True, bbox_to_anchor=(0.5, 1.3))
plt.show()
'''