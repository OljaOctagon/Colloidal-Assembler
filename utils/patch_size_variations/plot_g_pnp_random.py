import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
import matplotlib as mpl
import matplotlib.style as style

style.use('seaborn-ticks') 
mpl.rcParams['font.family'] = "sans-serif"
#sns.set_context('poster')
plt.rcParams['axes.axisbelow'] = True

plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['figure.titlesize'] = 15 

blue_c ='#9999FF'
red_c ='#FF9999'
purple_c='#C17DCB'
green_c='#61DCB7'


savepath = "../plots/constr/"
df = pd.read_pickle("../analysis_data/constr_mp.pkl")

###############

mus = np.sort(np.array(df.mu.unique()))
deltas = np.sort(np.array(df.delta.unique()))
radii = np.sort(np.array(df.r_patch.unique()))
len_g = len(df['g'].iloc[0])
print('all mus : ', mus)
print('all deltas : ', deltas)
print('all radii : ', radii)
mu_dict ={x:i for i,x in enumerate(mus)}
del_dict={x:i for i,x in enumerate(deltas)}
rad_dict={x:i for i,x in enumerate(radii)}

#################

maxbonds  = np.zeros((len(mus), len(radii), len(deltas)))
polybonds = np.zeros((len(mus), len(radii), len(deltas)))
polyratio = np.zeros((len(mus), len(radii), len(deltas)))
psi       = np.zeros((len(mus), len(radii), len(deltas), 160))
psi_avg   = np.zeros((len(mus), len(radii), len(deltas)))
psi_std   = np.zeros((len(mus), len(radii), len(deltas)))
bonds     = np.zeros((len(mus), len(radii), len(deltas), 2)) #rename
g         = np.zeros((len(mus), len(radii), len(deltas), 160, len_g))
g_avg     = np.zeros((len(mus), len(radii), len(deltas), len_g))
maxcluster_avg= np.zeros((len(mus), len(radii), len(deltas)))

##################

for properties, ensemble in df.groupby(['mu', 'delta', 'r_patch']):
    mu = properties[0]
    delta = properties[1]
    radius = properties[2]
    i_m = mu_dict[mu]
    i_d = del_dict[delta]
    print(del_dict[delta], delta)
    i_r = rad_dict[radius]
    
    maxbonds[i_m,i_r,i_d] = np.mean(ensemble['max_bonds'].values)
    
    psi[i_m,i_r,i_d]     = ensemble['psi'].values
    psi_avg[i_m,i_r,i_d] = np.mean(psi[i_m,i_r,i_d])
    psi_std[i_m,i_r,i_d] = np.std(psi[i_m,i_r,i_d])
    
    g[i_m,i_r,i_d] = np.stack(ensemble['g'].values)
    g_avg[i_m,i_r,i_d] = np.average(ensemble['g'].values, axis=0)

    bonds[i_m,i_r,i_d,0] = ensemble['psi'][ensemble['psi'] > 0.4].count()
    bonds[i_m,i_r,i_d,1] = ensemble['psi'][ensemble['psi'] < -0.3].count()

    #if(ensemble['psi'][ensemble['psi'] > 0.8].count() >41):
    #    print(ensemble['psi'])
    maxcluster_avg[i_m,i_r,i_d]=np.mean(ensemble['max_cluster'].values)

###################


#plt.figure(figsize=(7.0582, 7.0582*0.7))
x = np.linspace(0,len_g,601)/100

gp=np.zeros((len_g))
gn=np.zeros((len_g))
ga=g_avg[2,0,8]
pcount=0
ncount=0
for i_p, p in enumerate(psi[2,0,0]):
    if p == 1:
        print(p)
        gp += g[2,0,0,i_p]
        pcount += 1
    if p == -1:
        gn += g[2,0,0,i_p]
        ncount += 1
gp = gp/pcount
gn = gn/ncount

fig, ax = plt.subplots(3,1, figsize=(13,7),sharex=True,sharey=True)

ax[0].plot(x,ga, label=f'$\Delta=0.5$', c='#747270', lw=2 )
ax[0].legend(loc='best')
ax[1].plot(x,gp, label=f'$\Psi = 1$,$\Delta=0.2$', c='#f44b74',lw=2)
ax[1].legend(loc='best')
ax[2].plot(x,gn,label=f'$\Psi = -1$, $\Delta=0.2$', c='#0000ff',lw=2)
ax[2].legend(loc='best')

fig.text(0.02, 0.5, 'g(l)', va='center', rotation='vertical',fontsize=20, fontweight='bold')

plt.xlim((0.8,6.0))
plt.xlabel("l")
#plt.ylabel("g(l)")
#plt.tight_layout()
plt.savefig("../plots/gr_r_0.05.pdf")
plt.show()
