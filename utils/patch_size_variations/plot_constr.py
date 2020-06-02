import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import matplotlib as mpl
import matplotlib.style as style

#matplotlib.use('pgf')
#matplotlib.rcParams.update({
#        "pgf.texsystem": "pdflatex",
#        'font.family'  : 'serif',
#        'text.usetex'  : True,
#        'pgf.rcfonts'  :False,
#        })

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

savepath = "../plots/"

df = pd.read_pickle("../analysis_data/constr_mp.pkl")

f = ['#state',
	'#cp',
	'mu', 
	'delta',
	'r_patch', 
	'percolation', 
	'#p', 
	'#np', 
	'psi', 
	'max_cluster', 
	'max_p_domain', 
	'max_np_domain', 
	'polybonds/bonds', 
	'max_bonds',
	'theta_diff_list',
	'g',
	'cluster_sizes',
	'p_domain_sizes',
	'np_domain_sizes'
	] 

psi_mu_d_r = df.groupby(['mu', 'delta', 'r_patch'])['psi'].mean().to_numpy()
psi_std =  df.groupby(['mu', 'delta', 'r_patch'])['psi'].std().to_numpy()
deltas = df.delta.unique()
mus = df.mu.unique()
g_mu_d_r = df.groupby(['mu', 'delta', 'r_patch'])['g'].apply(np.mean).to_numpy()

fcounter=0
pcounter=0
npcounter=0

for i in range(len(psi_mu_d_r)):
	print(psi_mu_d_r[i])

x = np.linspace(0,np.shape(g_mu_d_r[0]),601)/100

for i in range(len(g_mu_d_r)):
	if(fcounter == 0 and np.abs(psi_mu_d_r[i])<0.1):
		plt.plot(x,g_mu_d_r[i], linewidth=0.5,label=f'$\Psi$ = {psi_mu_d_r[i]}')
		fcounter=1
		print("fluid", i)
	if(npcounter == 0 and psi_mu_d_r[i]<-0.7):
		plt.plot(x,g_mu_d_r[i], linewidth=0.5,label=f'$\Psi$ = {psi_mu_d_r[i]}')
		npcounter=1
		print("np", i)
	if(pcounter == 0 and psi_mu_d_r[i]>0.05):
		plt.plot(x,g_mu_d_r[i], linewidth=0.5,label=f'$\Psi$ = {psi_mu_d_r[i]}')
		pcounter=1
		print("p", i, psi_mu_d_r[i])

plt.xlabel("r")
plt.ylabel("g(r)")
plt.legend()
plt.savefig(savepath+"constr_g.pdf")
plt.clf()

for mu in range(len(mus)):
	plt.errorbar(deltas,psi_mu_d_r[mu*len(deltas):(mu+1)*len(deltas)],
                 yerr=psi_std[mu*len(deltas):(mu+1)*len(deltas)],
                 label = f'$\mu$ = {mus[mu]}',
                 linewidth=.5, capsize=1,
                 capthick=0.5)

plt.legend()
plt.ylim(-1,1)
plt.xlabel("$\Delta$")
plt.ylabel("$\Psi_{\mathrm{rand}}$")
plt.savefig(savepath+"constr_psis.pdf")
plt.clf()

anomaly = df[(df.mu == 0.2) & (df.delta == 0.35)]['psi'].mean()

print(anomaly)
#######errorbars
#####evtl split 1 / -1





#print(find_row("#state"))

# for i,row in df.iterrows():
# 	print(row[:])
# 	print()


# print(df)
# print(df.groupby(['mu','delta','r_patch']).mean())

# print(df.groupby(['mu','delta','r_patch']).mean()['Psi'][:])

# print(df.groupby(['mu','delta','r_patch'])['np_domain_sizes'].apply(lambda g : len(g)))


# print(df.groupby(['mu','delta','r_patch']).sum())

# print(df[df.delta == '0.45']['percolation'])

# sum = 0
# for i in df[df.delta == '0.45']['Psi']:
# 	print(sum)
# 	print(i)
# 	sum+=i
# 	print(sum)
# 	print()
# print(df.apply(np.sum, axis = 0, columns = 'cluster_sizes'))

# for i,row in df[df.delta == '0.45'].iterrows():
# 	print(row['delta'])
# 	print()
# df.groupby
# df.reindex
# df.apply
# df.to_frame











