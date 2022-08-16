import pandas as pd 
import numpy as np 
import argparse 
import warnings
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser()
parser.add_argument('-f')
#parser.add_argument('-l')

args = parser.parse_args()

df = pd.read_csv(args.f, delim_whitespace=True, names=["ptype",'x','y','z'])

Npoly=340
l_poly=30

a = np.arange(1, Npoly + 1)
b = np.arange(1, l_poly + 1)
f_to_number = lambda tup: (tup[0] - 1) * l_poly + (tup[1] - 1)
arr_combinations = np.array([np.meshgrid(a, b)]).T.reshape(-1, 2)
tuple_combinations = [(i, j) for [i, j] in arr_combinations]
map_to_number = list(map(f_to_number, tuple_combinations))
dict_to_number = dict(zip(tuple_combinations, map_to_number))
dict_to_tuple = dict(zip(map_to_number,tuple_combinations))


import networkx as nx 
df['index'] = df.reset_index().index

df_cross = df[df.ptype=='C']
arr_index = df_cross.index.values
arr_pos = df_cross[['x','y','z']].values 

linker_arr = []

L=np.ones(3)*30
sigma=1.5

for inum,index_i in enumerate(arr_index):
	for jnum, index_j in enumerate(arr_index):
		if inum<jnum:
			dist = arr_pos[inum] - arr_pos[jnum]
			dist = dist - L*np.rint(dist/L)
			dnorm = np.linalg.norm(dist)
			if dnorm<sigma:
				linker_arr.append([index_i,index_j])


G = nx.Graph()
G.add_edges_from(linker_arr)
domains = list(nx.connected_components(G))
domain_lengths = np.array([ len(domain) for domain in domains ])
d_id = np.argmax(domain_lengths)
particles_max_domain = np.array(list(domains[d_id]))
particles_max_domain = np.unique(particles_max_domain)
pid_max_domain = [ dict_to_tuple[num][0]-1 for num in particles_max_domain ]
pid_max_domain = np.unique(pid_max_domain)
#######################

poly_id=np.repeat(range(Npoly),l_poly)
df['poly_id'] = poly_id 

dg = df[df.poly_id.isin(pid_max_domain)]


print(len(pid_max_domain))
print(len(particles_max_domain))

#### 

length=len(dg)
rep = int(length/l_poly)
arr=np.repeat(np.array(range(0,rep)).astype(str), l_poly)
dg.ptype = dg.ptype + arr 
dg.ptype[dg.ptype.str.match('C')] = 'C'

#### assemble connected periodic image of crosslinker cluster 

def center_of_mass_crosslinker(arr,L):

	theta = (arr/L)*2*np.pi
	cos_theta = np.cos(theta)
	sin_theta = np.sin(theta)
	
	av_cos_theta = np.mean(cos_theta)
	av_sin_theta = np.mean(sin_theta)
	

	theta_bar = np.arctan2(-av_cos_theta, -av_sin_theta) + np.pi 

	center_mass = L*(theta_bar/(2*np.pi))

	return center_mass 

arr = dg[dg.index.isin(particles_max_domain)][['x','y','z']].values
center_mass = center_of_mass_crosslinker(arr,L)

list_poly = []
for pi in particles_max_domain:

	pl_id = dg[ dg.index == pi].poly_id.values[0]
	
	if pl_id not in list_poly: 
		list_poly.append(pl_id)
		i_cross_monomer = pi%l_poly

		cpos = dg[dg.index == pi][['x','y','z']].values 
		cdist = cpos - center_mass 
		cdist = cdist - L*np.rint(cdist/L)
		pos = np.zeros((30,3))
		pos[i_cross_monomer] = center_mass + cdist

		di = dg[dg.poly_id == pl_id]
		arr = di[['x','y','z']].values
		arr_index  = di.index.values	

		for i in range(i_cross_monomer):
			dist = arr[i_cross_monomer-i-1] - pos[i_cross_monomer-i]
			dist = dist - L*np.rint(dist/L)
			pos[i_cross_monomer-i-1] = pos[i_cross_monomer-i] + dist 

		for i in range(i_cross_monomer,l_poly-1):
			dist = arr[i+1] - pos[i]
			dist = dist - L*np.rint(dist/L)
			pos[i+1] = pos[i]+dist
			
		dg['x'][dg.poly_id == pl_id] =  pos[:,0]
		dg['y'][dg.poly_id == pl_id] =  pos[:,1]
		dg['z'][dg.poly_id == pl_id] =  pos[:,2]
	
'''
a=np.array([-30,0,30])
axis_combinations = np.array([np.meshgrid(a,a,a)]).T.reshape(-1,3)

dg_copy = dg.copy()

for ax in axis_combinations:
	dg1 = dg_copy.copy()
	dg1['x'] = dg1['x'] + ax[0]
	dg1['y'] = dg1['y'] + ax[1]
	dg1['z'] = dg1['z'] + ax[2]

	dg = dg.append(dg1)
'''

brr = dg[['ptype', 'x','y','z']].values
with open("cluster_pid_{}".format(args.f),'w') as f:
	for [mtype, x,y,z,] in brr:
		f.write("{}   {}   {}   {}\n".format(mtype,x,y,z))
