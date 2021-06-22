import numpy as np
import pandas as pd 
import os 
import argparse
import networkx as nx
from glob import glob

parser = argparse.ArgumentParser()
parser.add_argument('-lf', type=str, help='linker file destination')
parser.add_argument('-nbonds', type=int, help='number of monomers per polymer')
parser.add_argument("-time", type=int, help='time stamp')
args = parser.parse_args()

time_offset =  args.time + 5000 
file_id = int(np.ceil((time_offset/10)/2))

print(file_id)
poly_file = glob("PolymerPos*.T00{}.dat".format(file_id))[0]
linker_file = glob("LinkerPos.*.T00{}.dat".format(file_id))[0]

save_dir = "xyz-files-cluster-size"

if os.path.exists(save_dir):
    pass 
else:
    os.mkdir(save_dir)

df = pd.read_csv(poly_file, 
        delim_whitespace=True, names=['time','x','y','z'])

dg = pd.read_csv(linker_file,
        delim_whitespace=True, names=['time','pid', 'lid','x','y','z','var1','var2'])


arr = pd.read_csv(args.lf, delim_whitespace=True).values
arr = arr[arr[:, 1] != arr[:, 3]]

print("time 1000")
arr_i = arr[arr[:, 0] == args.time]
connections = np.column_stack((arr_i[:, 1], arr_i[:, 3]))

print("graph")
G = nx.Graph()
G.add_edges_from(connections)

domains = list(nx.connected_components(G))
domain_lengths = np.array([ len(domain) for domain in domains ])

print("domain element")
pid_dict = {}
for domain in domains:
    len_domain = len(domain)
    for elem in domain:
        pid_dict[elem] = len_domain


df_pid = pd.DataFrame(pid_dict.items(),columns=['pid','len_pid'])
print(df_pid.pid.min(), df_pid.pid.max())

df_pos = df[df.time == time_offset ]
print(len(df_pos),"len")

df_linker = dg[dg.time == time_offset ]
m_id =np.linspace(0,len(df_pos),len(df_pos))
df_pos['pid'] = np.floor(m_id/float(args.nbonds))+1 
print("len", len(df_pos))

df_enriched = df_pos.merge(df_pid, on='pid',how='left').fillna(0)

print("len enrich", len(df_enriched))
monomer_type = 'M'+ df_enriched['len_pid'].astype('int').astype('str').values

arr = ((df_linker.pid-1)*args.nbonds + (df_linker.lid-1)).values
monomer_type[arr] = 'C'+monomer_type[arr]
df_enriched['type'] = monomer_type 


brr = df_enriched[['type', 'x','y','z']].values
with open("{}/biogel_cluster_color_{}.xyz".format(save_dir, args.time),'w') as f:
    for [mtype, x,y,z,] in brr:
        f.write("{}   {}   {}   {}\n".format(mtype,x,y,z))


