import numpy as np
import pandas as pd 
import os 
import argparse
import networkx as nx
from glob import glob


def get_domains(arr):
    # get connections for indicated time 
    arr_i = arr[arr[:, 0] == args.time]
    connections = np.column_stack((arr_i[:, 1], arr_i[:, 3]))

    # make a graph form connectins 
    G = nx.Graph()
    G.add_edges_from(connections)

    domains = list(nx.connected_components(G))

return domains 


parser = argparse.ArgumentParser()
parser.add_argument('-lf', type=str, help='linker file destination')
parser.add_argument('-nbonds', type=int, help='number of monomers per polymer')
parser.add_argument("-time", type=int, help='time stamp')
parser.add_argument("-type", type=str, choices=['poly_id', 'cluster_id', 'cluster_size', 'largest_cluster'])
parser.add_argument('-pbc', type=str)
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


# select time in linker file 
df_linker = dg[dg.time == time_offset ]

# select time in pos file 
df_pos = df[df.time == time_offset ]
m_id =np.linspace(0,len(df_pos),len(df_pos))

# add column with polymer id starting with 1 to Npolymers 
df_pos['pid'] = np.floor(m_id/float(args.nbonds))+1 

# read crosslinker connections
arr = pd.read_csv(args.lf, delim_whitespace=True).values

# erase self bonds
# arr = arr[arr[:, 1] != arr[:, 3]]

df_enriched = pd.DataFrame(columns=['a','b'])


if args.type == 'poly_id':
    
    df_enriched = df_pos
    # monomer type is polymer id 
    monomer_type = 'M'+ df_enriched['pid'].astype('int').astype('str').values

    # denote crosslinkers 
    brr = ((df_linker.pid-1)*args.nbonds + (df_linker.lid-1)).values
    monomer_type[brr] = 'C'+monomer_type[brr]
    # add type to enriched dataframe 
    df_enriched['type'] = monomer_type 


if args.type == 'cluster_size':

    # get connected domains 
    domains = get_domains(arr)
    pid_dict = {}
    
    # get dictionary of pid and length of domain pid is element of 
    for domain in domains:
        len_domain = len(domain)
        for elem in domain:
            pid_dict[elem] = len_domain

    df_pid = pd.DataFrame(pid_dict.items(),columns=['pid','len_pid'])

    # join position frame and domain length frame 
    df_enriched = df_pos.merge(df_pid, on='pid',how='left').fillna(0)
    
    # denote monomer type according to its cluster size 
    monomer_type = 'M'+ df_enriched['len_pid'].astype('int').astype('str').values

    # denote crosslinkers 
    brr = ((df_linker.pid-1)*args.nbonds + (df_linker.lid-1)).values
    monomer_type[brr] = 'C'+monomer_type[brr]
    
    # add type to enriched dataframe 
    df_enriched['type'] = monomer_type 


if args.type == 'cluster_id':
    # get connected domains 
    domains = get_domains(arr)
    pid_dict = {}
    
    # get dictionary of pid and length of domain pid is element of 
    for di, domain in enumerate(domains):
        for elem in domain:
            pid_dict[elem] = di

    df_pid = pd.DataFrame(pid_dict.items(),columns=['pid','cluster_id'])

    # join position frame and domain length frame 
    df_enriched = df_pos.merge(df_pid, on='pid',how='left').fillna(0)
    
    # denote monomer type according to its cluster size 
    monomer_type = 'M'+ df_enriched['cluster_id'].astype('int').astype('str').values

    # denote crosslinkers 
    brr = ((df_linker.pid-1)*args.nbonds + (df_linker.lid-1)).values
    monomer_type[brr] = 'C'+monomer_type[brr]
    
    # add type to enriched dataframe 
    df_enriched['type'] = monomer_type 


if args.type == 'largest_cluster:'
    # get connected domains 
    domains = get_domains(arr)
    pid_dict = {}
    
    domain_lengths = np.array([ len(domain) for domain in domains ])
    max_domain = np.max(domain_lengths)

    # get dictionary of pid and length of domain pid is element of 
    for ni, domain in enumerate(domains):
        len_domain = len(domain)
        di = "0"
        if len_domain == max_domain:
            di == str(ni+1) 

        for elem in domain:
            pid_dict[elem] = di

    df_pid = pd.DataFrame(pid_dict.items(),columns=['pid','cluster_id'])

    # join position frame and domain length frame 
    df_enriched = df_pos.merge(df_pid, on='pid',how='left').fillna(0)
    
    # denote monomer type according to its cluster size 
    monomer_type = 'M'+ df_enriched['cluster_id'].astype('int').astype('str').values

    # denote crosslinkers 
    brr = ((df_linker.pid-1)*args.nbonds + (df_linker.lid-1)).values
    monomer_type[brr] = 'C'+monomer_type[brr]
    
    # add type to enriched dataframe 
    df_enriched['type'] = monomer_type 


# periodic boundaries 
if args.pbc=="on":
    Lb=30
    df_enriched['x'] = df_enriched['x'] - Lb*np.rint(df_enriched['x']/Lb)
    df_enriched['y'] = df_enriched['y'] - Lb*np.rint(df_enriched['y']/Lb)
    df_enriched['z'] = df_enriched['z'] - Lb*np.rint(df_enriched['z']/Lb)


brr = df_enriched[['type', 'x','y','z']].values
with open("{}/biogel_cluster_color_{}_{}.xyz".format(save_dir, args.time, args.type),'w') as f:
    for [mtype, x,y,z,] in brr:
        f.write("{}   {}   {}   {}\n".format(mtype,x,y,z))


