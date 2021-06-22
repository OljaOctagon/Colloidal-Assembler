import numpy as np
import pandas as pd 
import os 
import argparse
import networkx as nx
from glob import glob

def get_domains(arr):
    # get connections for indicated time 
    connections = np.column_stack((arr_i[:, 1], arr_i[:, 3]))

    # make a graph form connectins 
    G = nx.Graph()
    G.add_edges_from(connections)

    domains = list(nx.connected_components(G))

    return domains 


def get_xyz_file(time,nbonds,pbc,type_i,dir_i):

    time_offset =  time + 5000 
    file_id = int(np.ceil((time_offset/10)/2))

    print(file_id)
    poly_file = glob("{}/poly/PolymerPos*.T00{}.dat".format(dir_i, file_id))[0]
    linker_file = glob("{}/poly/LinkerPos.*.T00{}.dat".format(dir_ifile_id))[0]

    # read crosslinker connections
    lf = glob("{}/*link.dat".format(dir_i))[0]
    arr = pd.read_csv(lf, delim_whitespace=True).values
    # erase self bonds
    # arr = arr[arr[:, 1] != arr[:, 3]]

    dir_id = dir_i.split("/")[1]
    print("dir id", dir_id)
    save_dir = "/home/carina/git_repos/rhombi/utils/biogel/visual_inspection/{}".format(dir_id)

    if os.path.exists(save_dir):
    pass 
    else:
    os.mkdirs(save_dir)

    df = pd.read_csv(poly_file, 
        delim_whitespace=True, names=['time','x','y','z'])

    dg = pd.read_csv(linker_file,
        delim_whitespace=True, names=['time','pid', 'lid','x','y','z','var1','var2'])


    # select time in linker file 
    df_linker = dg[dg.time == time_offset ]

    # select time in pos file 
    df_pos = df[df.time == time_offset ]

    # add column with polymer id starting with 1 to Npolymers 
    m_id = np.array(range(0,len(df_pos)))
    df_pos['pid'] = np.floor(m_id/float(nbonds))+1 


    # initialize output dataframe 
    df_enriched = pd.DataFrame(columns=['a','b'])

    if type_i == 'poly_id':

    df_enriched = df_pos
    # monomer type is polymer id 
    monomer_type = 'M'+ df_enriched['pid'].astype('int').astype('str').values

    # denote crosslinkers 
    brr = ((df_linker.pid-1)*nbonds + (df_linker.lid-1)).values
    monomer_type[brr] = 'C'+monomer_type[brr]
    # add type to enriched dataframe 
    df_enriched['type'] = monomer_type 


    if type_i == 'cluster_size':

    # get connected domains 
    arr_i = arr[arr[:, 0] == time]
    domains = get_domains(arr_i)
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
    brr = ((df_linker.pid-1)*nbonds + (df_linker.lid-1)).values
    monomer_type[brr] = 'C'+monomer_type[brr]

    # add type to enriched dataframe 
    df_enriched['type'] = monomer_type 


    if type_i == 'cluster_id':
    # get connected domains
    arr_i = arr[arr[:, 0] == time]
    domains = get_domains(arr_i)
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
    brr = ((df_linker.pid-1)*nbonds + (df_linker.lid-1)).values
    monomer_type[brr] = 'C'+monomer_type[brr]

    # add type to enriched dataframe 
    df_enriched['type'] = monomer_type 


    if type_i == 'largest_cluster':
    # get connected domains 
    arr_i = arr[arr[:, 0] == time]
    domains = get_domains(arr_i)
    pid_dict = {}

    domain_lengths = np.array([ len(domain) for domain in domains ])
    max_domain = np.max(domain_lengths)
    print(max_domain)

    # get dictionary of pid and length of domain pid is element of 
    for ni, domain in enumerate(domains):
        len_domain = len(domain)
        print(len_domain)
        di = "0"
        if len_domain == max_domain:
            di = str(ni+1) 

        for elem in domain:
            pid_dict[elem] = di

    df_pid = pd.DataFrame(pid_dict.items(),columns=['pid','cluster_id'])

    # join position frame and domain length frame 
    df_enriched = df_pos.merge(df_pid, on='pid',how='left').fillna(0)

    # denote monomer type according to its cluster size 
    monomer_type = 'M'+ df_enriched['cluster_id'].astype('int').astype('str').values

    # denote crosslinkers 
    brr = ((df_linker.pid-1)*nbonds + (df_linker.lid-1)).values
    monomer_type[brr] = 'C'+monomer_type[brr]

    # add type to enriched dataframe 
    df_enriched['type'] = monomer_type 


    # periodic boundaries 
    if pbc=="on":
    Lb=30
    df_enriched['x'] = df_enriched['x'] - Lb*np.rint(df_enriched['x']/Lb)
    df_enriched['y'] = df_enriched['y'] - Lb*np.rint(df_enriched['y']/Lb)
    df_enriched['z'] = df_enriched['z'] - Lb*np.rint(df_enriched['z']/Lb)


    brr = df_enriched[['type', 'x','y','z']].values
    with open("{}/biogel_cluster_color_time_{}_pbc_{}_type_{}.xyz".format(save_dir, time, pbc, type_i),'w') as f:
    for [mtype, x,y,z,] in brr:
        f.write("{}   {}   {}   {}\n".format(mtype,x,y,z))


if __name__ == "__main__": 

names=['high_flexibility_kb0','short_chains','plink_variation','std_conditions']

pbc=["on","off"]
types=['poly_id', 'cluster_id', 'cluster_size', 'largest_cluster']

time==1000 
for name in names:
    dirs=glob("{}/*/pdf1".format(name))
    for dir_i in dirs:
        print("dir_i", dir_i)
        nbonds=30
        if name == "short_chains":
            nbonds = 10 

        for pbc_i in pbc:
            for type_i in types:
                get_xyz_file(time,nbonds,pbc_i,type_i,dir_i)


