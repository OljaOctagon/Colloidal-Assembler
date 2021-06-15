import numpy as np
import pandas as pd 
import os 
from glob import glob 
import argparse
import networkx as nx


parser = argparse.ArgumentParser()
parser.add_argument('-lf', type=str, help='linker file destination')
parser.add_argument('-nbonds', type=int, help='number of monomers per polymer')
args = parser.parse_args()

poly_files = glob("PolymerPos*.dat")

save_dir = "xyz-files-cluster-size"

if os.path.exists(save_dir):
    pass 
else:
    os.mkdir(save_dir)

print ("enter loop")
for file in poly_files:
    df = pd.read_csv(file, 
            delim_whitespace=True, names=['time','x','y','z'])

    ids  = file.split(".")
    linker_id = ids[2]
    run_id = ids[1]

    dg = pd.read_csv("LinkerPos.{}.{}.dat".format(run_id, linker_id), 
            delim_whitespace=True, names=['time','pid', 'lid','x','y','z','var1','var2'])
    
    Ntimes = np.sort(df.time.unique())
    
    for time in Ntimes:
        if time==1000:

            arr = pd.read_csv(args.lf, delim_whitespace=True).values
            arr = arr[arr[:, 1] != arr[:, 3]]

            print("time 1000")
            arr_i = arr[arr[:, 0] == 2*time+5000]
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

            df_pos = df[df.time == time]
            print(len(df_pos),"len")

            df_linker = dg[dg.time == time]
            m_id =np.linspace(0,len(df_pos),len(df_pos))
            df_pos['pid'] = np.floor(m_id/float(args.nbonds))+1 
            print("len", len(df_pos))

            df_enriched = df_pos.merge(df_pid, on='pid',how='left').fillna(0)

            print("len enrich", len(df_enriched))
            monomer_type = 'M'+ df_enriched['len_pid'].astype('int').astype('str').values
            
            arr = ((df_linker.pid-1)*args.nbonds + (df_linker.lid-1)).values
            monomer_type[arr] = 'C'
            df_enriched['type'] = monomer_type 


            brr = df_enriched[['type', 'x','y','z']].values
            with open("{}/biogel_cluster_color_{}.xyz".format(save_dir, time),'w') as f:
                for [mtype, x,y,z,] in brr:
                    f.write("{}   {}   {}   {}\n".format(mtype,x,y,z))


