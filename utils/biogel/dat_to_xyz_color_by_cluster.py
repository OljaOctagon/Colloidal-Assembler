import numpy as np
import pandas as pd 
import os 
from glob import glob 
import argparse

nbonds=30

parser = argparse.ArgumentParser()
parser.add_argument('-lf', type='str', help='linker file destination')
parser.add_argument('-nbonds', type='str', help='number of monomers per polymer')
args = parser.parse_args()


poly_files = glob("PolymerPos*.dat")

save_dir = "xyz-files"
os.mkdir(save_dir)

for file in poly_files:
    df = pd.read_csv(file, 
            delim_whitespace=True, names=['time','x','y','z'])

    ids  = file.split(".")
    linker_id = ids[2]
    run_id = ids[1]

    dg = pd.read_csv("LinkerPos.{}.{}.dat".format(run_id, linker_id), 
            delim_whitespace=True, names=['time','pid', 'lid','x','y','z','var1','var2'])
    
    Ntimes = np.sort(df.time.unique())
    
    arr = pd.read_csv(linker_file, delim_whitespace=True).values
    arr = arr[arr[:, 1] != arr[:, 3]]
    
    for time in Ntimes:
        if time==1000:

            arr_i = arr[arr[:, 0] == time_i]
            connections = np.column_stack((arr_i[:, 1], arr_i[:, 3]))
            G = nx.Graph()
            G.add_edges_from(connections)

            domains = list(nx.connected_components(G))
            domain_lengths = np.array([ len(domain) for domain in domains ])
            
            pid_dict = {}
            for domain in domains:
                len_domain = len(domain)
                for elem in domains:
                    pid_dict[elem] = len_domain

            
            df_pid = pd.DataFrame(pid_dict.items(),columns=['pid','len_pid'])


            dh = df[df.time == time]
            ds = dg[dg.time == time]
            m_id =int(np.linspace(0,len(dh),len(dh))
            dh['pid'] = np.floor(m_id/nbonds)

            df_enriched = dh.merge(df_pid, on='pid')
            monomer_type = 'M'+str(df_enriched['len_pid'].values)
            
            arr = ((ds.pid-1)*nbonds + (ds.lid-1)).values
            monomer_type[arr] = 'C'
            df_enriched['type'] = monomer_type 


            brr = df_enriched[['type', 'x','y','z']].values
            with open("{}/biogel_{}.xyz".format(save_dir, time),'w') as f:
                for [mtype, x,y,z,] in brr:
                    f.write("{}   {}   {}   {}\n".format(mtype,x,y,z))


