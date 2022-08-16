import numpy as np
import pandas as pd 
import os 
from glob import glob 
nbonds=30

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
    for time in Ntimes:
        if time==1000:

            dh = df[df.time == time]
            ds = dg[dg.time == time]
    
            arr = ((ds.pid-1)*nbonds + (ds.lid-1)).values
            monomer_type = np.array((['M']*len(dh)))
            monomer_type[arr] = 'C'
            dh['type'] = monomer_type 

            brr = dh[['type', 'x','y','z']].values
            with open("{}/biogel_{}.xyz".format(save_dir, time),'w') as f:
                for [mtype, x,y,z,] in brr:
                    f.write("{}   {}   {}   {}\n".format(mtype,x,y,z))

