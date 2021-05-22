import numpy as np
import pandas as pd 
import os 
from glob import glob 
nbonds=30

odir="/home/zoettl/Results/data/EllPoly/2019b/visual_inspection/std_conditions/00111/poly/pos-files"

poly_files = glob("{}/PolymerPos.00001.T*.dat".format(odir))

save_dir = "/home/zoettl/Results/data/EllPoly/2019b/visual_inspection/std_conditions/00111/poly/xyz-files"
os.mkdir(save_dir)

for file in poly_files:
    df = pd.read_csv(file, 
            delim_whitespace=True, names=['time','x','y','z'])

    linker_id = file.split(".")[2]
    dg = pd.read_csv("{}/LinkerPos.00001.{}.dat".format(odir, linker_id), 
            delim_whitespace=True, names=['time','pid', 'lid','x','y','z','var1','var2'])
    

    Ntimes = np.sort(df.time.unique())
    for time in Ntimes:
        if time%100 == 0:

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

