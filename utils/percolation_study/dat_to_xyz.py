import numpy as np
import pandas as pd 

file='PolymerPos.00001.T00001.dat'
df = pd.read_csv(file, 
        delim_whitespace=True, names=['time','x','y','z'])

dg = pd.read_csv("LinkerPos.00001.T00001.dat", 
        delim_whitespace=True, names=['time','pid', 'lid','x','y','z','var1','var2'])
nbonds=30
arr = ((dg.pid-1)*nbonds + (dg.lid-1)).values
monomer_type = np.array((['N']*len(df)))
monomer_type[arr] = 'C'

df['type'] = monomer_type 

Ntimes = np.sort(df.time.unique())
print(Ntimes)
for time in Ntimes:
    arr =  df[df.time==time][['type', 'x','y','z']].values
    with open("biogel_{}.xyz".format(time),'w') as f:
        for [mtype, x,y,z,] in arr:
            f.write("{}   {}   {}   {}\n".format(mtype,x,y,z))

