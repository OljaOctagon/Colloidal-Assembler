import pandas as pd 
import numpy as np 

import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('-f')

args = parser.parse_args()

df = pd.read_csv(args.f, delim_whitespace=True, names=["ptype",'x','y','z'])

length=len(df)
l_poly=30
rep = int(length/l_poly)
arr=np.repeat(np.array(range(0,rep)).astype(str), l_poly)
df.ptype = df.ptype + arr 
df.ptype[df.ptype.str.match('C')] = 'C'

Lh=15
box_arr = np.array([
	[-Lh,-Lh,-Lh],
	[Lh,-Lh,-Lh],
	[Lh,Lh,-Lh],
	[-Lh,Lh,-Lh],
	[-Lh,-Lh,Lh],
	[Lh,-Lh,Lh],
	[Lh,Lh,Lh],
	[-Lh,Lh,Lh]]
	)
df_box = pd.Dataframe(data=box_arr, columns=['x','y','z'])
df_box['ptype'] = ['B']*8

print(df_box)

df.append(df_box)

brr = df[['ptype', 'x','y','z']].values
with open("color_pid_{}".format(args.f),'w') as f:
	for [mtype, x,y,z,] in brr:
		f.write("{}   {}   {}   {}\n".format(mtype,x,y,z))


