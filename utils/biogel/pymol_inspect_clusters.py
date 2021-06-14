from glob import glob
from pymol import cmd
import pandas as pd

id='00360'
name='visual_inspection/high_flexibility'

directory=
bfile = "{}/{}/{}/biogel_cluster_color_1000.0.xyz".format(directory, name,id)
cmd.load(bfile,"e0")

df = pd.read_csv(bfile, names=['mtype','x','y','z'], delim_whitespace=True)
unique_val = df.mtype.unique()

cmd.show_as("spheres")
cmd.alter('elem C', "vdw=0.7")
cmd.alter('elem C', "color=30")
for i,mi in enumerate(unique(val)):
    cmd.alter('elem'+str(mi), "color="+str(i)
    cmd.alter('elem'+str(mi), "vdw=0.5")

cmd.bg_color("white")
cmd.set("orthoscopic")
cmd.rebuild()
