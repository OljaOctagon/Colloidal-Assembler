from glob import glob
from pymol import cmd

id='00384'
name='short_chains'

bdir="/home/carina/Documents/biogels/visual_inspection/"
subdir="{}/{}/pdf1/poly/xyz-files/".format(name,id)
lst=glob("{}{}/*xyz".format(bdir,subdir))
lst.sort()
for f in lst: 
    cmd.load(f,"e0")

cmd.show_as("spheres")
cmd.alter('elem M', "vdw=0.5")
cmd.alter('elem C', "vdw=0.7")
cmd.alter('elem C', "color=30")
cmd.alter('elem M', "color=11")
cmd.bg_color("white")
cmd.set("orthoscopic")
cmd.rebuild()
