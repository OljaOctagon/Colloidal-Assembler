from glob import glob
from pymol import cmd
import os 
import argparse

#parser = argparse.ArgumentParser()
#parser.add_argument("-f")
#args = parser.parse_args()

cwd=os.getcwd()
file_path="poly/xyz-files/biogel_9900.0.xyz"
fs="00364"
cmd.load("{}/{}/{}".format(cwd,fs,file_path), fs)

cmd.show_as("spheres")
cmd.alter('elem M', "vdw=0.5")
cmd.alter('elem C', "vdw=0.7")
cmd.alter('elem C', "color=30")
cmd.alter('elem M', "color=11")
cmd.bg_color("white")
cmd.set("orthoscopic")
cmd.rebuild()
