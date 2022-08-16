from glob import glob
from pymol import cmd

cmd.load("/Users/ada/Documents/Code_Development_2020/biogel/results/snapshot_analysis/short_chains_e5.xyz", "e5")

cmd.show_as("spheres")
cmd.alter('elem M', "vdw=0.5")
cmd.alter('elem C', "vdw=0.7")
cmd.alter('elem C', "color=30")
cmd.alter('elem M', "color=11")
cmd.bg_color("white")
cmd.set("orthoscopic")
cmd.rebuild()
