from glob import glob
from pymol import cmd
import os 
import argparse
import matplotlib.pylab as plt


#parser = argparse.ArgumentParser()
#parser.add_argument("-f")
#args = parser.parse_args()

cwd=os.getcwd()
file_path="poly/xyz-files/color_pid_biogel_9900.0.xyz"
fs="00369"

ncolors=10
cmap = plt.cm.get_cmap('viridis', ncolors)

for i in range(ncolors): cmd.set_color("mycol{}".format(i), list(cmap(i))[:3])

cmd.load("{}/{}/{}".format(cwd,fs,file_path), fs)

cmd.show_as("spheres")
cmd.alter('elem C', "vdw=0.6")
#cmd.alter('elem C', "color=30")
cmd.color("red", "elem C")
n_poly = 340 
for i in range(340):

	color_i= i%ncolors
	cmd.color("mycol{}".format(color_i), "elem M{}".format(i) )
	#cmd.alter('elem M'+str(i), "color={}".format(color_i))
	cmd.alter('elem M'+str(i), "vdw=0.5")


cmd.bg_color("white")
cmd.set("orthoscopic")
cmd.rebuild()
