from glob import glob
from pymol import cmd
import os 
import argparse
import cmasher as cmr
import matplotlib.pyplot as plt 

cwd=os.getcwd()
file_path="poly/xyz-files/color_pid_biogel_9900.0.xyz"
#file_path="poly/xyz-files/cluster_pid_biogel_9900.0.xyz"

scenes = ["00111","00112","00113",'00114',"00115"]

for fs in scenes: 
	n_monomers=10200
	lpoly=30
	ncolors=10
	n_poly = 340 
	
	cmap = plt.get_cmap('cmr.lavender', ncolors)
	for i in range(ncolors): cmd.set_color("mycol{}".format(i), list(cmap(i))[:3])
	cmd.load("{}/{}/{}".format(cwd,fs,file_path), fs)
	cmd.unbond("all", "all")

	# every polymer has different color 
	for i in range(n_poly):
		color_i= i%ncolors
		cmd.color("mycol{}".format(color_i), "elem M{}".format(i))
		cmd.alter('elem M'+str(i), "vdw=0.1")
	
	# bond monomers for every polymer 
	for imono in range(1,n_monomers):
		a=int((imono-1)//lpoly)
		b=int(imono//lpoly)
		if a==b:
			bonds = cmd.find_pairs('index {}'.format(imono-1), 'index {}'.format(imono), cutoff=1.5)
			for pair in bonds:
				cmd.bond('index {}'.format(imono-1), 'index {}'.format(imono))

	# draw polymers as sticks 
	cmd.set("stick_radius", 0.15)
	cmd.show("sticks")
		
	# draw box 	
	bonds = cmd.find_pairs('n. b', 'n. b', cutoff=31)
	for pair in bonds:
		pair1=pair[0][1]
		pair2=pair[1][1]

		if pair[0][1]>pair[1][1]:
			pair1 = pair[1][1]
			pair2 = pair[0][1]

		cmd.bond('index {}'.format(pair1), 'index {}'.format(pair2))
		cmd.set("stick_radius", 0.2, "elem B")

	cmd.color("black", "elem B")
	cmd.alter("elem B", "vdw=0.1")

	# draw corsslinkers as speres 
	cmd.show("spheres", "elem C")
	cmd.alter('elem C', "vdw=0.5")
	cmd.set_color("ccolor",[105/255, 62/255, 227/255])	
	cmd.color("red", "elem C")

	# rotate and zoom view 
	cmd.rotate("[0,-1,0]", "-30", "all")
	cmd.rotate("[1,0,0]", "20", "all")
	cmd.zoom( "all", "3")

	# set background color and orthoscopic view 
	cmd.bg_color("white")
	cmd.set("orthoscopic")

	cmd.rebuild()

	# ray tracing 
	cmd.set("ray_texture", "1")
	cmd.set("ray_opaque_background", 'on')
	cmd.set("ray_trace_gain", '8')
	cmd.set("depth_cue", '0')
	util.ray_shadows('soft')

	# save output 
	cmd.png("{}_lavender.png".format(fs),width=3200,height=3200, dpi=1000, ray=1)
	cmd.delete("all")

cmd.quit()
