from glob import glob
from pymol import cmd
import pandas as pd
import seaborn as sns
from natsort import natsorted

directory="/Users/ada/Documents/git_repos/phd/rhombi/utils/biogel/visual_inspection"
#names=['high_flexibility_kb0','short_chains','plink_variation','std_conditions']
names=['std_conditions']
print("hi")
for name in names:
    print(name)
    filelist=glob("{}/{}/*/*.xyz".format(directory, name))
    for i,fi in enumerate(filelist):
        cmd.load(fi,"e"+str(i))
        df = pd.read_csv(fi, names=['mtype','x','y','z'], delim_whitespace=True)
        unique_val = df.mtype.unique()
        cmd.show_as("spheres")
        cmd.alter('all', "vdw=0.4")

        num_shades = len(unique_val)

        unique_val=natsorted(unique_val)
        color_list = sns.color_palette("viridis", num_shades)

        color_dict=dict(zip(unique_val, color_list))

        for i,mi in enumerate(unique_val):
            rgb_val = color_dict[mi]
            cmd.set_color( str(mi), rgb_val)
            cindex = cmd.get_color_index(str(mi))

            cmd.alter('elem '+str(mi), "color="+str(cindex))
            cmd.alter('elem '+str(mi), "vdw=0.5")

            cmd.alter('elem C'+str(mi), "color = 8")
            cmd.alter('elem C'+str(mi), "vdw=0.5")

        cmd.bg_color("white")
        cmd.set("orthoscopic")
        cmd.rebuild()
        
        savefile=fi.split(".")[0]+".png"
        cmd.png(savefile, dpi=300)
        cmd.delete("e"+str(i))

