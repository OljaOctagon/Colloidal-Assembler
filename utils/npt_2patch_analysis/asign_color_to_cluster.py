import numpy as np
from itertools import cycle 


color_array = np.zeros(1000)
current_color = 0

time_id = 0

with open("All_Clusters_info.dat",'r') as f:
    with open("cluster_colors.dat",'w') as g:
        for line in f.readlines():
            line = line.split(" ")
            if len(line) == 3:
                current_color = np.random.randint(0,11) 
                if time_id != line[0]:
                    if time_id != 0:
                        for i in range(1000):
                            g.write("{}\n".format(color_array[i]))

                    time_id = line[0]
                    g.write("1000\n")
                    g.write("Cluster colors of frame {}\n".format(time_id))

            elif len(line) == 1:
                line = line[0].split("\n")
                color_array[int(line[0])] = current_color 
            
