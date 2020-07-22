import numpy as np
import numpy as np
import matplotlib.pyplot as plt 
import glob 

systems=[dma_as1,dmo_as1,dmo_s1,dmo_s2]
deltas=[0.2,0,3,0.4,0.5]

sys_dict = {"dma_as1":"double_manta_asymm_1",
            "dmo_as1":"double_mouse_asymm_1",
            "dmo_s1":"double_mouse_symm_1",
            "dmo_s2":"double_mouse_symm_2"}

fig,ax = plot.subplots(ncol=len(deltas),nrow=len(systems),figsize=(12,12))

for si, system in enumerate(system):
    for di, delta in enumerate(deltas):

        double_mouse_symm_1_phi_0.1_delta_0.2_temp_0.05/
        ndir="{}/{}_phi_*_delta_{}_temp_*".format(system,sys_dict[system],delta)
        files = glob.glob(ndir)
        print(files)
