import numpy as np
import h5py
import glob
import configparser

# create hd5f file
f = h5py.File("test.h5", 'w')

# read all bin binary files of directory
subdirs = glob.glob("double*")

for dir in subdirs:
    print(dir)
    grp = f.create_group(dir)
    grp.attrs['id'] = dir

    config = configparser.ConfigParser()
    config.read("{}/para.ini".format(dir))

    grp.attrs['number_particles'] = config['System']['Number_of_Particles']
    grp.attrs['packing_fraction'] = config['System']['Packing_Fraction']
    grp.attrs['temperature'] = config['System']['Temperature']
    grp.attrs['number_cells'] = config['System']['Number_of_Cells']

    grp.attrs['sigma_translation'] = config['Monte_Carlo_parameters']['Sigma_Translation']
    grp.attrs['sigma_rotation'] = config['Monte_Carlo_parameters']['Sigma_Rotation']

    grp.attrs['ptype'] = config['Rhombus']['rhombus_type']
    grp.attrs['energy'] = config['Rhombus']['Energy_Level']
    grp.attrs['patch_delta'] = config['Rhombus']['patch_delta']
    grp.attrs['patch_size'] = config['Rhombus']['patch_size']

    grp.attrs['ensemble'] = 'NVT'

    checkpoints = glob.glob("{}/positions_*bin".format(dir))
    check_point_values = np.sort(
        [int(point.split("_")[-1].split(".")[0]) for point in checkpoints]
    )
    val = check_point_values[0]
    box = np.fromfile("{}/Box_{}.bin".format(dir, val))
    grp.attrs['box_xcenter'] = box[0]
    grp.attrs['box_ycenter'] = box[1]
    grp.attrs['box_lx'] = box[3]
    grp.attrs['box_ly'] = box[4]

    N = int(grp.attrs['number_particles'])
    for j, val in enumerate(check_point_values):

        x_extent = N
        pos_dim = 2
        orient_dim = 1
        y_extent = pos_dim + orient_dim
        sgrp = f.create_group("{}/check_point_{}".format(dir, j))

        patch_i = np.fromfile(
            "{}/patch_energy_{}.bin".format(dir, val), dtype=np.int32)

        sgrp.attrs['MC time'] = val
        sgrp.attrs['patch_type_0'] = patch_i[0]
        sgrp.attrs['patch_type_1'] = patch_i[1]
        sgrp.attrs['patch_type_2'] = patch_i[2]
        sgrp.attrs['patch_type_3'] = patch_i[3]
        sgrp.attrs['random_state-0'] = 'xxx'
        sgrp.attrs['random_state-1'] = 'xxx'

        dset_pos = sgrp.create_dataset(
            "positions_{}".format(j), (x_extent, pos_dim), dtype='f')

        dset_orient = sgrp.create_dataset(
            "orientations_{}".format(j), (x_extent), dtype='f')

        pos_i = np.fromfile("{}/positions_{}.bin".format(dir, val))
        pos_i = np.reshape(pos_i, (-1, 3))
        pos_i = pos_i[:, :2]
        orient_i = np.fromfile("{}/orientations_{}.bin".format(dir, val))
        orient_i = np.reshape(orient_i, (-1, 5))[:, 4]

        dset_pos[...] = pos_i
        dset_orient[...] = orient_i

f.close()
