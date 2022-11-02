import numpy as np
import h5py
import glob
import configparser

# create hd5f file
f = h5py.File("mytestfile.hdf5", 'a')

# read all bin binary files of directory
subdirs = glob.glob("double*")

for dir in subdirs:
    grp = f.create_group(dir)
    grp.attrs['id'] = dir

    config = configparser.ConfigParser()

    grp.attrs['number_particles'] = configparser['System']['Number_of_Particles']
    grp.attrs['packing_fraction'] = configparser['System']['Packing_Fraction']
    grp.attrs['temperature'] = configparser['System']['Temperature']
    grp.attrs['number_cells'] = configparser['System']['Number_of_Cells']

    grp.attrs['calculation_frequency'] = configparser['Output']['Calculation_Frequency']
    grp.attrs['checkpoint_frequency'] = configparser['Output']['Checkpoint_Frequency']
    grp.attrs['sigma_translation'] = configparser['Monte_Carlo_Parameters']['Sigma_Translation']
    grp.attrs['sigma_rotation'] = configparser['Monte_Carlo_Parameters']['Sigma_Rotation']

    grp.attrs['ptype'] = configparser['Rhombus']['rhombus_type']
    grp.attrs['energy'] = configparser['Rhombus']['Energy_Level']
    grp.attrs['patch_delta'] = configparser['Rhombus']['patch_delta']
    grp.attrs['patch_size'] = configparser['Rhombus']['patch_size']

    grp.attrs['ensemble'] = 'NVT'

    box = np.fromfile("Box_{}.bin".format(val))
    grp.attrs['box_xcenter'] = box[0]
    grp.attrs['box_ycenter'] = box[1]
    grp.attrs['box_lx'] = box[3]
    grp.attrs['box_ly'] = box[4]

    checkpoints = glob.glob("{}/positions.*bin".format(dir))
    check_point_values = np.sort(
        [int(point.split("_")[-1].split(".")[0]) for point in checkpoints]
    )

    N = int(grp.atrrs['number_particles'])
    for j, val in enumerate(check_point_values):

        x_extent = N
        pos_dim = 2
        orient_dim = 1
        y_extent = pos_dim + orient_dim
        dset = grp.create_dataset("checkpoint_{}".format(
            j), (x_extent, y_extent), dtype='f')

        patch_i = np.fromfile("patch_energy_{}.bi".format(val))
        dset.attrs['patch_type_0'] = patch_i[0]
        dset.attrs['patch_type_1'] = patch_i[1]
        dset.attrs['patch_type_2'] = patch_i[2]
        dset.attrs['patch_type_3'] = patch_i[3]

        dset.attrs['random_state-0'] = 'xxx'
        dset.attrs['random_state-1'] = 'xxx'

        pos_i = np.fromfile("positions_{}.bin".format(val))
        pos_i = np.reshape(pos_i, (-1, 3))
        pos_i = pos_i[:, :2]
        orient_i = np.fromfile("orientations_{}.bin".format(val))
        orient_i = np.reshape(orient_i, (-1, 5))[:, 4]

        orient_i = orient_i.reshape(-1, 1)
        config_i = np.concatenate((pos_i, orient_i), axis=1)
        dset[...] = config_i

f.close()
