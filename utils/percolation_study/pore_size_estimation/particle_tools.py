import numpy as np 


def rotation_matrix(theta):
    rot_mat = np.zeros((2,2))

    rot_mat[0,0] = np.cos(theta) 
    rot_mat[0,1] = -np.sin(theta)
    rot_mat[1,0] = np.sin(theta)
    rot_mat[1,1] = np.cos(theta)

    return rot_mat


def get_orient(v, rot_mat):
    return rot_mat.dot(v)

def read_bonds(filen):
    first_line_pair = [0,0,0,0]
    cut=False
    with open(filen, 'r') as f:
        network_list = []
        for line in f:
            if "#" in line:
                network_list.append([])
                first_line_pair = [0,0,0,0]
                cut=False

            else:
                line_counter=len(network_list[-1])
                pairs = list(map(int, line.split(" ")))
                if pairs == first_line_pair or cut==True:
                    cut=True
                else:
                    network_list[-1].append(np.array(pairs))

                if line_counter == 0:
                    first_line_pair = pairs
    network_list = [ np.array(item) for item in network_list]

    return network_list


def read_config(fdir, val):
    pos_i = np.fromfile("{}/positions_{}.bin".format(fdir, val))
    pos_i = np.reshape(pos_i, (-1,3))
    pos_i = pos_i[:,:2]

    orient_i = np.fromfile("{}/orientations_{}.bin".format(fdir, val))
    orient_i = np.reshape(orient_i, (-1,5))[:,4]


    return pos_i, orient_i