import numpy as np
import itertools
def read_bonds(filen):
    first_line_pair = (0,0,0,0)
    cut=False
    with open(filen, 'r') as f:
        network_list = []
        for line in f:
            if "#" in line:
                network_list.append([])
                first_line_pair = (0,0,0,0)
                cut=False
            else:
                line_counter=len(network_list[-1])
                pairs = list(map(int, line.split(" ")))
                if pairs == first_line_pair or cut==True:
                    cut=True
                else:
                    network_list[-1].append((pairs[0], pairs[1], pairs[2], pairs[3]))

                if line_counter == 0:
                    first_line_pair = pairs
    #network_list = [ np.array(item) for item in network_list]
    return network_list

def runs_of_ones_list(seq):
    return [sum(g) for b, g in itertools.groupby(seq) if b==1]

def calculate_cluster_stability(network_list):
    # add padding so we can cast network_list to array
    lengths = [ len(item) for item in network_list]
    max_length = max(lengths)
    padding_line = (-1,-1,-1,-1)
    for item_list in network_list:
        padding = [padding_line]*(max_length-len(item_list))
        item_list.extend(padding)

    network_arr = np.array(network_list, dtype = 'i,i,i,i')
    dimx,dimy = network_arr.shape
    uniques,inverse = np.unique(network_arr,return_inverse=True)
    inverse = np.reshape(inverse, (dimx,dimy))
    print(inverse)
    lengths = []
    for i in range(len(uniques)):
        sequence = np.where(inverse == i)[0]
        diff_seq = np.diff(sequence)
        length_runs = runs_of_ones_list(diff_seq)
        if length_runs:
            lengths.append(length_runs)

    return uniques, lengths

def calc_stats_per_type(m,n, uniques, lengths):

    flatten = lambda l: [item for sublist in l for item in sublist]
    f_type = lambda i,m,n: ((uniques[i][2],uniques[i][3]) == (m,n)) or ((uniques[i][2],uniques[i][3]) == (n,m))

    #print(lengths)
    c_type_lengths = [ lengths[i] for i in range(len(lengths)) if f_type(i,m,n) ]
    c_type_lengths = np.array(flatten(c_type_lengths))

    c_mean = 0
    c_std = 0

    if c_type_lengths.size:
        c_mean = np.mean(c_type_lengths)
        c_std = np.std(c_type_lengths)

    return c_mean, c_std


if __name__ == '__main__':
    pn_file = "patch_network.dat"
    network_list = read_bonds("patch_network.dat")

    # cluster stability
    uniques, lengths = calculate_cluster_stability(network_list)
    p_mean, p_std = calc_stats_per_type(1,1, uniques, lengths)
    s_mean, s_std = calc_stats_per_type(0,0, uniques, lengths)
    l_mean, l_std = calc_stats_per_type(0,1, uniques, lengths)

    print("p-type", p_mean, p_std)
    print("s-type", s_mean, s_std)
    print("l-type", l_mean, l_std)

    # cluster distribution of size 2 clusters
    flatten = lambda l: [item for sublist in l for item in sublist]

    network_list = read_bonds("patch_network.dat")
    flat_network_list = flatten(network_list)

    N_total = len(flat_network_list)
    arr_flat = np.array(flat_network_list)
    
    N_p = len( arr_flat[np.where((arr_flat[:,2] == 1) & (arr_flat[:,3] == 1))[0]] )
    print("N_p", N_p/N_total)

    N_s = len( arr_flat[np.where((arr_flat[:,2] == 0) & (arr_flat[:,3] == 0))[0]] )
    print("N_s", N_s/N_total)

    N_l1 = len( arr_flat[np.where((arr_flat[:,2] == 0) & (arr_flat[:,3] == 1))[0]] )
    N_l2 = len( arr_flat[np.where((arr_flat[:,2] == 1) & (arr_flat[:,3] == 0))[0]] )
    print("N_l", (N_l1+N_l2)/N_total)
