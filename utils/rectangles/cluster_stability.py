import numpy as np

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

if __name__ == '__main__':
    pn_file = "patch_network.dat"
    network_list = read_bonds("patch_network.dat")

    # add padding so we can cast network_list to array
    lengths = [ len(item) for item in network_list]
    max_length = max(lengths)
    padding_line = (-1,-1,-1,-1)
    for item_list in network_list:
        padding = [padding_line]*(max_length-len(item_list))
        item_list.extend(padding)

    network_arr = np.array(network_list, dtype = 'i,i,i,i')
    print(network_arr.shape)
    uniques,indices,inverse, counts  = np.unique(network_arr,return_index=True, return_inverse=True, return_counts=True)
    print(inverse.shape)
    #print(indices[0])
