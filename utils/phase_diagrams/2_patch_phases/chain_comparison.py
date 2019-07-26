import numpy as np
import pandas as pd
import argparse
import glob
import os
import re
import networkx as nx

cluster_type_dict = {
    '6-star': 6,
    '5-star': 5,
    'box': 3,
    'open-box':3,
    'open-6-star':6,
    'open-5-star': 5}

patch_dict = {
    'mouse': [0,2],
    'manta': [0,1]
}


# calulate director of a chain. We take the director as the normed distance vector
# from the first to last particle 
def calculate_director(Gc_sorted,pos):
    vec = pos[Gc_sorted[-1]] - pos[Gc_sorted[0]]
    n_vec = np.linalg.norm(vec)
    vec = vec/n_vec
    return vec[:2]


def parse_config(directory):
    """returns positions, box lengths of max_checkpoint as arrays. """

    filelist = glob.glob(directory+"/positions_*.bin")
    g = lambda x: x.split("_")[-1].split(".")[0]

    checkpoints = [ int(g(l)) for l in filelist]
    max_cp = max(checkpoints)

    pos = np.fromfile("{}/{}_{}.bin".format(directory, 'positions', max_cp))
    pos = np.reshape(pos, (-1,3))

    orient = np.fromfile("{}/{}_{}.bin".format(directory,
                                                 'orientations', max_cp))
    orient = np.reshape(orient, (-1,5))[:,4]

    lbox = np.fromfile("{}/{}_{}.bin".format(directory, 'Box', max_cp))
    box_l = np.array([lbox[3], lbox[4]])

    return pos,orient,box_l



def order_chain(Gc, subarr):
    uniques, counts = np.unique(subarr[:,:2], return_counts=True)
    chain_ends  = np.where(counts == 1)[0]
    cluster_type='chain'

    if len(chain_ends) == 2:
        Gc_sorted = list(nx.all_simple_paths(Gc,uniques[chain_ends[0]], uniques[chain_ends[1]]))[0]
        cluster_type='chain'

    elif len(chain_ends) == 0:
        paths = list(nx.all_simple_paths(Gc,subarr[0,0], subarr[0,1]))
        for path in paths:
            if len(path) > 2:
                Gc_sorted = path

        cluster_type='loop'

    else:
        raise ValueError("branched graph! Something is wrong")

    return Gc_sorted, cluster_type

def calculate_flexibility(Gc_sorted,arr,orient):
    parallel_angles=np.array([0, np.pi, 2*np.pi])
    non_parallel_angles = np.array([np.pi/3, 2*np.pi/3, 4*np.pi/3, 5*np.pi/3])
    delta_angle = np.fabs(orient[Gc_sorted[:-1]] - orient[Gc_sorted[1:]])

    seq=[]
    delta_angle_p = []
    delta_angle_np = []

    for i,angle in enumerate(delta_angle):
        pid = Gc_sorted[i]
        nid = Gc_sorted[i+1]
        connections = arr[
            np.where(((arr[:,0]==pid)&(arr[:,1]==nid))|((arr[:,0]==nid)&(arr[:,1]==pid)))][0][2:]

        # parallel bond
        if connections[0] == connections[1]:
            rel_angle = np.fabs(parallel_angles - angle)
            delta_angle[i]=np.min(rel_angle)
            delta_angle_p.append(delta_angle[i])
            seq.append('p')
        # NON paralle bond, kink.
        else:
            rel_angle = np.fabs(non_parallel_angles - angle)
            delta_angle[i]=np.min(rel_angle)
            delta_angle_np.append(delta_angle[i])
            seq.append('np')

    flex =np.average(delta_angle)
    flex_p = np.average(delta_angle_p)
    #flex_np = 0
    flex_np = np.average(delta_angle_np)
    return flex, flex_p, flex_np, seq


def get_kink_density(arr):
    nkinks = len(
        [ item for item in arr if item[2] != item[3]])
    kink_density  = nkinks / len(arr)

    return kink_density

if __name__ == '__main__':

    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-ptype', type=str, default='mouse', choices=['manta','mouse', 'parallel'])
    parser.add_argument('-tsize', type=int, default=3)
    parser.add_argument('-o', type=str, default='m_chains.csv')
    args  = parser.parse_args()

    filen = "patch_network.dat"
    filenpt = 'NPT_OUT.txt'
    dirlist = glob.glob("mu_*")
    pwd = os.getcwd()


    features = ['mu', 'energy', 'topology',
                'delta', 'percolation_coeff',
                'packing_fraction', 'kink_density',
                'bend','bend_parallel', 'bend_non_parallel', 'sequence']

    df = pd.DataFrame(columns=features)


    features_1 = ['mu', 'energy', 'topology', 'delta', 'nematic_op']
    dg = pd.DataFrame(columns=features_1)


    for dir in dirlist:
        numbers = re.findall(r"[-+]?\d*\.\d+|\d+", dir)
        mu = numbers[0]
        energy = numbers[1]
        delta = numbers[2]
        topology = 'symm'
        p = 'Asymm'
        if re.search(p, dir):
            topology = 'Asymm'

        arr = pd.read_csv('{}/{}'.format(dir, filen),
                            header=None,
                            delim_whitespace=True).values

        arr_npt = pd.read_csv('{}/{}'.format(dir, filenpt),
                                header=None,
                                delim_whitespace=True).values

        N_particles = arr_npt[-1,2]
        packing_fraction = arr_npt[-1,4]

        # parse config
        pos, orient, boxl  = parse_config(dir)

        # cluster size distribution, find biggest cluster and percolation

        G = nx.Graph()
        G.add_edges_from(arr[:,:2])
        domains = nx.connected_components(G)
        chains=list(domains)
        domain_length = np.array([ len(domain) for domain in chains ])
        d_id = np.argmax(domain_length)
        # list of particles in longest chain 
        particles_longest_chain = np.array(chains[d_id])

        # size of largest cluster 
        N_largest = len(chains[d_id])
        percolation_coeff = N_largest/N_particles

        #TODO treat cycles as cycles and not as chains!

        bend_per_particle = np.zeros(1000)
        # from previous calcs
        max_bend = 29.782761453318788


        dir_u = []

        for domain in nx.connected_component_subgraphs(G):
            chain = list(domain)
            if(len(chain)>args.tsize):
                # kink density
                subarr=[]
                for pid in chain:
                    sarr = arr[np.where( (arr[:,0] == pid) | (arr[:,1]== pid))]
                    for item in sarr:
                        if item.size > 0:
                                subarr.append(item)



                subarr = np.unique(subarr, axis=0)
                kink_density = get_kink_density(subarr)

                # sort chain graph: Gc_sorted: sorted list of chain_i
                Gc_sorted, cluster_type = order_chain(domain,subarr)

                # calculate director
                dir_u.append(calculate_director(Gc_sorted))

                cluster_size = len(Gc_sorted)
                # bend, bend_between_kink, straight_distribution
                bend, bend_p, bend_np, sequence = calculate_flexibility(
                    Gc_sorted, subarr,orient)

                bend_per_particle[chain] = np.rint(((bend*360/(2*np.pi))/max_bend)*5)
                df = pd.concat([df,pd.DataFrame({'mu': [mu],
                                                'energy': [energy],
                                                'topology': [topology],
                                                'delta': [delta],
                                                'percolation_coeff': [percolation_coeff],
                                                'packing_fraction': [packing_fraction],
                                                 'cluster_size': [cluster_size],
                                                 'cluster_type': [cluster_type],
                                                 'kink_density': [kink_density],
                                                'bend': [bend],
                                                'bend_parallel': [bend_p],
                                                'bend_non_parallel': [bend_np],
                                            'sequence': [sequence]})])


        nematic_op = calculate_nematic_order(np.array(dir_u))

        dg = pd.concat([dg,pd.DataFrame({'mu': [mu],
                                        'energy': [energy],
                                        'topology': [topology],
                                         'delta': [delta],
                                         'nematic_op': [nematic_op]})])

        with open(dir+'/bend_op.dat', 'w') as f:
            f.write(str(1000)+'\n')
            f.write("Particles of frame\n")
            for i in range(1000):
                f.write(str(bend_per_particle[i]+10)+'\n')

    with open(args.o, 'w') as f:
        df.to_csv(f)

    with open("nematic_op.csv", 'w') as f:
        dg.to_csv(f)

