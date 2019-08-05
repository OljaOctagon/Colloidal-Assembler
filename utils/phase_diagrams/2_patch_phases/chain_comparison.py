import numpy as np
import pandas as pd
import argparse
import glob
import os
import re
import networkx as nx
from scipy import linalg as la
import matplotlib.pyplot as plt
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

S_threshold = 0.4
Alignment_threshold = 0.85
pos_threshold = np.power(1.5,2) + (2+2*np.cos(np.pi/3.))


def calculate_chain_neighbours(i, chain, all_chains, pos, boxl):
    chain_neighbours = []
    for j, chain_j in enumerate(all_chains):
        for elem_i in chain:
            for elem_j in chain_j:
                dist_v = (pos[elem_i] - pos[elem_j])[:2]
                dist_v = dist_v - boxl*np.rint(dist_v/boxl)
                dist = np.power(dist_v[0],2) + np.power(dist_v[1],2)

                if dist < pos_threshold and i!=j:
                    chain_neighbours.append(j)
                    break

            else:
                continue
            break

    return chain_neighbours

def calculate_local_nematic_op(i, chain_neighbours, dir_u):

    if len(chain_neighbours) > 0:
        denom = 1/len(chain_neighbours)

    else:
        denom = 1 
    f = lambda i,j: 2*np.power(np.dot(dir_u[i], dir_u[j]),2) -  1
    sum_terms = [ f(i,neigh) for neigh in chain_neighbours]
    local_nematic = denom*np.sum(sum_terms)

    return local_nematic


def get_alignment(i, j, dir_u):
    is_aligned = False
    if np.fabs(np.dot(dir_u[i], dir_u[j])) > Alignment_threshold:
        is_aligned = True

    return is_aligned


def get_length_type(l_chain):

    if l_chain < 5:
        length_type = 10

    if l_chain >=5 and l_chain<10:
        length_type = 11

    if l_chain >=10 and l_chain<15:
        length_type = 12

    if l_chain >=15 and l_chain<20:
        length_type = 13

    if l_chain >=20:
        length_type = 14

    if l_chain >=30:
        length_type = 15

    return length_type


# calulate director of a chain. We take the director as the normed distance vector
# from the first to last particle 
def calculate_director(Gc_sorted,pos, boxl):

    vec = np.zeros((len(Gc_sorted),2))

    for i in range(1,len(Gc_sorted)):
        distance = (pos[Gc_sorted[i]] - pos[Gc_sorted[i-1]])[:2]
        distance = distance - boxl*np.rint(distance/boxl)

        vec[i] = vec[i-1] + distance

    full_vec = vec[-1]
    pos_zero = pos[Gc_sorted[0]]

    n_vec = np.linalg.norm(vec[-1])
    dir_vec = vec[-1]/n_vec
    return dir_vec, full_vec[:2], pos_zero[:2]



def calculate_bond_vector_angles(Gc_sorted, pos, boxl):

    L = len(Gc_sorted)
    bond_vec = np.zeros((L-1,2))
    bond_angle = np.zeros((L-2,1))

    for i in range(1,L):
        bond_vec[i-1] = (pos[Gc_sorted[i]] - pos[Gc_sorted[i-1]])[:2]
        bond_vec[i-1] = bond_vec[i-1] - boxl*np.rint(bond_vec[i-1]/boxl)

    N_bonds = L-1
    for j in range(1,N_bonds):

        ab_scalar = np.dot(bond_vec[j-1], bond_vec[j])
        norm_a = np.linalg.norm(bond_vec[j-1])
        norm_b = np.linalg.norm(bond_vec[j])

        bond_angle[j-1] = np.arccos(ab_scalar/(norm_a*norm_b))

    return bond_vec, bond_angle


def calculate_nematic_order(dir_u):
    # calculate nematic tensor 
    Nm = len(dir_u)
    Q_ab = np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            Q_ab[i,j] = (1/Nm) * np.sum(2*dir_u[:,i]*dir_u[:,j] - np.eye(2)[i,j])


    eigen_val, eigen_vec = la.eig(Q_ab)
    print('eigen_val', eigen_val)
    nematic_op = max(eigen_val).real
    return nematic_op 

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
                'bend','bend_parallel', 'bend_non_parallel',
                'end_to_end_distance'
                'sequence','bond_vec', 'bond_angle']

    df = pd.DataFrame(columns=features)

    features_1 = ['mu', 'energy', 'topology', 'delta',
                  'nematic_op','local_nematic_op',
                  'fraction_largest_nematic']

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
        full_vec = []
        pos_zero = []

        id_length_type = np.ones(1000)*10


        all_chains = []
        all_bond_angles = []
        

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

                all_chains.append(Gc_sorted)

                # calculate director
                dir_u_i, full_vec_i, pos_zero_i = calculate_director(Gc_sorted,pos, boxl)
                dir_u.append(dir_u_i)
                full_vec.append(full_vec_i)
                pos_zero.append(pos_zero_i)

                # asign length category to particle ids
                length_type = get_length_type(len(chain))
                id_length_type[Gc_sorted] = length_type


                cluster_size = len(Gc_sorted)
                # bend, bend_between_kink, straight_distribution
                bend, bend_p, bend_np, sequence = calculate_flexibility(
                    Gc_sorted, subarr,orient)

                end_to_end_distance = np.sqrt(np.power(full_vec_i[0],2) + np.power(full_vec_i[1],2))

                bond_vec, bond_angle = calculate_bond_vector_angles(Gc_sorted, pos, boxl)

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
                                                'end_to_end_distance': [end_to_end_distance],
                                                 'sequence': [sequence],
                                                 'bond_vec':[bond_vec],
                                                 'bond_angle':[bond_angle]})])

        dir_u = np.array(dir_u)
        nematic_op = calculate_nematic_order(dir_u)
       
        is_nematic = []
        NN_list = []

        is_nematic = np.zeros(len(all_chains), dtype=bool)

        for i, chain in enumerate(all_chains):
   
            # get all neigbours of chains 
            chain_neighbours = calculate_chain_neighbours(i, chain, all_chains, pos, boxl)
            NN_list.append(chain_neighbours)
            # local order parameters
            S_cid = calculate_local_nematic_op(i, chain_neighbours,dir_u)

            if S_cid > S_threshold:
                is_nematic[i] = True


        N_graph_list = []
        for i, chain in enumerate(all_chains):
            for j in NN_list[i]:
                if is_nematic[i] and is_nematic[j]:
                    is_aligned = get_alignment(i, j, dir_u)
                    if is_aligned:
                        N_graph_list.append([i,j])

        N_largest = 1
        local_nematic_op = 0
        if(len(N_graph_list)>0):
            H = nx.Graph()
            H.add_edges_from(N_graph_list)
            domains = nx.connected_components(H)
            clusters= [ list(domain) for domain in list(domains)]
            domain_length = np.array([ len(domain) for domain in clusters ])
            d_id = np.argmax(domain_length)
            # size of largest cluster 
            N_largest = len(clusters[d_id])
            local_nematic_op = calculate_nematic_order(dir_u[clusters[d_id]])

        fraction_largest = N_largest/len(all_chains)
        print(mu, energy, topology, delta, fraction_largest, local_nematic_op, nematic_op)


        dg = pd.concat([dg,pd.DataFrame({'mu': [mu],
                                        'energy': [energy],
                                        'topology': [topology],
                                         'delta': [delta],
                                         'nematic_op': [nematic_op],
                                         'local_nematic_op': [local_nematic_op],
                                         'fraction_largest_nematic': [fraction_largest]})])

        with open(dir+'/bend_op.dat', 'w') as f:
            f.write(str(1000)+'\n')
            f.write("Particles of frame\n")
            for i in range(1000):
                f.write(str(bend_per_particle[i]+10)+'\n')

        with open(dir+'/length_op.dat', 'w') as f:
            f.write(str(1000)+'\n')
            f.write("Particles of frame\n")
            for i in range(1000):
                f.write(str(id_length_type[i])+'\n')

        #for chain in all_chains:
        #    plt.scatter(pos[chain][:,0], pos[chain][:,1])
        #    plt.plot(pos[chain][:,0], pos[chain][:,1])

        #plt.savefig(dir+"/director_particle_scatter.png")
        #plt.close()

    with open(args.o, 'w') as f:
        df.to_csv(f)

    with open("nematic_op.csv", 'w') as f:
        dg.to_csv(f)

