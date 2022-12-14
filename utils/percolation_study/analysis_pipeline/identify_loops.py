import numpy as np
import gel_tools as gt
import networkx as nx
from itertools import islice
from itertools import cycle
import pandas as pd 
import configparser
import multiprocessing
import glob
import re 
import argparse
import h5py
from collections import defaultdict 
import gel_tools as gt 

def generator_from_fsys(dir_string):

    fsys_iterator = glob.glob(dir_string)
    for dir_i in fsys_iterator:
        config = configparser.ConfigParser()
        print(dir_i)
        config.read('{}para.ini'.format(dir_i))

        N = int(config['System']['Number_of_Particles'])
        phi = float(config['System']['Packing_Fraction'])
        temperature = float(config['System']['Temperature'])
        ptype = config['Rhombus']['rhombus_type']
        delta = config['Rhombus']['patch_delta']
        patch_size = config['Rhombus']['patch_size']

        pos_files = glob.glob('{}positions_*.bin'.format(dir_i))

        # get the last value from the string
        def g(x): return int(re.findall(r'\d+', x)[-1])

        mc_times = list(map(g, pos_files))

        last_time = np.max(mc_times)

        pos_file = "{}positions_{}.bin".format(dir_i, last_time)
        pos = np.fromfile(pos_file)
        pos = np.reshape(pos, (-1, 3))
        pos = pos[:, :2]
        orient_file = "{}orientations_{}.bin".format(dir_i, last_time)
        orient = np.fromfile(orient_file)
        orient = np.reshape(orient, (-1, 5))[:, 4]

        box_file = "{}Box_{}.bin".format(dir_i, last_time)
        box = np.fromfile(box_file)
        fid = dir_i
        
        file_patch_network = "{}patch_network.dat".format(dir_i)
        connections = gt.read_bonds(file_patch_network)[-1]

        yield (fid, ptype, phi, temperature, delta, last_time, pos, orient, box, connections)


def get_bond_type_domains(G, bond_type):
    '''input: 
        bond_type: int, 1: parallel, -1: non parallel 
        particle_did: np.array(), shape: (N_particles,2), 
        row entries: [int bond_type, int domain_id]

        output:
        domains: list of lists, sublist contains particle ids of domain i
        particle_did: updated particle_did, see input 
    '''

    S = nx.Graph(((source, target, attr) for source, target,
                attr in G.edges_iter(data=True) if attr['bond_type'] == bond_type))

    domains = list(nx.connected_components(S))

    return domains 

def get_domains(G, N_particles):

    domains_p = get_bond_type_domains(G,1)
    domains_np = get_bond_type_domains(G, -1)
    domains = domains_p + domains_np 
    l_p = len(domains_p)
    l_np = len(domains_np)
    
    domain_sizes = [ len(domain) for domain in domains]
    domain_types = np.concatenate((np.ones(l_p), -1*np.ones(l_np)))

    domain_sizes = np.reshape(np.array(domain_sizes),(-1,1))
    domain_types = np.reshape(np.array(domain_types),(-1,1))
    domain_sizes_types = np.concatenate((domain_sizes, domain_types), axis=1)

    particle_did = np.empty((N_particles, 2))
    for di, domain in enumerate(domains):
        for pi in domain:
            particle_did[pi] = domain_sizes_types[di]

    return domain_sizes_types, particle_did   


def get_loops(G,N_particles):
    DG = nx.Digraph(G)
    DG.remove_edges_from(nx.selfloop_edges(DG))

    loops = list(nx.simple_cycles(DG))
    loop_types = np.empty(len(loops))
    loop_sizes = np.empty(len(loops))

    particle_loop_types = np.empty((N_particles,2))

    for li, loop in enumerate(loops):
        l=len(loop)
        node_cycles = [ ( loop[i], loop[(i+1)%l] ) for i in range(l)]
        bt = [ G[n1][n2]["bond_type"] for (n1,n2) in node_cycles ]
        loop_sizes[li] = l 
    
        sum_bt=sum(bt)
        # parallel 
        if sum_bt == l:
            loop_types[li] = 1
            particle_loop_types[loop] = [1,l]
        if sum_bt == -l: 
            loop_types[li] = -1  
            particle_loop_types[loop] = [-1,l]

        if np.abs(sum_bt) < l:
            loop_types[li] = 0
            particle_loop_types[loop] = [0,l]

    loop_sizes = np.reshape(loop_sizes, (-1,1))
    loop_types = np.reshape(loop_types, (-1,1))
    loop_sizes_types = np.concatenate((loop_sizes,loop_types),axis=1)

    return loop_sizes_types, particle_loop_types


def identify_loops(connections, N_particles):
    '''
    connections: list of np.arrays(), l
    len(connections): total number of particle bonds
    per entry: len(connections[ientry] = 4
    connections[ientry][0]: connecting particle id1 
    connections[ientry][1]: connecting particle id2 
    connections[ientry][2]: connecting patch of id2 
    connections[ientry][3]: connecting patch of id2
    '''

    G = nx.Graph()
    G.add_edges_from(connections[:, :2])
    attrs = {}
    '''
    1: parallel
    -1: non parallel 
    '''

    orient_dict = {(0, 0): 1, (0, 1): -1, (0, 2): -1, (0, 3): 1,
                   (1, 0): -1, (1, 1): 1, (1, 2): 1, (1, 3): -1,
                   (2, 0): -1, (2, 1): 1, (2, 2): 1, (2, 3): -1,
                   (3, 0): 1, (3, 1): -1, (3, 2): -1, (3, 3): 1}

    bond_type_str = {0: "P", 1: "NP"}

    for entry in connections:
        attrs[(entry[0], entry[1])] = orient_dict[(entry[2], entry[3])]

    nx.set_edge_attributes(G, attrs, name="bond_type")
   
    domain_sizes_types, pdt = get_domains(G, N_particles)
    loop_sizes_types, plt = get_loops(G, N_particles)

    return domain_sizes_types, loop_sizes_types, pdt, plt


def calculate(vals):
    fid, ptype, phi, temperature, delta, last_time, pos, orient, box, connections = vals

    N_particles = len(pos)
    domain_sizes_types, loop_sizes_types, pdt, plt = identify_loops(connections, N_particles)

    meta = {}
    fid = fid.split("/")[1]
    print("get for fid: {}".format(fid))

    meta["fid"] = "{}_{}".format(fid, last_time)
    meta["ptype"] = ptype
    meta["phi"] = phi
    meta["temperature"] = temperature
    meta["delta"] = delta
    meta["last_time"] = last_time
    meta["run_id"] = 0
    meta["N_particles"] = N_particles
    
    results = {}
    results['positions'] = pos 
    results['orientations'] = orient 
    results['box'] = box 
    results['domain_sizes_types'] = domain_sizes_types
    results['loop_sizes_types'] = loop_sizes_types
    results['particle_domain_id'] = pdt
    results['particle_loop_id'] = plt 

    return meta,results

if __name__ == '__main__':

    # read data either through files system via glob or via db
    parser = argparse.ArgumentParser()
    parser.add_argument('-run_id', type=str)
    parser.add_argument('-ptype', type=str)
    parser.add_argument('-delta', type=float)
    parser.add_argument('-temperature', type=str)
    parser.add_argument('-phi', type=float)
    parser.add_argument('-ncores', type=int)

    args = parser.parse_args()
    para = defaultdict(lambda: '*')
    for arg in vars(args):
        if arg is not None:
            para[arg] = getattr(args,arg) 

    dir_string = "{}/{}_phi_{}_delta_{}_temp_{}/".format(
        para['ptype'], para['ptype'], para['phi'],
        para['delta'], para['temperature'],)

    print("files", dir_string)
    
    gen_fsys = generator_from_fsys(dir_string)

    N_CORES = int(args.ncores)
    N_CORES_MAX = 8

    if N_CORES > 1 and N_CORES <= N_CORES_MAX:
        print("Multiprocessing with {} cores".format(N_CORES))
        pool = multiprocessing.Pool(N_CORES)
        results = pool.map(calculate, gen_fsys)
        pool.close()
        pool.join()

    if N_CORES == 1:
        results=[]
        print("single core job")
        for vals in gen_fsys:
            results.append(calculate(vals))

    if N_CORES > N_CORES_MAX:
        print("Too many cores allocated, please do not use more than {} cores".format(
            N_CORES_MAX))
        exit()

    # create hd5f file
    f = h5py.File("domains_and_loops_{}.h5".format(args.run_id), 'w')
    for i, res in enumerate(results):
        grp = f.create_group(res[0]['fid'])

        for key in res[0]:
            grp.attrs[key] = res[0][key]

        for key in res[1]:
            dset = grp.create_dataset(key, res[1][key].shape, dtype='f')
            dset[...] = res[1][key]

    f.close()

