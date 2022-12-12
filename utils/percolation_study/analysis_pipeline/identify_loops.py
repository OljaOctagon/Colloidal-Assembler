import numpy as np
import gel_tools as gt
import networkx as nx
from itertools import islice
from itertools import cycle


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

    def get_bond_type_domains(bond_type, particle_did):
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

        for i, domain_i in enumerate(domains):
            for particle_j in domain_i:
                particle_did[particle_j] = np.array([bond_type, i])

        return domains, particle_did

    DG = nx.Digraph(G)
    loops = list(nx.simple_cycles(DG))

    def get_cycles(cycle_size, cylce_type, loops):
        '''
        cycle_size: int, size of the cycle length to be probed 
        cycle_type: int, cycle type to be probed: 
        -1:np-cycles, 0: mixed cycles, 1: parallel cycles
        loops: list of lists, where particles in a loop are in same sublist 
        '''

        len_loops = [len(loop) for loop in loops]
        cluster = np.array(
            [loop for loop in loops if len(loop) == cycle_size]).flatten()

    particle_did = np.zeros((N_particles, 2))

    domain_p, particle_did = get_bond_type_domains(1, particle_did)
    domain_np, particle_did = get_bond_type_domains(-1, particle_did)

    # write

    # write out particle_did
