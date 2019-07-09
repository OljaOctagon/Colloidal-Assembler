import numpy as np
import networkx as nx
import pandas as pd


def filter_bonds(data,pid1,pid2):
    fdata1 = data[np.logical_or(data[:,2]==pid1, data[:,2]==pid2)]
    fdata1 = fdata[np.logical_or(fdata[:,3]==pid1, fdata[:,3]==pid2)]
    fdata = np.concatenate((fdata1,fdata2))
    return fdata


def filter_bonds_parallel_dma(data,pid1,pid2):
    # for dma-as2: that's parallel bonds 
    fdata1 = data[np.logical_and(data[:,2]==pid1, data[:,3]==pid1)]
    fdata2 = data[np.logical_and(data[:,2]==pid2, data[:,3]==pid2)]
    fdata = np.concatenate((fdata1,fdata2))
    return fdata

def filter_bonds_non_parallel_dma(data,pid1,pid2):
    # for dma-as2: that's parallel bonds 
    fdata1 = data[np.logical_and(data[:,2]==pid1, data[:,3]==pid2)]
    fdata2 = data[np.logical_and(data[:,2]==pid2, data[:,3]==pid1)]
    fdata = np.concatenate((fdata1,fdata2))
    return fdata

def filter_non_parallel_dmo(data):
    fdata1 = data[np.logical_and(data[:,2]==0, data[:,3]==2)]
    fdata2 = data[np.logical_and(data[:,2]==2, data[:,3]==0)]
    fdata3 = data[np.logical_and(data[:,2]==1, data[:,3]==3)]
    fdata4 = data[np.logical_and(data[:,2]==3, data[:,3]==1)]

    fdata5 = data[np.logical_and(data[:,2]==0, data[:,3]==1)]
    fdata6 = data[np.logical_and(data[:,2]==1, data[:,3]==0)]
    fdata7 = data[np.logical_and(data[:,2]==2, data[:,3]==3)]
    fdata8 = data[np.logical_and(data[:,2]==3, data[:,3]==2)]
    fdata = np.concatenate((fdata1,fdata2,fdata3, fdata4,fdata5,fdata6,fdata7,fdata8))
    return fdata

def filter_dma_asymm1(data, kind):
    if kind == 'np_blue':
        fdata1 = data[np.logical_and(data[:,2]==0, data[:,3]==1)]
        fdata2 = data[np.logical_and(data[:,2]==1, data[:,3]==0)]
        fdata = np.concatenate((fdata1,fdata2))

    if kind == 'np_red':
        fdata1 = data[np.logical_and(data[:,2]==3, data[:,3]==2)]
        fdata2 = data[np.logical_and(data[:,2]==2, data[:,3]==3)]
        fdata = np.concatenate((fdata1,fdata2))
    return fdata

def filter_dmo_symm1(data, kind):
    if kind == 'np_blue':
        fdata1 = data[np.logical_and(data[:,2]==0, data[:,3]==2)]
        fdata2 = data[np.logical_and(data[:,2]==2, data[:,3]==0)]
        fdata = np.concatenate((fdata1,fdata2))

    if kind == 'np_red':
        fdata1 = data[np.logical_and(data[:,2]==1, data[:,3]==3)]
        fdata2 = data[np.logical_and(data[:,2]==3, data[:,3]==1)]
        fdata = np.concatenate((fdata1,fdata2))
    return fdata

def filter_dmo_symm2(data, kind):

    if kind == 'np_red':
        fdata1 = data[np.logical_and(data[:,2]==1, data[:,3]==3)]
        fdata2 = data[np.logical_and(data[:,2]==3, data[:,3]==1)]
        fdata = np.concatenate((fdata1,fdata2))

    if kind == 'p_red':
        fdata1 = data[np.logical_and(data[:,2]==1, data[:,3]==1)]
        fdata2 = data[np.logical_and(data[:,2]==3, data[:,3]==3)]
        fdata = np.concatenate((fdata1,fdata2))

    if kind == 'np_blue':
        fdata1 = data[np.logical_and(data[:,2]==0, data[:,3]==2)]
        fdata2 = data[np.logical_and(data[:,2]==2, data[:,3]==0)]
        fdata = np.concatenate((fdata1,fdata2))

    if kind == 'p_blue':
        fdata1 = data[np.logical_and(data[:,2]==0, data[:,3]==0)]
        fdata2 = data[np.logical_and(data[:,2]==2, data[:,3]==2)]
        fdata = np.concatenate((fdata1,fdata2))
    return fdata

def calculate_domains(nid,T):
    G=nx.Graph()
    G.add_edges_from(nid)
    domains = list(nx.connected_components(G))
    particles=[]
    for domain in domains:
        domain = np.array(list(domain))
        if len(domain)>=T:
            particles.extend(domain)
    return particles



def make_ziqzaqs(nid):

    arr_orient = pd.read_csv("color_orient.dat", header=None, delim_whitespace=True).values

    def select_orientations(arr_orient,nid,kick_out_id):
        idx = np.argwhere(arr_orient!=kick_out_id)
        print(idx)
        print(nid.shape)
        sub_nid = np.array([ elem for elem in nid if elem[0] in idx and elem[1] in idx ])
        print(sub_nid)
        return(sub_nid)

    particles=[]
    for o_i in range(1,4):
        sub_nid = select_orientations(arr_orient,nid, o_i)
        part=calculate_domains(sub_nid,T)
        particles.extend(part)

    return particles

def make_cycles(nid,size):
    G = nx.Graph()
    G.add_edges_from(nid)
    #loops = nx.cycle_basis(G)
    DG = nx.DiGraph(G)
    loops = list(nx.simple_cycles(DG))

    len_loops = [ len(loop) for loop in loops]
    cluster  = np.array([ loop for loop in loops if len(loop)==size ]).flatten()

    return cluster

def dmo_asymm1_route(data):

    mixed_c6 = make_cycles(data[:,:2],6)
    mixed_c7 = make_cycles(data[:,:2],7)
    mixed_c8 = make_cycles(data[:,:2],8)
    mixed_c9 = dimake_cycles(data[:,:2],9)
    mixed_c10 = make_cycles(data[:,:2],10)

    color=np.zeros(1000)
    color[mixed_c6] = 3
    color[mixed_c7] = 3
    color[mixed_c8] = 3
    color[mixed_c9] = 3
    color[mixed_c10] = 3
    return color

def dma_asymm2_route(data):
    filtered_blue = filter_bonds_non_parallel_dma(data,0,1)
    filtered_red = filter_bonds_non_parallel_dma(data,2,3)

    blue_c3= make_cycles(filtered_blue[:,:2],3)
    blue_c6= make_cycles(filtered_blue[:,:2],6)
    red_c3 = make_cycles(filtered_red[:,:2],3)
    red_c6 = make_cycles(filtered_red[:,:2],6)

    filter_parallel_blue = filter_bonds_parallel_dma(data,0,1)
    filter_parallel_red = filter_bonds_parallel_dma(data,2,3)

    parallel_blue = calculate_domains(filter_parallel_blue[:,:2],2)
    parallel_red = calculate_domains(filter_parallel_red[:,:2],2)
    
    non_parallel_blue = calculate_domains(filtered_blue[:,:2],2)
    non_parallel_red = calculate_domains(filtered_red[:,:2],2)

    color = np.zeros(1000)

    color[parallel_blue] = 3 
    color[parallel_red] = 3 
    color[non_parallel_red] = 1
    color[non_parallel_blue] = 1
    #color[red_c6] = 2
    #color[blue_c6] = 2

    return color

def dmo_asymm2_route(data):

    filtered_blue = filter_bonds_non_parallel_dma(data,0,2)
    filtered_red = filter_bonds_non_parallel_dma(data,1,3)

    non_parallel_blue = calculate_domains(filtered_blue[:,:2],2)
    non_parallel_red = calculate_domains(filtered_red[:,:2],2)

    filter_parallel_blue = filter_bonds_parallel_dma(data,0,2)
    filter_parallel_red = filter_bonds_parallel_dma(data,1,3)

    parallel_blue = calculate_domains(filter_parallel_blue[:,:2],2)
    parallel_red = calculate_domains(filter_parallel_red[:,:2],2)
    
    color = np.zeros(1000)

    color[parallel_blue] = 3 
    color[parallel_red] = 3 
    
    color[non_parallel_blue] = 1 
    color[non_parallel_red] = 1

    return color

def dmo_symm1_route(data):
    # stars with satisfied and unsatisfied bonds
    filtered_blue = filter_bonds(data,0,2)
    filtered_red = filter_bonds(data,1,3)

    blue_c5= make_cycles(filtered_blue[:,:2],5)
    blue_c6= make_cycles(filtered_blue[:,:2],6)
    red_c5 = make_cycles(filtered_red[:,:2],5)
    red_c6 = make_cycles(filtered_red[:,:2],6)

    np_data = filter_non_parallel_dmo(data)
    mixed_c6 = make_cycles(np_data[:,:2],6)


    filter_dmo_b = filter_dmo_symm1(data, 'np_blue')
    ziqzaq1 = calculate_domains(filter_dmo_b[:,:2],2)
    filter_dmo_r = filter_dmo_symm1(data, 'np_red')
    ziqzaq2 = calculate_domains(filter_dmo_r[:,:2],2)
    

    color = np.zeros(1000)

    color[ziqzaq1] = 3 
    color[ziqzaq2] = 3 
    color[blue_c5] = 1
    color[mixed_c6] = 2
    color[blue_c6] = 2
    color[red_c5] = 1
    color[red_c6] = 2
    

    return color

def dma_asymm1_route(data):
    filtered_p = filter_dmo_symm2(data,'p_red')
    parallel = calculate_domains(filtered_p[:,:2],3)

    filtered_npb = filter_dma_asymm1(data,'np_blue')
    filtered_npr= filter_dma_asymm1(data,'np_red')
    c3_blue = make_cycles(filtered_npb[:,:2],3)
    c3_red = make_cycles(filtered_npr[:,:2],3)
    color = np.zeros(1000)
    color[parallel] =  2
    color[c3_blue] = 1
    color[c3_red] = 1 

    return color

def dmo_symm2_route(data):
    filtered_p = filter_dmo_symm2(data,'p_red')
    parallel = calculate_domains(filtered_p[:,:2],3)
    filtered_np = filter_dmo_symm2(data,'np_red')
    n_parallel = calculate_domains(filtered_np[:,:2],3)

    filtered_pr = filter_dmo_symm2(data,'np_red')
    c6= make_cycles(filtered_pr[:,:2],6)

    color = np.zeros(1000)
    color[parallel] =  2
    color[n_parallel] = 3
    color[c6] = 1

    return color

# patch_network data:
# particle id1 (rid1), particle id2 (rid2), patch id (pid1), patch id (pid2)
data = pd.read_csv("patch_network.dat", header=None, delim_whitespace=True).values

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-type')
args=parser.parse_args()

color = np.zeros(1000)

if args.type == 'dmo-symm1':
    color = dmo_symm1_route(data)

if args.type == 'dmo-asymm1':
    color = dmo_asymm1_route(data)

if args.type == 'dmo-symm2':
    color = dmo_symm2_route(data)

if args.type == 'dma-asymm1':
    color = dma_asymm1_route(data)

if args.type == 'dma-asymm2':
    color = dma_asymm2_route(data)

if args.type == 'dmo-asymm2':
    color = dmo_asymm2_route(data)

with open("color_op.dat", 'w') as fhandle:
    fhandle.write("1000\n")
    fhandle.write("Particles of frame X\n")
    for col in color:
        fhandle.write(str(int(col))+'\n')

