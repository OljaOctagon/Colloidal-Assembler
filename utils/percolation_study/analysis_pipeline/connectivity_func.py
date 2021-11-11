import numpy as np
import pandas as pd
import argparse
import gel_tools as gt
import networkx as nx
import glob
import os
from itertools import islice
from itertools import cycle 


def get_connectivity(connections):

    connections = connections[:,:2]
    G=nx.Graph()
    G.add_edges_from(connections)
    
    degree_sequence = [d for n,d in G.degree()]

    average_degree = np.mean(degree_sequence)
    np.hist(degree_sequence)

    # save average degree with run_id together with spanning 
    # save degree in pandas dataframe with run_id and 0,1,2,3,4 as columns ready for plotting 

    largest_cc = max(nx.connected_components(G), key=len)
    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]
   
    # get edge connectivity 
    n_to_disconnect = nx.edge_connectivity(S)

    # take larget connected component and take that one 


    communities = label_propagation_communities(G_starWars)
    print([community for community in communities])
