
# import network x 
import network as nx
import numpy # we need it later 

# input: list of all connected pairs
# output: particles and size of largest cluster


# that's how you initialize a graph. It's empty now
 G = nx.Graph()

# you can add all edges at once by just passing a list 'edge_list'
G.add_edges_from(edge_list)

# edge_list is a list of all the bonded particles
# Format is
# particle i particle j
#  int        int
# you need to generate it beforehand. 


# this gives you all connected domains as node-list, the data type is generator ( Generators are common data constructs in python).

domains = nx.connected_components(G)

# to transform it to a list, do list(domains). Now we have a list of lists of # the nodes in all domains.
domain_list =list(domains)

# use list comprehension to get the length of the domains 
domain_length = np.array([ len(domain) for domain in domain_list ])

# get the id of the largest domain 
d_id = np.argmax(domain_length)

# get the  particles in the largest domain  (That's what you want)
particles_biggest_cluster = np.array(domain_list[d_id])

# get the size of the largest size of largest cluster
N_largest = len(domain_list[d_id])


# tada :D 
