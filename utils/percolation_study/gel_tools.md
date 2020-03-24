# Gel Observalbes and their meaning 

### Measures (static=s, dynamic=d)


#### Most relevant 

* Final degree distribution (s)
* Final cluster size distribution (s)
* Size of largest cluster (s,d)
* Average degrees (s,d)
* Fraction of occupied bonds pb(t) (s,d)
* fraction of p/np bonds (s,d)

#### Less relevant (but maybe useful later)
* Final domain size distribution of p/np domains (s)
* Largest domain size distribution of p/np domains (s,d)


### Why are these measures relevant

##### Final degree distribution (s)
The final degree distribution varies between 0 and 4 (integer values)
A node of degree=1 is an endpoint, a node of degree 2 is a chain-point, and nodes of degree 3 and 4 are branching points.
But note that many nodes of degree 4 also indicate a transistion from a gel-like to a crystal-like assembly. 

##### Final cluster size distribuiton (s)
The cluster size distribution indicates how far along in the percolation transition the system is. 
If there is only one cluster left, the system is definitely percolated. 

Question to be answered: what if the cluster size is exponentially etc distributed? 

##### Size of largest cluster (s,d)
The size of the largest cluster as function of time indicates how far along in the peroclation transition the system is.
A widely fluctating largest cluster size at after a long simulation runs, indicates closeness to percolation. 

##### Average degrees (s,d)
The growth of the average degrees with time indicates how the gel evolves.

##### Fraction of occupied bonds (s,d)
The growth of the pb with time indicates how the gel evolves. If pb does not reach a plateau, this indicates that the gel 
is stuck/ not an euquilibrium gel anymore. 

##### Fraction of p/np bonds (s,d)

The fraction of p/np bonds shows the local ordering of the gel. For example at phi=0.5 and T=0.15 we see many np-bonds, that indicate the formation of the boxes that are building blocks of the T=0.15 gel, while at T=0.1, we observe many parallel bonds, inicative of the local harmonica ordering. 

TBC
