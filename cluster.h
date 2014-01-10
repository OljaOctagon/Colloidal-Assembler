#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include <stdbool.h>

typedef struct {

	bool	isactive;
	int		clusterid;

} particle_cluster_info;

// particle cluster info für jedes Particle/Würfel -> *op_cluster info

extern particle_cluster_info *op_cluster_info;

//op_cluster_build weist jedem Cluster eine id zu
void op_cluster_build(void);
//initialisiert op_cluster info
void op_cluster_init(double cutoff1ss);

void op_cluster_free(void);
int op_cluster_get_largest_cluster_size(void);

#endif
