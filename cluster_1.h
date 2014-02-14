#include <cstdio>
#include <gsl/gsl_math.h>
#include "order_parameter_2.h"
#include "particles.h"
#include "box.h"


		
	class particle_cluster_info{

	public:

	bool isactive;
	int  cluster_id;
	
	
	};

	class cluster{

	public:

	particle_cluster_info *op_cluster_info;

	cluster(box* Box);
	~cluster();
	
	int Size_R;
	int cnums;
	int list_j;
	
	void HarvestClusterRecursion(int cluster_number, int id_j, particles& Particles);
	void op_cluster_build(box* Box, particles& Particles);
	void op_cluster_free();
	int  op_cluster_get_largest_cluster_size(int cnums, box* Box);
	void Calculate(order_parameter& Order_Parameter, particles& Particles, box* Box);
	void Reset_cluster(box* Box);
	
    };





