#include <cstdio>
#include <gsl/gsl_math.h>
#include "cluster.h"
#include "nl.h"

particle_cluster_info *op_cluster_info;

/* local variables */
static int cnums;
static double rcut_1ss;

/* local functions*/
static void HarvestClusterRecursion(int clusternumber,int ipart);

static void HarvestClusterRecursion(int clusternumber,int ipart){

	int i,j;
 
 //ipart = index des momentanen Teilchens
 // nl_neigh (in nl.h) das Struct für die Nachbarliste
 // nl_neigh.nneigh = Anzahl der Nachbarn dieses Teilchens
 // nl_neigh.ids[i] = der tatsächliche Index des iten Nachbars des Teilchens
 // nl_neigh.dnorm[i] Abstand zum iten Nachbarn
 
 
 
 
	for (i=0; i<nl_neigh[ipart].nneigh; i++){
		
		j = nl_neigh[ipart].ids[i];

		// only use neighbors of 1 solvation shell
		if (nl_neigh[ipart].dnorm[i] <= rcut_1ss){
			if ( op_cluster_info[j].isactive == true && op_cluster_info[j].clusterid == -1){
				op_cluster_info[j].clusterid = clusternumber;
				HarvestClusterRecursion(clusternumber, j);
			}	
		}
	}

}

void op_cluster_build(void){

	int i;
	int clusternumber=0;
		
	for (i=0; i<nl_natoms; i++){
		if (op_cluster_info[i].isactive == true && op_cluster_info[i].clusterid == -1){
			clusternumber++;
			op_cluster_info[i].clusterid = clusternumber;
			HarvestClusterRecursion(clusternumber, i);
		}
	}
	
	cnums = clusternumber;

}

void op_cluster_init(double cutoff1ss){

	int i;
	rcut_1ss = cutoff1ss;
	op_cluster_info = (particle_cluster_info*)malloc(nl_natoms*sizeof(particle_cluster_info));
	cnums = 0;
	for(i=0; i<nl_natoms; i++) {
		op_cluster_info[i].isactive = false;
        op_cluster_info[i].clusterid = -1;
    }

}

void op_cluster_free(void){

	free(op_cluster_info);

}

int op_cluster_get_largest_cluster_size(void){

	int *a;
	int i;
	int ret=0;
	
	a = (int*)calloc(cnums, sizeof(int));
	
	for (i=0; i<nl_natoms; i++) {
        if (op_cluster_info[i].isactive) {
	    	a[op_cluster_info[i].clusterid-1]++;
        }
    }
	
	for (i=0; i<cnums; i++){
		ret = GSL_MAX(ret, a[i]);
	}
	
	free(a);
	
	return ret;

}


/*
  
   .) ClusterListe init und setzen
  
  
  
  
    .) // op = BCtWm
    if(calc_BCtWm == 1) {
        op_cluster_init(cutoff);
        for(int i=0; i<natom; i++) {
            if(tWm[i] == 1.0) {    // Kriterum für Cubus fest/flüssig...
                op_cluster_info[i].isactive = true; /if fest
            }
            else {
                op_cluster_info[i].isactive = false; /if flüssig
            }
        }
      .)  op_cluster_build();
        //BCtWm = double(op_cluster_get_largest_cluster_size()); // Groesse des Grossten Clusters
        BCtWm = op_cluster_get_largest_cluster_size_weighted();
        for(int i=0; i<natom; i++) {
            tWmCID[i] = op_cluster_info[i].clusterid;
        }
       .) op_cluster_free();
    }


*/
