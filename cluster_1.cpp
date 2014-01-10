#include "cluster_1.h"



	cluster::cluster(box* Box){
		
		op_cluster_info = new particle_cluster_info[Box->N];

		cnums = 0;
		for(int i=0; i<Box->N; i++) {
			op_cluster_info[i].isactive = false;
			op_cluster_info[i].cluster_id = -1;
		}

	}
	
	void cluster::Reset_cluster(box* Box){
		
		cnums = 0;
		for(int i=0; i<Box->N; i++) {
			op_cluster_info[i].isactive = false;
			op_cluster_info[i].cluster_id = -1;
		}
		
	}


	cluster::~cluster(){
		
		delete [] op_cluster_info;	
		
	}


	void cluster::HarvestClusterRecursion(int cluster_number, int id_j, particles& Particles){

		for(int j=0;j<Particles.Collision_List[id_j].Nm;j++){
			 
				if(Particles.Collision_List[id_j].Elements[j].is_neighbour == 1){
									
				list_j = Particles.Collision_List[id_j].Elements[j].nl_id;
			
			    if( op_cluster_info[list_j].isactive == true && op_cluster_info[list_j].cluster_id == -1){
				op_cluster_info[list_j].cluster_id = cluster_number;
				HarvestClusterRecursion(cluster_number, list_j, Particles);
		
				}


			}

		}	
    }



	void cluster::op_cluster_build(box* Box, particles& Particles){

		int cluster_number=0;

		for(int id=0;id<Box->N;id++){
		   if((op_cluster_info[id].isactive == true) && (op_cluster_info[id].cluster_id == -1 )){
				cluster_number++;
				op_cluster_info[id].cluster_id = cluster_number;
				HarvestClusterRecursion(cluster_number,id, Particles);
			}
		}

		cnums=cluster_number;

	}


	int cluster::op_cluster_get_largest_cluster_size(int cnums, box* Box){

		int *a;
		int ret=0;
		
		
		
		a = new int[cnums];
		for(int j=0;j<cnums;j++){
		   a[j]=0;	
		}	
		
		int l = 0;	
		
		//cout<<"op_cluster_info[id].cluster_id: "<<op_cluster_info[4].cluster_id<<endl;	
			
		for(int id=0;id<Box->N;id++){
			if(op_cluster_info[id].isactive){
				
			   l = l+1;
			   a[op_cluster_info[id].cluster_id-1] = a[op_cluster_info[id].cluster_id-1] + 1;
			}

		
		}

		//cout<<"l: "<<l<<endl;

		for(int j=0;j<cnums;j++){
			
			ret=GSL_MAX(ret,a[j]);	
			
		}	
		//cout<<"ret: "<<ret<<endl;

		delete [] a;

		return ret;	

	}



 void cluster::Calculate(order_parameter& Order_Parameter, particles& Particles, box* Box){

	
		for(int id=0; id<Box->N; id++) {
            if(Order_Parameter.is_nPhase[id] == 1.0) {    // Kriterum für Cubus fest/flüssig
                op_cluster_info[id].isactive = true; //if fest
            }
            else {
                op_cluster_info[id].isactive = false; //if flüssig
            }
        }

		op_cluster_build(Box, Particles);
      
        Size_R = op_cluster_get_largest_cluster_size(cnums, Box);
        
        Reset_cluster(Box);
        
        //cout <<"Size_R: "<<Size_R<<endl;
       
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













