#include "cluster_1.h"



	cluster::cluster(box* Box){
		
		int Resize;

		Resize=2000;

		op_cluster_info = new particle_cluster_info[Resize];
		Largest_Cluster_List= new int[Resize];


		cnums = 0;
		for(int i=0; i<Resize; i++) {
			op_cluster_info[i].isactive = false;
			op_cluster_info[i].cluster_id = -1;
			
		}

	}
	
	void cluster::Reset_cluster(box* Box){
		
		cnums = 0;
		for(int i=0; i<Box->N; i++) {
			op_cluster_info[i].isactive = false;
			op_cluster_info[i].cluster_id = -1;
			Largest_Cluster_List[i] = -1;
		}
		
	}


	cluster::~cluster(){
		
		delete [] op_cluster_info;	
		delete [] Largest_Cluster_List;
		
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


	int cluster::op_cluster_get_largest_cluster_size(int cnums, box* Box, int time_s){

		int *a;
		int ret=0;
		int l_cid;
		l_cid=0;
		
		
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

	
		string out_name;
		out_name="All_Clusters.dat";

		ofstream r_out( out_name.c_str(), ios::out | ios::app);

		for(int j=0;j<cnums;j++){
			r_out<<time_s<<" "<<a[j]<<endl;	
			if(a[j]>ret){
				ret = a[j];
				l_cid=j;	
			}	
		}	


		if(cnums==0){
            r_out<<time_s<<" "<<0<<endl;

       	}

       	r_out.close();


       //	cout<<"Cluster List Calc ...."<<endl;

		int k=0;
		for (int id=0;id<Box->N;id++){
			
			if((op_cluster_info[id].cluster_id -1)==l_cid){
				Largest_Cluster_List[k]= id;
				
				k=k+1;
			
			}
		}	
		//cout<<"N: "<<k<<endl; 

		out_name="All_Clusters_info.dat";
		ofstream r1_out(out_name.c_str(), ios::out | ios::app);	
	
			for (int j=0;j<cnums;j++){
				r1_out<<time_s<<" "<<j<<" "<<a[j]<<endl;
				for (int id=0;id<Box->N;id++){
					if((op_cluster_info[id].cluster_id-1)==j){
						r1_out<<id<<endl;
					}	
		
				}	
			}

		r1_out.close();	

		
		delete [] a;

		return ret;	

	}


int cluster::op_cluster_get_largest_cluster_size(int cnums, box* Box, particles& Particles, int time_s, bool is_mass_calc){

		int *a;
		int ret=0;
		int l_cid;
		l_cid=0;
		
		
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

	
		string out_name;
		out_name="All_Clusters.dat";

		ofstream r_out( out_name.c_str(), ios::out | ios::app);

		for(int j=0;j<cnums;j++){
			r_out<<time_s<<" "<<a[j]<<endl;	
			if(a[j]>ret){
				ret = a[j];
				l_cid=j;	
			}	
		}	


		if(cnums==0){
            r_out<<time_s<<" "<<0<<endl;

       	}

       	r_out.close();


       //	cout<<"Cluster List Calc ...."<<endl;

		int k=0;
		for (int id=0;id<Box->N;id++){
			
			if((op_cluster_info[id].cluster_id -1)==l_cid){
				Largest_Cluster_List[k]= id;
				
				k=k+1;
			
			}
		}	
		//cout<<"N: "<<k<<endl; 

		out_name="All_Clusters_info.dat";
		ofstream r1_out(out_name.c_str(), ios::out | ios::app);	
	
			for (int j=0;j<cnums;j++){
				r1_out<<time_s<<" "<<j<<" "<<a[j]<<endl;
				for (int id=0;id<Box->N;id++){
					if((op_cluster_info[id].cluster_id-1)==j){
						r1_out<<id<<endl;
					}	
		
				}	
			}

		r1_out.close();	

		//Calculate center of mass and write them for 2D gr


		////////////////////////////////////////////////////////////////////////////////////////

		double theta_x, theta_y, theta_av_x, theta_av_y;
		double xp1, xp2, yp1, yp2;
		int N_List;
		m_vector center_mass;


		if (is_mass_calc==true){

			out_name="cluster_center_of_mass.dat";
			ofstream r1_out(out_name.c_str(), ios::out | ios::app);	
		
			int cnums2;
			cnums2=0;
			for (int j=0;j<cnums;j++){
				if (a[j]==6){
					cnums2+=1;

				}

			}	

			int k=0;

 	 	 	for (int j=0;j<cnums;j++){
 	 	 		xp1=0;
				xp2=0;
				yp1=0;
				yp2=0;
				N_List = a[j];
				if (a[j]==6){

					for (int id=0;id<Box->N;id++){
						if((op_cluster_info[id].cluster_id-1)==j){
							theta_x = (Particles.N_Particle[id]->x_center/Box->Lx)*2*M_PI;
	 	 	 				theta_y = (Particles.N_Particle[id]->y_center/Box->Ly)*2*M_PI;

	 	 	 				xp1 = xp1 + (Box->Lx/(2*M_PI))*cos(theta_x);
	 	 	 				xp2 = xp2 + (Box->Lx/(2*M_PI))*sin(theta_x);

	 	 	 				yp1 = yp1 + (Box->Ly/(2*M_PI))*cos(theta_y);
	 	 	 				yp2 = yp2 + (Box->Ly/(2*M_PI))*sin(theta_y);
	 	 	 			}	
					}

					xp1=xp1/double(N_List);
					xp2=xp2/double(N_List);
					yp1=yp1/double(N_List);
					yp2=yp2/double(N_List);
					N_List=int(N_List);

					theta_av_x =atan2(-xp2, -xp1) + M_PI;
					theta_av_y =atan2(-yp2, -yp1) + M_PI;

					center_mass.x = (Box->Lx/(2*M_PI))*(theta_av_x);
					center_mass.y = (Box->Ly/(2*M_PI))*(theta_av_y);
					center_mass.z = Box->Lz/2.0;	

					r1_out<<time_s<<" "<<cnums2<<" "<<Box->Lx<<" "<<Box->Ly<<" "<<k<<" "<<center_mass.x<<" "<<center_mass.y<<endl;

					k=k+1;

				}	
			}	

		}

		/////////////////////////////////////////////////////////////////////////////////////////

		//cout<<"ret: "<<ret<<endl;

		delete [] a;

		return ret;	

	}







void cluster::Calculate(order_parameter& Order_Parameter, particles& Particles, box* Box, int time_s){

		cout<<"hello cluster calculate"<<endl;
	
		Reset_cluster(Box);

		for(int id=0; id<Box->N; id++) {
            if(Order_Parameter.is_nPhase[id] == 1.0) {    // Kriterum für Cubus fest/flüssig
                op_cluster_info[id].isactive = true; //if fest
            }
            else {
                op_cluster_info[id].isactive = false; //if flüssig
            }
        }

		op_cluster_build(Box, Particles);
      
        Size_R = op_cluster_get_largest_cluster_size(cnums, Box, time_s);

      
        cout <<"Size_R: "<<Size_R<<endl;
       
    }

void cluster::Calculate(order_parameter& Order_Parameter, particles& Particles, box* Box, int time_s, bool is_mass_calc){


		cout<<"hello cluster calculate"<<endl;
	
		Reset_cluster(Box);

		for(int id=0; id<Box->N; id++) {
            if(Order_Parameter.is_nPhase[id] == 1.0) {    // Kriterum für Cubus fest/flüssig
                op_cluster_info[id].isactive = true; //if fest
            }
            else {
                op_cluster_info[id].isactive = false; //if flüssig
            }
        }

		op_cluster_build(Box, Particles);
      
        Size_R = op_cluster_get_largest_cluster_size(cnums, Box, Particles, time_s, is_mass_calc);

        cout <<"Size_R: "<<Size_R<<endl;
       
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













