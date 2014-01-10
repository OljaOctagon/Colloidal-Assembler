
#include "move.h"



	double move::Calculate_Cross_Product_x(m_vector a, m_vector b){
		double cx;
		cx = a.y*b.z - a.z*b.y;			
		return cx;	
		
	}
	
	double move::Calculate_Cross_Product_y(m_vector a, m_vector b){	
		double cy;
		cy = a.z*b.x - a.x*b.z;
		return cy;	
	}	
	
	double move::Calculate_Cross_Product_z(m_vector a, m_vector b){
		double cz;
		cz = a.x*b.y - a.y*b.x;
		return cz;	
	
	}	


	int move::Calculate_Separating_Axis_GENERAL(particles& Particles, int id, int j, m_vector* L_axis){
		
	
		
		int N_indep, N_cross;
		double norm_L;
		int N_separating_axis;
		
		N_indep = Particles.N_Particle[id].N_independent_faces;
		N_cross = Particles.N_Particle[id].N_cross_edges;
		
		N_separating_axis =  2*N_indep + N_cross*N_cross;
		//cout<<"N_separating_axis: "<<N_separating_axis<<endl;
		
		for (int k=0; k< N_indep; k++){
		
			L_axis[k].x = Particles.N_Particle[id].facenormal[k].x;
			L_axis[k].y = Particles.N_Particle[id].facenormal[k].y;
			L_axis[k].z = Particles.N_Particle[id].facenormal[k].z;
			
		}	
	
		for(int k=N_indep; k<(2*N_indep);k++){
			
			L_axis[k].x = Particles.N_Particle[j].facenormal[k-N_indep].x;
			L_axis[k].y = Particles.N_Particle[j].facenormal[k-N_indep].y;
			L_axis[k].z = Particles.N_Particle[j].facenormal[k-N_indep].z;
			
		
		}
		
		int k;
		
		k = 2*N_indep;
		
		for(int s=0; s<N_cross;s++){
			for(int m=0;m<N_cross;m++){
			
			
				L_axis[k].x = Calculate_Cross_Product_x(Particles.N_Particle[id].edges[s], Particles.N_Particle[j].edges[m]);
				L_axis[k].y = Calculate_Cross_Product_y(Particles.N_Particle[id].edges[s], Particles.N_Particle[j].edges[m]);
				L_axis[k].z = Calculate_Cross_Product_z(Particles.N_Particle[id].edges[s], Particles.N_Particle[j].edges[m]);
				
				norm_L = L_axis[k].norm();
				
				L_axis[k].x = L_axis[k].x/norm_L;
				L_axis[k].y = L_axis[k].y/norm_L;
				L_axis[k].z = L_axis[k].z/norm_L;
				
				
				k=k+1;
				
			}
		}
		
	   return N_separating_axis;
			
			
	}

/*
	void move::Calculate_Separating_Axis( particles& Particles, int id, int j, m_vector* L_axis ){
				
		double norm_L;		
				
				
		 L_axis[0].x = Particles.N_Particle[id].ax_1.x;
		 L_axis[0].y = Particles.N_Particle[id].ax_1.y;
		 L_axis[0].z = Particles.N_Particle[id].ax_1.z;

			
		 L_axis[1].x = Particles.N_Particle[id].ax_2.x;
		 L_axis[1].y = Particles.N_Particle[id].ax_2.y;
		 L_axis[1].z = Particles.N_Particle[id].ax_2.z;

         L_axis[2].x = Particles.N_Particle[id].ax_3.x;
		 L_axis[2].y = Particles.N_Particle[id].ax_3.y;
		 L_axis[2].z = Particles.N_Particle[id].ax_3.z;


		 L_axis[3].x = Particles.N_Particle[j].ax_1.x;
		 L_axis[3].y = Particles.N_Particle[j].ax_1.y;
		 L_axis[3].z = Particles.N_Particle[j].ax_1.z;


		 L_axis[4].x = Particles.N_Particle[j].ax_2.x;
		 L_axis[4].y = Particles.N_Particle[j].ax_2.y;
		 L_axis[4].z = Particles.N_Particle[j].ax_2.z;


		 L_axis[5].x = Particles.N_Particle[j].ax_3.x;
		 L_axis[5].y = Particles.N_Particle[j].ax_3.y;
		 L_axis[5].z = Particles.N_Particle[j].ax_3.z;


		//Cross-Product of L_axis 0 - 5 with each other

		 L_axis[6].x = L_axis[0].y*L_axis[3].z - L_axis[0].z*L_axis[3].y;
		 L_axis[6].y = L_axis[0].z*L_axis[3].x - L_axis[0].x*L_axis[3].z;
		 L_axis[6].z = L_axis[0].x*L_axis[3].y - L_axis[0].y*L_axis[3].x;
			
		
		 norm_L = L_axis[6].norm();
			
		 L_axis[6].x = L_axis[6].x/norm_L;
		 L_axis[6].y = L_axis[6].y/norm_L;
		 L_axis[6].z = L_axis[6].z/norm_L;
			
			
		 L_axis[7].x = L_axis[0].y*L_axis[4].z - L_axis[0].z*L_axis[4].y;
		 L_axis[7].y = L_axis[0].z*L_axis[4].x - L_axis[0].x*L_axis[4].z;
		 L_axis[7].z = L_axis[0].x*L_axis[4].y - L_axis[0].y*L_axis[4].x;


		 norm_L = L_axis[7].norm();	
		 
		 L_axis[7].x = L_axis[7].x/norm_L;
		 L_axis[7].y = L_axis[7].y/norm_L;
		 L_axis[7].z = L_axis[7].z/norm_L;



		 L_axis[8].x = L_axis[0].y*L_axis[5].z - L_axis[0].z*L_axis[5].y;
		 L_axis[8].y = L_axis[0].z*L_axis[5].x - L_axis[0].x*L_axis[5].z;
		 L_axis[8].z = L_axis[0].x*L_axis[5].y - L_axis[0].y*L_axis[5].x;

		 norm_L = L_axis[8].norm();	


		 L_axis[8].x = L_axis[8].x/norm_L;
		 L_axis[8].y = L_axis[8].y/norm_L;
		 L_axis[8].z = L_axis[8].z/norm_L;


		 L_axis[9].x = L_axis[1].y*L_axis[3].z - L_axis[1].z*L_axis[3].y;
		 L_axis[9].y = L_axis[1].z*L_axis[3].x - L_axis[1].x*L_axis[3].z;
		 L_axis[9].z = L_axis[1].x*L_axis[3].y - L_axis[1].y*L_axis[3].x;


		norm_L = L_axis[9].norm();	

		 L_axis[9].x = L_axis[9].x/norm_L;
		 L_axis[9].y = L_axis[9].y/norm_L;
		 L_axis[9].z = L_axis[9].z/norm_L;


		 L_axis[10].x = L_axis[1].y*L_axis[4].z - L_axis[1].z*L_axis[4].y;
		 L_axis[10].y = L_axis[1].z*L_axis[4].x - L_axis[1].x*L_axis[4].z;
		 L_axis[10].z = L_axis[1].x*L_axis[4].y - L_axis[1].y*L_axis[4].x;


		norm_L = L_axis[10].norm();	
			
		 L_axis[10].x = L_axis[10].x/norm_L;
		 L_axis[10].y = L_axis[10].y/norm_L;
		 L_axis[10].z = L_axis[10].z/norm_L;



		 L_axis[11].x = L_axis[1].y*L_axis[5].z - L_axis[1].z*L_axis[5].y;
		 L_axis[11].y = L_axis[1].z*L_axis[5].x - L_axis[1].x*L_axis[5].z;
		 L_axis[11].z = L_axis[1].x*L_axis[5].y - L_axis[1].y*L_axis[5].x;

		norm_L = L_axis[11].norm();	
			
		 L_axis[11].x = L_axis[11].x/norm_L;
		 L_axis[11].y = L_axis[11].y/norm_L;
		 L_axis[11].z = L_axis[11].z/norm_L;


		 L_axis[12].x = L_axis[2].y*L_axis[3].z - L_axis[2].z*L_axis[3].y;
		 L_axis[12].y = L_axis[2].z*L_axis[3].x - L_axis[2].x*L_axis[3].z;
		 L_axis[12].z = L_axis[2].x*L_axis[3].y - L_axis[2].y*L_axis[3].x;

			
		norm_L = L_axis[12].norm();		
			
		 L_axis[12].x = L_axis[12].x/norm_L;
		 L_axis[12].y = L_axis[12].y/norm_L;
		 L_axis[12].z = L_axis[12].z/norm_L;



		 L_axis[13].x = L_axis[2].y*L_axis[4].z - L_axis[2].z*L_axis[4].y;
		 L_axis[13].y = L_axis[2].z*L_axis[4].x - L_axis[2].x*L_axis[4].z;
		 L_axis[13].z = L_axis[2].x*L_axis[4].y - L_axis[2].y*L_axis[4].x;

        norm_L = L_axis[13].norm();	

		 L_axis[13].x = L_axis[13].x/norm_L;
		 L_axis[13].y = L_axis[13].y/norm_L;
		 L_axis[13].z = L_axis[13].z/norm_L;


		 L_axis[14].x = L_axis[2].y*L_axis[5].z - L_axis[2].z*L_axis[5].y;
		 L_axis[14].y = L_axis[2].z*L_axis[5].x - L_axis[2].x*L_axis[5].z;
		 L_axis[14].z = L_axis[2].x*L_axis[5].y - L_axis[2].y*L_axis[5].x;

		norm_L = L_axis[14].norm();	

		 L_axis[14].x = L_axis[14].x/norm_L;
		 L_axis[14].y = L_axis[14].y/norm_L;
		 L_axis[14].z = L_axis[14].z/norm_L;



	}  

*/

	double move::Scalar_Product(m_vector a, m_vector b){

		 m_scalar_product = a.x*b.x + a.y*b.y + a.z*b.z;	
         return m_scalar_product;	
		
        }  
         
         
    
    void move::Collision_Test( particles& Particles, box* Box, int id, collision_list* Collision_List){
		
		int cp_id = 0;
		int N_it;
		N_it = 0;	
			
		//cout<<"beginning exit_status: "<<exit_status<<endl;	
		Particles.N_Particle[id].Calculate_Face_Normals();
			
		do{
					
			 collision_partner = Collision_List[id].Elements[cp_id].nl_id;	
					
			 if(collision_partner!=id){
					
					
				 	
				 Particles.N_Particle[collision_partner].Calculate_Face_Normals();	
					
			     N_it = Calculate_Separating_Axis_GENERAL( Particles, id, collision_partner, L_axis);
			   
                
                 int j;	
				 j=-1;	 

				 f_exit_MAX = -3.0;
				 f_exit = 0.0;
							
				 do{	

						
				     j= j+1;
				     
				     
				     //cout<<"L_axis.x  "<<j<<" "<<L_axis[j].x<<endl;
					 //cout<<"L_axis.y  "<<j<<" "<<L_axis[j].y<<endl;
					 //cout<<"L_axis.z  "<<j<<" "<<L_axis[j].z<<endl;
					 //cout<<"L_axis.norm(): "<<j<<" "<<L_axis[j].norm()<<endl;
				     
				     
					 L_axis_norm2 = Scalar_Product(L_axis[j],L_axis[j]);
							
					 if(L_axis_norm2>0.5){
								
					     dL = Scalar_Product(Collision_List[id].Elements[cp_id].distance, L_axis[j]);	
						 dL = fabs(dL);
					
					
						 R_id = Particles.N_Particle[id].Calculate_Projection_to_Separating_Axis(L_axis[j]);
						 R_collision_partner = Particles.N_Particle[collision_partner].Calculate_Projection_to_Separating_Axis(L_axis[j]);
						 	
								 
				
						 f_exit = dL - R_id - R_collision_partner;
                             
                         if(f_exit>f_exit_MAX){
							 f_exit_MAX = f_exit;
							}
									
					    }	
							
							
					//}while((j<N_it)&&(f_exit_MAX<0.0));	
				}while((j!=43)&&(f_exit_MAX<0.0));
				
				 //cout<<"f_exit_MAX:"<<f_exit_MAX<<endl;						
								
			     if(f_exit_MAX<0.0){
									
				     col_count = 1;
				     
				     //cout<<"collision! particle "<<id<< " with particle "<<collision_partner<<endl;
				     
								
					}
								
			     if(f_exit_MAX>=0.0){
											
					 //cout<<"no collsioon... Separating_axis: "<<j<<endl;										
					 col_count = 0;
							 
				    }
							
				 exit_status = exit_status + col_count;
						
				
			    }
				
				cp_id=cp_id+1; 
				//cout<<"exit_status: "<<exit_status<<endl; 
				
				//cout<<"Collision_List[id].Nm "<<id<<" "<<Collision_List[id].Nm<<endl;
				
				 

		} while((exit_status==0)&&(cp_id<Collision_List[id].Nm));	
				
		
	}
     
	
	
	/*
	void move::Collision_Test( particles& Particles, box* Box, int id, collision_list* Collision_List){
		
		int cp_id = 0;
			
		//cout<<"beginning exit_status: "<<exit_status<<endl;	
			
		do{
					
			 collision_partner = Collision_List[id].Elements[cp_id].nl_id;	
					
			 if(collision_partner!=id){
					
			     Calculate_Separating_Axis( Particles, id, collision_partner, L_axis);
                 
         
                 int j;	
				 j=-1;	 

				 f_exit_MAX = -3.0;
				 f_exit = 0.0;
							
				 do{	

				     j= j+1;
					 L_axis_norm2 = Scalar_Product(L_axis[j],L_axis[j]);
							
					 if(L_axis_norm2!=0){
								
					     dL = Scalar_Product(Collision_List[id].Elements[cp_id].distance, L_axis[j]);	
						 dL = fabs(dL);
								

						 R_id = (fabs(Scalar_Product(Particles.N_Particle[id].ax_1, L_axis[j])) + fabs(Scalar_Product(Particles.N_Particle[id].ax_2, L_axis[j])) + fabs(Scalar_Product(Particles.N_Particle[id].ax_3, L_axis[j])))/2.0;
						 R_collision_partner = (fabs(Scalar_Product(Particles.N_Particle[collision_partner].ax_1, L_axis[j])) + fabs(Scalar_Product(Particles.N_Particle[collision_partner].ax_2, L_axis[j])) + fabs(Scalar_Product(Particles.N_Particle[collision_partner].ax_3, L_axis[j])))/2.0;
						 f_exit = dL - R_id - R_collision_partner;
                             
                         if(f_exit>f_exit_MAX){
								
							 f_exit_MAX = f_exit;
									
										
							}
									
							
					    }	
							
							
					}while((j!=14)&&(f_exit_MAX<0.0));	
					
				
								
			     if(f_exit_MAX<0.0){
									
				     col_count = 1;
				     
				     //cout<<"collision! particle "<<id<< " with particle "<<collision_partner<<endl;
				     
								
					}
								
			     if(f_exit_MAX>=0.0){
											
					 //cout<<"no collsioon... Separating_axis: "<<j<<endl;										
					 col_count = 0;
							 
				    }
							
				 exit_status = exit_status + col_count;
						
				
			    }
				
				cp_id=cp_id+1; 
				//cout<<"exit_status: "<<exit_status<<endl; 
				 

		}while((exit_status==0)&&(cp_id<Collision_List[id].Nm));	
				
		
		
		
		
	}
	*/
