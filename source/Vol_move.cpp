#include "pmove.h"

    void pmove::Iso_Vol_Change (particles& Particles, box* Box, fileio& Fileio, int mc_time){
		
    	double dU_old, dU_new;	
		double delta_U;
		double Total_Halo_Energy;
	
		dU_old=0;
		dU_new=0;

		rand_V =  gsl_rng_uniform(r01);
		rand_V = dmax_V*(rand_V-0.5);
		Box->V_old = Box->V;

		Box->Lx_old = Box->Lx;
		Box->Ly_old = Box->Ly;
		Box->Lz_old = Box->Lz;
		
		Box->Lx = Box->Lx + ((Box->Lx)/(Box->Lx + Box->Ly))*rand_V;
		Box->Ly = Box->Ly + ((Box->Ly)/(Box->Lx + Box->Ly))*rand_V;

		Box->V = Box->Lx*Box->Ly*Box->Lz;
		Box->V_rel = double(Box->V)/double(Box->V_old);

		Box->Lx_scale = Box->Lx/Box->Lx_old;
		Box->Ly_scale = Box->Ly/Box->Ly_old;
		Box->Lz_scale = Box->Lz/Box->Lz_old;
		  
        b_factor_pre = exp(-1.0*Box->P_sigma*(Box->V - Box->V_old) + Box->N*log(Box->V_rel));
		b_factor = minimum(1,b_factor_pre);
		  
		XI = gsl_rng_uniform(r01);
		  
	    //b_factor<XI: Trial move not accepted
		if(b_factor<XI){
				
		     Box->V = Box->V_old;		
		     Box->Lx = Box->Lx_old;
		     Box->Ly = Box->Ly_old;
			 Box->Lz = Box->Lz_old;
			 //cout<<"no favorable volume"<<endl;
		
		}
		
		if(b_factor>=XI){
			
			
				Box->edges_from_center();	
			
				for(int id=0;id<Box->N;id++){  
					  
					trans_vec.x = Box->x_center + (Particles.N_Particle[id]->x_center - Box->x_center)*Box->Lx_scale - Particles.N_Particle[id]->x_center;
					trans_vec.y = Box->y_center + (Particles.N_Particle[id]->y_center - Box->y_center)*Box->Ly_scale - Particles.N_Particle[id]->y_center;
					//trans_vec.z = Box->z_center + (Particles.N_Particle[id]->z_center - Box->z_center)*Box->Lz_scale - Particles.N_Particle[id]->z_center;
						
					Trans_Update_Positions(Particles, id, trans_vec);	
					
					Particles.Check_Periodic_CM(id, Box);
				    Update_Periodic_Positions(Particles, Box, id);	
					
					
				}    

				exit_status = 0;
				col_count = 0;
				id = -1;
				
	     		do{
				
				    id = id +1;
					
				    Particles.Collision_List[id].Calculate_OP(Box, id, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);
					Collision_Test( Particles, Box, id, Particles.Collision_List);
					Total_Halo_Energy += Halo_Energy;

											
				    } while((exit_status==0)&&(id!=Box->N-1));
						
					
				 id=0;
				

				if(exit_status>=1){

					 //cout<<"overlaps"<<endl;
				
					for(int id=0;id<Box->N;id++){				
						Reset_Positions(Particles, id);
						
					}	
							
							
					Box->V = Box->V_old;
					
					Box->Lx = Box->Lx_old;
					Box->Ly = Box->Ly_old;
					Box->Lz = Box->Lz_old;
								
					Box->edges_from_center();				
				}	  
	

				 if(exit_status == 0){

				 	for (id=0;id<Box->N;id++){
				 		Particles.N_Particle[id]->Calculate_Axis();
						Particles.N_Particle[id]->Calculate_Patch_Position();

					}

		 			dU_old = Particles.Total_Energy;

		 			Calculate_Pair_Potential( Particles, Box);		
					dU_new = Particles.Total_Energy;
					delta_U = dU_new - dU_old;

				 	//MC Step for Potential
				 	
				 	b_factor_pre = exp(-1.0*beta*delta_U);
					b_factor = minimum(1,b_factor_pre);
				  
					XI = gsl_rng_uniform(r01);
				 
				 	  if (XI > b_factor){
				 	  	 //cout<<"no favorable energy"<<endl;
				
				
						for(int id=0;id<Box->N;id++){				
							Reset_Positions(Particles, id);
							Particles.N_Particle[id]->Calculate_Axis();
						 	Particles.N_Particle[id]->Calculate_Patch_Position();
						}	
								
						Box->V = Box->V_old;
						
						Box->Lx = Box->Lx_old;
						Box->Ly = Box->Ly_old;
						Box->Lz = Box->Lz_old;
									
						Box->edges_from_center();		
						Particles.Total_Energy = dU_old;

				 	  }

				 	  // Accept move if XI <= b_factor 

				 	if (XI <= b_factor){

				 		//cout<<"accepted"<<endl;

				 	 	for(int id=0;id<Box->N;id++){	
						Set_Positions(Particles, id);
						}
					    
						Box->packing_fraction = Particles.N_Particle[0]->V*double(Box->N)/Box->V;
				 	}


				}
		 
			}	
		}	
						
			
