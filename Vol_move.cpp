#include "move.h"

    void move::Iso_Vol_Change (particles& Particles, box* Box, fileio& Fileio, int mc_time){

		Particles.Set_Cell_List(Box);	
		 
		rand_V =  gsl_rng_uniform(r01);
		 
		Box->V_old = Box->V;
		Box->V = Box->V_old + dmax_V*(rand_V-0.5);
		
		Box->V_rel = double(Box->V)/double(Box->V_old);
		Box->VL_rel = pow(Box->V_rel,(1./3.));
		 
		Box->Lx_old = Box->Lx;
		Box->Ly_old = Box->Ly;
		Box->Lz_old = Box->Lz;
		 
		Box->Lx = Box->Lx*Box->VL_rel;
		Box->Ly = Box->Ly*Box->VL_rel;
		Box->Lz = Box->Lz*Box->VL_rel;
		  
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
		
		}
			
		
			//if Box.V greater than Box.V_old: move accepted, no collision test needed
		if(b_factor>=XI){
			
			
			if(Box->V > Box->V_old){	
			
			
			
				Particles.Cell[0].V = Box->V/double(Particles.number_of_cells);
					 
				Particles.Cell[0].Lx = pow(Particles.Cell[0].V,1./3.);
				Particles.Cell[0].Ly = Particles.Cell[0].Lx;
				Particles.Cell[0].Lz = Particles.Cell[0].Lx;
	
	
				//cout<<"Particles.Cell[0].Lx"<<Particles.Cell[0].Lx<<endl;	
				
				Box->edges_from_center();
				
				for(int c_id = 0; c_id<Particles.number_of_cells; c_id++){
						
					Particles.Cell[c_id].x_center = Box->x[0] + Particles.Cell[0].Lx/2.0 + Particles.Cell[0].Lx*double(c_id%Particles.N_c);
					Particles.Cell[c_id].y_center = Box->y[0] + Particles.Cell[0].Ly/2.0 + Particles.Cell[0].Ly*double((c_id/Particles.N_c)%Particles.N_c);	
					Particles.Cell[c_id].z_center = Box->z[0] + Particles.Cell[0].Lz/2.0 + Particles.Cell[0].Lz*double(c_id/(Particles.N_c*Particles.N_c));				
					Particles.Cell[c_id].edges_from_center();
				
				}
			  
				Box->packing_fraction = Particles.N_Particle[0].V*double(Box->N)/Box->V;
					
					
				for(int id=0;id<Box->N;id++){  
					  
					trans_vec.x = Box->x_center + (Particles.N_Particle[id].x_center - Box->x_center)*Box->Lx_scale - Particles.N_Particle[id].x_center;
					trans_vec.y = Box->y_center + (Particles.N_Particle[id].y_center - Box->y_center)*Box->Ly_scale - Particles.N_Particle[id].y_center;
					trans_vec.z = Box->z_center + (Particles.N_Particle[id].z_center - Box->z_center)*Box->Lz_scale - Particles.N_Particle[id].z_center;
						
					Trans_Update_Positions(Particles, id, trans_vec);	
					
					//Particles.Check_Periodic_CM(id, Box);
					//Update_Periodic_Positions(Particles, Box, id);	
					
					Set_Positions(Particles, id); 
						 
						
				}    
				
				Particles.Make_Cell_List(Box);  
				  
				accept_iso_vol = accept_iso_vol +1;
						
					
					
					
			}	  
	
				
		    // Box.V smaller than Box.V_old: Collision Test needed
		     
			if(Box->V<=Box->V_old){
			
					
			    //Save old Cell Volumina and Lengths
				Particles.Cell_old[0].V = Particles.Cell[0].V;
						
				Particles.Cell_old[0].Lx = Particles.Cell[0].Lx; 
				Particles.Cell_old[0].Ly = Particles.Cell[0].Ly;
				Particles.Cell_old[0].Lz = Particles.Cell[0].Lz;
						
				//Calculate new Cell Volumina and Lengths
				Particles.Cell[0].V = Box->V/double(Particles.number_of_cells);
			
				Particles.Cell[0].Lx = pow(Particles.Cell[0].V,1./3.);
				Particles.Cell[0].Ly = Particles.Cell[0].Lx;
				Particles.Cell[0].Lz = Particles.Cell[0].Lx;
				
  		        Box->edges_from_center();
						
				for(int c_id = 0; c_id<Particles.number_of_cells; c_id++){
											
					Particles.Cell[c_id].x_center = Box->x[0] + Particles.Cell[0].Lx/2.0 + Particles.Cell[0].Lx*double(c_id%Particles.N_c);
					Particles.Cell[c_id].y_center = Box->y[0] + Particles.Cell[0].Ly/2.0 + Particles.Cell[0].Ly*double((c_id/Particles.N_c)%Particles.N_c);	
					Particles.Cell[c_id].z_center = Box->z[0] + Particles.Cell[0].Lz/2.0 + Particles.Cell[0].Lz*double(c_id/(Particles.N_c*Particles.N_c));				
					Particles.Cell[c_id].edges_from_center();
				
				}
		  
				for(id=0;id<Box->N; id++){		
						
					trans_vec.x = Box->x_center + (Particles.N_Particle[id].x_center - Box->x_center)*Box->Lx_scale - Particles.N_Particle[id].x_center;
					trans_vec.y = Box->y_center + (Particles.N_Particle[id].y_center - Box->y_center)*Box->Ly_scale - Particles.N_Particle[id].y_center;
					trans_vec.z = Box->z_center + (Particles.N_Particle[id].z_center - Box->z_center)*Box->Lz_scale - Particles.N_Particle[id].z_center;
						
			        Trans_Update_Positions(Particles, id, trans_vec);	
						
					//Check if center of mass is outside of the box
					Particles.Check_Periodic_CM(id, Box);

					Update_Periodic_Positions(Particles, Box, id);
							
				}	
					
				Particles.Make_Cell_List(Box);
				 
				exit_status = 0;
				col_count = 0;
				id = -1;
				
	     		do{
				
				    id = id +1;
							
					Particles.Collision_List[id].Calculate(Box, id, Particles.Id_Cell_List, Particles.Cell_List, Particles.Cell, Particles.N_Particle, Particles.N_Particle[0].cut_off, Particles.MAX_coll_p);
					Collision_Test( Particles, Box, id, Particles.Collision_List);	
											
				    } while((exit_status==0)&&(id!=Particles.max_id));
						
					
				 id=0;
				
				
									
				if(exit_status>=1){
				
					for(int id=0;id<Box->N;id++){
										
						Reset_Positions(Particles, id);
												
					}	
							
							
					Box->V = Box->V_old;
					
					Box->Lx = Box->Lx_old;
					Box->Ly = Box->Ly_old;
					Box->Lz = Box->Lz_old;
								
					Box->edges_from_center();
								
					Particles.Cell[0].V = Particles.Cell_old[0].V;
						
					Particles.Cell[0].Lx = Particles.Cell_old[0].Lx; 
					Particles.Cell[0].Ly = Particles.Cell_old[0].Ly;
					Particles.Cell[0].Lz = Particles.Cell_old[0].Lz;
								
					Box->edges_from_center();
					
							
					for(int c_id = 0; c_id<Particles.number_of_cells; c_id++){
								
							Particles.Cell[c_id].x_center = Box->x[0] + Particles.Cell[0].Lx/2.0 + Particles.Cell[0].Lx*double(c_id%Particles.N_c);
							Particles.Cell[c_id].y_center = Box->y[0] + Particles.Cell[0].Ly/2.0 + Particles.Cell[0].Ly*double((c_id/Particles.N_c)%Particles.N_c);	
							Particles.Cell[c_id].z_center = Box->z[0] + Particles.Cell[0].Lz/2.0 + Particles.Cell[0].Lz*double(c_id/(Particles.N_c*Particles.N_c));				
							Particles.Cell[c_id].edges_from_center();
						
						    }
				  
					Particles.Reset_Cell_List(Box);
							
								
				}
				
				
			    if(exit_status == 0){
						
					for(int id=0;id<Box->N;id++){	
						Set_Positions(Particles, id);
							
					}
							
				    accept_iso_vol = accept_iso_vol +1;
			        Box->packing_fraction = Particles.N_Particle[0].V*double(Box->N)/Box->V;
						    
				}
							
						
            }
                
		}


        N_iso_vol_moves = N_iso_vol_moves + 1;
            
    }        
        
        
    void move::Aniso_Vol_Change (particles& Particles, box* Box, fileio& Fileio, int mc_time){
	
	     Particles.Set_Cell_List(Box);	
				
		 rand_Lx =  gsl_rng_uniform(r01);
		 rand_Ly =  gsl_rng_uniform(r01);
		 rand_Lz =  gsl_rng_uniform(r01);
			
		 Box->Lx_old = Box->Lx;
		 Box->Ly_old = Box->Ly;
		 Box->Lz_old = Box->Lz;
				
		 Box->V_old = Box->V;
				
			
		 Box->Lx = Box->Lx + dmax_L*(rand_Lx - 0.5);
		 Box->Ly = Box->Ly + dmax_L*(rand_Ly - 0.5);
		 Box->Lz = Box->Lz + dmax_L*(rand_Lz - 0.5);

			
		 // Trial move: Scale Volume isotropically
		 Box->V = Box->Lx*Box->Ly*Box->Lz;
		 Box->V_rel = Box->V/Box->V_old;
		 Box->VL_rel = pow(Box->V_rel,(1./3.));
				  
		 Box->Lx_scale = Box->Lx/Box->Lx_old;
		 Box->Ly_scale = Box->Ly/Box->Ly_old;
		 Box->Lz_scale = Box->Lz/Box->Lz_old;	
		
		 b_factor_pre = exp(-1.0*Box->P_sigma*(Box->V - Box->V_old) + Box->N*log(Box->V_rel));
		 b_factor = minimum(1,b_factor_pre);
		
		 XI = gsl_rng_uniform(r01);
			  
		 Box->edges_from_center();
		
		 if(b_factor<XI){
					
			 Box->V = Box->V_old;
						
			 Box->Lx = Box->Lx_old;
			 Box->Ly = Box->Ly_old;
			 Box->Lz = Box->Lz_old;
				
					
		    }
				
			  
		 if(b_factor>=XI){
			  
			 if((Box->Lx>Box->Lx_old) && (Box->Ly> Box->Ly_old) && (Box->Lz>Box->Lz_old)){
				
			     Particles.Cell[0].Lx =  Particles.Cell[0].Lx*Box->Lx_scale;
				 Particles.Cell[0].Ly =  Particles.Cell[0].Ly*Box->Ly_scale;
				 Particles.Cell[0].Lz =  Particles.Cell[0].Lz*Box->Lz_scale;
				
				 Particles.Cell[0].V = Particles.Cell[0].Lx*Particles.Cell[0].Ly*Particles.Cell[0].Lz;   
			     Particles.Make_Cell_List(Box);
										
			     for(int id=0;id<Box->N;id++){  
					  
					 trans_vec.x = Box->x_center + (Particles.N_Particle[id].x_center - Box->x_center)*Box->Lx_scale - Particles.N_Particle[id].x_center;
					 trans_vec.y = Box->y_center + (Particles.N_Particle[id].y_center - Box->y_center)*Box->Ly_scale - Particles.N_Particle[id].y_center;
					 trans_vec.z = Box->z_center + (Particles.N_Particle[id].z_center - Box->z_center)*Box->Lz_scale - Particles.N_Particle[id].z_center;
						
					 Trans_Update_Positions(Particles, id, trans_vec);	
							
					 Set_Positions(Particles, id);
						
					}    
				  

				}	  
			
		
			
			 else{
					
				 for(id=0;id<Box->N; id++){		
							
				     trans_vec.x = Box->x_center + (Particles.N_Particle[id].x_center - Box->x_center)*Box->Lx_scale - Particles.N_Particle[id].x_center;
				     trans_vec.y = Box->y_center + (Particles.N_Particle[id].y_center - Box->y_center)*Box->Ly_scale - Particles.N_Particle[id].y_center;
				     trans_vec.z = Box->z_center + (Particles.N_Particle[id].z_center - Box->z_center)*Box->Lz_scale - Particles.N_Particle[id].z_center;
							
				     Trans_Update_Positions(Particles, id, trans_vec);	

				     // Check if center of mass is outside of the box
					 Particles.Check_Periodic_CM(id, Box);

					 Update_Periodic_Positions(Particles, Box, id);
						 
					}	
					
							
				 Particles.Cell[0].Lx =  Particles.Cell[0].Lx*Box->Lx_scale;
				 Particles.Cell[0].Ly =  Particles.Cell[0].Ly*Box->Ly_scale;
				 Particles.Cell[0].Lz =  Particles.Cell[0].Lz*Box->Lz_scale;
								
				 Particles.Cell[0].V = Particles.Cell[0].Lx*Particles.Cell[0].Ly*Particles.Cell[0].Lz;   
								
				 Particles.Make_Cell_List(Box);
										
				 exit_status = 0;
				 col_count = 0;
				 id = -1;
					
				 do {
					
					 id = id +1;
							
					 Particles.Collision_List[id].Calculate(Box, id, Particles.Id_Cell_List, Particles.Cell_List, Particles.Cell, Particles.N_Particle, Particles.N_Particle[0].cut_off, Particles.MAX_coll_p);
					 Collision_Test( Particles, Box, id, Particles.Collision_List);		
						
													
				    } while((exit_status==0) || (id!=Particles.max_id));
							
							
				 id = 0;
			
				 if(exit_status>=1){
							
				     Reset_Positions(Particles, id);
									
				     Particles.Reset_Cell_List(Box);
			
				     Box->V = Box->V_old;
							
				     Box->Lx = Box->Lx_old;
				     Box->Ly = Box->Ly_old;
				     Box->Lz = Box->Lz_old;
									
				     Particles.Cell[0].V = Particles.Cell_old[0].V;
							
				     Particles.Cell[0].Lx = Particles.Cell_old[0].Lx; 
				     Particles.Cell[0].Ly = Particles.Cell_old[0].Ly;
				     Particles.Cell[0].Lz = Particles.Cell_old[0].Lz;
									
								
					}
							
				 if(exit_status ==0){
								
					 Set_Positions(Particles, id);
								
					}
							
				}
			}

        }        
        
        
        
        
  
