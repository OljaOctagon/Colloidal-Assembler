


void move::change_to_fractional(particles& Particles, box* Box){
	 
	 double lx,ly,lz;
	 double eta;
	 double alpa_b, beta_b, gamma_b;
	 
	 lx = Box->Lx;
	 ly = Box->Ly;
	 lz = Box->Lz;
	 
	 /*alpha_b = Box->alpha;
	 beta_b  = Box->beta;
	 gamma_b = Box->gamma;
	 
	 eta=(cos(alpha_b)- cos(gamma_b)*cos(beta_b))/(sin(gamma_b));
	  
	 M[0] = lx;
	 M[1] = ly*cos(gamma_b);
	 M[2] = lz*cos(beta_b);
	 
	 M[3] = 0;
	 M[4] = ly*sin(gamma_b);
	 M[5] = lz*eta;
	 
	 M[6] = 0;
	 M[7] = 0;
	 M[8] = lz*sqrt(2-cos(beta_b)*cos(beta_b) - eta*eta);
	 */
	 
	
	 M[0] = Box.v1.x; 
	 M[1] = Box.v1.y;
	 M[2] = Box.v1.z;
	 
	 M[3] = Box.v2.x;
	 M[4] = Box.v2.y;
	 M[5] = Box.v2.z;
	 
	 M[6] = Box.v3.x;
	 M[7] = Box.v3.y;
	 M[8] = Box.v3.z;
	 
	
	 for (id=0;id<Box->N;id++){
		 
		 Particles.N_Particle[id].x_center_scaled = M[0]*Particles.N_Particles[id].x_center + M[1]*Particles.N_Particles[id].y_center + M[2]*Particles.N_Particles[id].z_center;
		 Particles.N_Particle[id].y_center_scaled = M[3]*Particles.N_Particles[id].x_center + M[4]*Particles.N_Particles[id].y_center + M[5]*Particles.N_Particles[id].z_center;
		 Particles.N_Particle[id].z_center_scaled = M[6]*Particles.N_Particles[id].x_center + M[7]*Particles.N_Particles[id].y_center + M[8]*Particles.N_Particles[id].z_center;
	 	
		 
	 }
	 
	 
 
	}	 


void move::change_to_cartesian(){
	
	 double lx,ly,lz;
	 double eta;
	 double alpa_b, beta_b, gamma_b;
	 
	 lx = Box->Lx;
	 ly = Box->Ly;
	 lz = Box->Lz;
	 
	 /*alpha_b = Box->alpha;
	 beta_b  = Box->beta;
	 gamma_b = Box->gamma;
	 
	 eta=(cos(alpha_b)- cos(gamma_b)*cos(beta_b))/(sin(gamma_b));
	  
	 double cot_g, csc_g;
	 
	 cot_g = cos(gamma_b)/sin(gamma_b);
	 csc_g = 1./sin(gamma_b)
	  
	 st=sqrt(1-eta*eta - cos(beta_b)*cos(beta_b); 
	  
	 M[0] = 1.0/lx;
	 M[1] = -cot_g/a;
	 M[2] = (-cos(beta_b) + eta*cot_(gamma_b))/(a*st);
	 
	 M[3] = 0;
	 M[4] = csc_g/b;
	 M[5] = (eta*csc_g)/(b*st);
	 
	 M[6] = 0;
	 M[7] = 0;
	 M[8] = 1./(c*st);
	 */
	
	  for (id=0;id<Box->N;id++){
		 
		 Particles.N_Particle[id].x_center = M[0]*Particles.N_Particles[id].x_center_scaled + M[1]*Particles.N_Particles[id].y_center_scaled + M[2]*Particles.N_Particles[id].z_center_scaled;
		 Particles.N_Particle[id].y_center = M[3]*Particles.N_Particles[id].x_center_scaled + M[4]*Particles.N_Particles[id].y_center_scaled + M[5]*Particles.N_Particles[id].z_center_scaled;
		 Particles.N_Particle[id].z_center = M[6]*Particles.N_Particles[id].x_center_scaled + M[7]*Particles.N_Particles[id].y_center_scaled + M[8]*Particles.N_Particles[id].z_center_scaled;
	 	
		
	 
	 }

	} 

 void move::Vol_deform (particles& Particles, box* Box, fileio& Fileio, int mc_time){
				
				
		//SAVE OLD COORDINATES
		
		
		//Particles.Set_Cell_List(Box);	
		Box->V_old = Box->V;
		 
		Box->Lx_old = Box->Lx;
		Box->Ly_old = Box->Ly;
		Box->Lz_old = Box->Lz;
		
		Box->alpha_old = Box->alpha;
		Box->gamma_old = Box->gamma; 
		
		Box->v1_old.x = Box->v1.x;
		Box->v1_old.y = Box->v1.y;
		Box->v1_old.z = Box->v1.z;
		
		Box->v2_old.x = Box->v2.x;
		Box->v2_old.y = Box->v2.y;
		Box->v2_old.z = Box->v2.z;
		
		Box->v3_old.x = Box->v3.x;
		Box->v3_old.y = Box->v3.y;
		Box->v3_old.z = Box->v3.z;
		
		

		//CHANGE BOX: ADD A SMALL PERTURBATION TO ONE OF THE COORDINATES OF ONE OF THE VECTORS
		
		
		rand_V=gsl_ran_flat (r01,0,9);
		rand_v=gsl_rng_uniform(r01);
		
		
		if (rand_V<3){
			if(rand_V<1){
				box_ax_1.x = box_ax1.x + dmax_v*(rand_v-0.5);
			}
			if(randV >=1 and rand_V<2){
				box_ax_1.y = box_ax1.y + dmax_v*(rand_v-0.5);
			}
			if(rand_V>=2){
				box_ax_1.z = box_ax1.z + dmax_v*(rand_v-0.5);
			}			
		 }
		 
		if (rand_V>=3 and rand_V <6)  {
			if(rand_V<4){
				box_ax_2.x = box_ax2.x + dmax_v*(rand_v-0.5);
			}
			if(randV >=4 and rand_V<5){
				box_ax_2.y = box_ax2.y + dmax_v*(rand_v-0.5);
			}
			if(rand_V>=5){
				box_ax_2.z = box_ax2.z + dmax_v*(rand_v-0.5);
			}			
			
		}
		
		if (rand_V>=6){	
			if(rand_V<7){
				box_ax_3.x = box_ax3.x + dmax_v*(rand_v-0.5);
			}
			if(randV >=7 and rand_V<8){
				box_ax_3.y = box_ax3.y + dmax_v*(rand_v-0.5);
			}
			if(rand_V>=8){
				box_ax_3.z = box_ax3.z + dmax_v*(rand_v-0.5);
			}			
		
		}
		
        Box->Set_Lengths(Box->box_ax1, Box->box_ax2, Box->box_ax3);
		
		Box->V_rel = double(Box->V)/double(Box->V_old);
		
		b_factor_pre = exp(-1.0*Box->P_sigma*(Box->V - Box->V_old) + Box->N*log(Box->V_rel));
		b_factor = minimum(1,b_factor_pre);
		  
		XI = gsl_rng_uniform(r01);
		  
	    
	    //b_factor<XI: Trial move not accepted
		if(b_factor<XI){
				
		     Box->V = Box->V_old;
					
		     Box->Lx = Box->Lx_old;
		     Box->Ly = Box->Ly_old;
			 Box->Lz = Box->Lz_old;
			 
			 Box->alpha = Box->alpha_old;
			 Box->beta = Box->beta_old;
			 
			 	
			Box->ax_1.x = Box->ax_1_old.x;
			Box->ax_1.y = Box->ax_1_old.y;
			Box->ax_1.z = Box->ax_1_old.z;
			 
			Box->ax_2.x = Box->ax_2_old.x;
			Box->ax_2.y = Box->ax_2_old.y;
			Box->ax_2.z = Box->ax_2_old.z;
			
			Box->ax_3.x = Box->ax_3_old.x;
			Box->ax_3.y = Box->ax_3_old.y;
			Box->ax_3.z = Box->ax_3_old.z; 
			 
		
		}
			
		if(b_factor>=XI{
		
			//CHANGE TO FRACTIONAL COORDINATES
			change_to_fractional();
						
						
			//do trickts (periodic images, lattice reduction....)

			//CHANGE BACK TO CATESIAN
			change_to_cartesian();

			// calculate possible overlaps
			
			exit_status = 0;
				col_count = 0;
				id = -1;
				
	     		do{
				
				    id = id +1;
					//cout<<"Before"<<endl;		
					Particles.Collision_List[id].Calculate(Box, id, Particles.Id_Cell_List, Particles.Cell_List, Particles.Cell, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);
					//cout<<"After"<<endl;	
					Collision_Test( Particles, Box, id, Particles.Collision_List);	
											
				    } while((exit_status==0)&&(id!=Particles.max_id));
						
					
				 id=0;
			
								
			if(exit_status>=1){
				
				for(int id=0;id<Box->N;id++){
					
					//Set Positions somewhere
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
						
						    Particles.Cell[c_id].Lx = Particles.Cell_old[0].Lx; 
							Particles.Cell[c_id].Ly = Particles.Cell_old[0].Ly;
							Particles.Cell[c_id].Lz = Particles.Cell_old[0].Lz;
						
								
							Particles.Cell[c_id].x_center = Box->x[0] + Particles.Cell[0].Lx/2.0 + Particles.Cell[0].Lx*double(c_id%Particles.N_c);
							Particles.Cell[c_id].y_center = Box->y[0] + Particles.Cell[0].Ly/2.0 + Particles.Cell[0].Ly*double((c_id/Particles.N_c)%Particles.N_c);	
							Particles.Cell[c_id].z_center = Box->z[0] + Particles.Cell[0].Lz/2.0 + Particles.Cell[0].Lz*double(c_id/(Particles.N_c*Particles.N_c));				
							Particles.Cell[c_id].edges_from_center();
						
						    }
				  
				  
				  
					Particles.Reset_Cell_List(Box);
					Particles.Set_Cell_List(Box);
					
				
								
				}
				
				
			    if(exit_status == 0){
						
					for(int id=0;id<Box->N;id++){	
						Set_Positions(Particles, id);
							
					}
					
						
					
									
				    accept_iso_vol = accept_iso_vol +1;
			        Box->packing_fraction = (Particles.N_Particle[0]->V*double(Box->N))/Box->V;
			        Particles.Set_Cell_List(Box);
			        
			        
						    
				}
											


		
		
		
			}

//change to fractional coordinates



//accept or reject move acoording to exit of overlap detection


// to back to real coordinates


//The end
