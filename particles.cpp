# include "particles.h"


	particles::~particles(){
		
			
	}

    particles::particles(){
		  
		N_Particle = 0;
		N_Particle_old = 0;
		Cell = 0;
		Cell_old = 0;
		 
		  
    }


	particles::particles(int number_of_cells_in, int size, int MAX_coll_p_in, int MAX_fshell_in){
  
	     //Allocation of vars		
         number_of_cells = number_of_cells_in;
		 max_id = size-1;
		 MAX_coll_p = MAX_coll_p_in; 
		 MAX_fshell_p = MAX_fshell_p;	
		
		
		 //N_Particle = new cube[size];
		 //N_Particle_old = new cube[size];
		 
		 
		 //N_Particle = new octahedron[size];
		 //N_Particle_old = new octahedron[size];
		 
		 N_Particle = new truncated_cube[size];
		 N_Particle_old = new truncated_cube[size];
		 
		 
		 
		
		 Cell = new cell[number_of_cells];
		 Cell_old = new cell[number_of_cells];


		//Initialization of N_Particle	
  
		for(int id =0;id<size;id++){
	 
			
			
			N_Particle[id].Set_Initial_Quaternion();
			N_Particle_old[id].Set_Initial_Quaternion();
			
		
			N_Particle[id].Set_Initial_Axis();
			N_Particle_old[id].Set_Initial_Axis();
			
			
		    N_Particle[id].Set_Lengths();
		    N_Particle_old[id].Set_Lengths();
		    
		    cout<<"N_Particle[0].V: "<<N_Particle[0].V<<endl;
		    cout<<"N_Particle[0].Lx: "<<N_Particle[0].Lx<<endl;
		    cout<<"N_Particle[0].cut_off: "<<N_Particle[0].cut_off<<endl;
			
			
		}	
			
	}
   
	
    void particles::Startconfig(box* Box){
 
		
		 N_p = rint(pow(double(Box->N),1./3.));
		 N_p_float = double(N_p);  
		  
		 l_diff = (Box->Lx - N_p_float)/N_p_float;   
		 Box->N = int(Box->N); 	

		 Cell[0].V = Box->V/double(number_of_cells);
		 Cell[0].Lx = pow(Cell[0].V,1./3.);
		 Cell[0].Ly = Cell[0].Lx;
		 Cell[0].Lz = Cell[0].Lx; 
	   
	     Cell_old[0].V = Box->V/double(number_of_cells);
		 Cell_old[0].Lx = pow(Cell[0].V,1./3.);
		 Cell_old[0].Ly = Cell[0].Lx;
		 Cell_old[0].Lz = Cell[0].Lx; 
	    
	     for(int c_id = 1;c_id<number_of_cells;c_id++){
	    
	        
			 Cell[c_id].V = Cell[0].V;
			 Cell[c_id].Lx = Cell[0].Lx;
			 Cell[c_id].Ly = Cell[0].Ly;
			 Cell[c_id].Lz = Cell[0].Lz; 
			
			 Cell_old[c_id].V = Cell[0].V;
			 Cell_old[c_id].Lx = Cell[0].Lx;
			 Cell_old[c_id].Ly = Cell[0].Ly;
			 Cell_old[c_id].Lz = Cell[0].Lz; 
			
	        }
		
		MAX_cell_members = 5*int(ceil(Cell[0].V)/double(N_Particle[0].V));
		
 
		Collision_List = new collision_list[Box->N];
		for(int id=0; id<Box->N; id++){
				Collision_List[id].Set(MAX_coll_p);
		}
		
		
	
		Cell_List = new int*[number_of_cells];
		for(int j=0;j<number_of_cells;j++){
			 Cell_List[j] = new int[MAX_cell_members]; 
			}   
			
		 Cell_List_old = new int*[number_of_cells];
		 for(int j=0;j<number_of_cells;j++){
			 Cell_List_old[j] = new int[MAX_cell_members]; 
			}   	
			    
		 Id_Cell_List = new int[Box->N]; 
		
		 Id_Cell_List_old = new int[Box->N];


	     for (int id=0; id<Box->N; id++){

			 N_Particle[id].Set_Start_Lattice_Position(id, Box->Lx, Box->N);
			 
			 /*		  
			 N_Particle[id].x_center = 0.5 +  l_diff/2.0 + double(id %N_p + l_diff*(id % N_p));
			 N_Particle[id].y_center = 0.5 +  l_diff/2.0 + double((id/N_p)%N_p + l_diff*((id/N_p)%N_p)); 
			 N_Particle[id].z_center = 0.5 +  l_diff/2.0 + double(id/int(pow((double)N_p,2)) + l_diff*(id/int((double)pow(N_p,2))));
		     */
		 
			 N_Particle_old[id].x_center = N_Particle[id].x_center;
			 N_Particle_old[id].y_center = N_Particle[id].y_center;
			 N_Particle_old[id].z_center = N_Particle[id].z_center;
			 
			 N_Particle[id].copy_count = 0;
			
			 N_Particle[id].edges_from_center();
			 N_Particle_old[id].edges_from_center();
			
				 
			 N_Particle[id].trans_periodic[0] = 0.0;
			 N_Particle[id].trans_periodic[1] = 0.0;
			 N_Particle[id].trans_periodic[2] = 0.0;
			
			
			 //calculate ax vectors
		
			 N_Particle[id].Calculate_Axis();
			 
			 	
			 N_Particle[id].cm_left_count = 0;
			 N_Particle[id].cm_right_count = 0;
			
			 N_Particle[id].cm_top_count = 0;
			 N_Particle[id].cm_bottom_count = 0;
			
			 N_Particle[id].cm_back_count = 0;
			 N_Particle[id].cm_front_count = 0;
			
			 N_Particle[id].cm_out = 0;
			
				 
	        }
	    
		 //Make_Cell_List(Box, Cell, N_Particle);
		 //Set_Cell_List(Box, Cell, Cell_old, N_Particle);
		 //Make_Cell_Neighbour_List(Cell); 

    }   

	void particles::Set_former_Config(box* Box){
		
			for(int id=0;id<Box->N;id++){
				
			
			 Rot_mat_INIT[0] = 1.0 - 2.0*N_Particle[id].q.y*N_Particle[id].q.y - 2.0*N_Particle[id].q.z*N_Particle[id].q.z;
			 Rot_mat_INIT[1] = 2.0*N_Particle[id].q.x*N_Particle[id].q.y - 2.0*N_Particle[id].q.w*N_Particle[id].q.z;
			 Rot_mat_INIT[2] = 2.0*N_Particle[id].q.x*N_Particle[id].q.z + 2.0*N_Particle[id].q.w*N_Particle[id].q.y;
				
			 Rot_mat_INIT[3] = 2.0*N_Particle[id].q.x*N_Particle[id].q.y + 2.0*N_Particle[id].q.w*N_Particle[id].q.z;
			 Rot_mat_INIT[4] = 1.0 - 2.0*N_Particle[id].q.x*N_Particle[id].q.x - 2.0*N_Particle[id].q.z*N_Particle[id].q.z;
			 Rot_mat_INIT[5] = 2.0*N_Particle[id].q.y*N_Particle[id].q.z - 2.0*N_Particle[id].q.w*N_Particle[id].q.x;
				
			 Rot_mat_INIT[6] = 2.0*N_Particle[id].q.x*N_Particle[id].q.z - 2.0*N_Particle[id].q.w*N_Particle[id].q.y;
			 Rot_mat_INIT[7] = 2.0*N_Particle[id].q.y*N_Particle[id].q.z + 2.0*N_Particle[id].q.w*N_Particle[id].q.x;
			 Rot_mat_INIT[8] = 1.0 - 2.0*N_Particle[id].q.x*N_Particle[id].q.x - 2.0*N_Particle[id].q.y*N_Particle[id].q.y;
			 
			
			
				
			 N_Particle[id].edges_from_center();
			 N_Particle[id].distance_from_center();
		   
			 for(int j=0;j<N_Particle[id].edge_N;j++){
			  
				 N_Particle[id].new_dist_x[j] = Rot_mat_INIT[0]*N_Particle[id].dist_x[j] +Rot_mat_INIT[1]*N_Particle[id].dist_y[j] + Rot_mat_INIT[2]*N_Particle[id].dist_z[j]; 
				 N_Particle[id].new_dist_y[j] = Rot_mat_INIT[3]*N_Particle[id].dist_x[j] +Rot_mat_INIT[4]*N_Particle[id].dist_y[j] + Rot_mat_INIT[5]*N_Particle[id].dist_z[j]; 
				 N_Particle[id].new_dist_z[j] = Rot_mat_INIT[6]*N_Particle[id].dist_x[j] +Rot_mat_INIT[7]*N_Particle[id].dist_y[j] + Rot_mat_INIT[8]*N_Particle[id].dist_z[j];
						
				 N_Particle[id].x[j] = N_Particle[id].x_center + N_Particle[id].new_dist_x[j]; 
				 N_Particle[id].y[j] = N_Particle[id].y_center + N_Particle[id].new_dist_y[j];
				 N_Particle[id].z[j] = N_Particle[id].z_center + N_Particle[id].new_dist_z[j]; 
						   
				}  
				
				
			 N_Particle_old[id].x_center = N_Particle[id].x_center;
			 N_Particle_old[id].y_center = N_Particle[id].y_center;
			 N_Particle_old[id].z_center = N_Particle[id].z_center;
				
			 for(int j=0;j<N_Particle[id].edge_N;j++){
					
				 N_Particle_old[id].x[j] = N_Particle[id].x[j];
				 N_Particle_old[id].y[j] = N_Particle[id].y[j];
				 N_Particle_old[id].z[j] = N_Particle[id].z[j];
					
				}
		
				
			N_Particle_old[id].q.x = N_Particle[id].q.x;
			N_Particle_old[id].q.y = N_Particle[id].q.y;
			N_Particle_old[id].q.z = N_Particle[id].q.z;
			N_Particle_old[id].q.w = N_Particle[id].q.w;
					
			
						
			N_Particle[id].Calculate_Axis();	
			
					

		}
		 
			
	}

	
    void particles::Check_Periodic_CM(int id, box* Box){
	 
		 N_Particle[id].cm_left_count = 0; 
		 N_Particle[id].cm_right_count = 0;
		 N_Particle[id].cm_front_count = 0;
		 N_Particle[id].cm_back_count = 0;
		 N_Particle[id].cm_top_count = 0;
		 N_Particle[id].cm_bottom_count = 0;
		
		 N_Particle[id].cm_out = 0;
		
	  
		 if(N_Particle[id].x_center > Box->x[1]){ 
		     N_Particle[id].cm_right_count = 1; 
		    }  

		 if(N_Particle[id].x_center < Box->x[0]){
		     N_Particle[id].cm_left_count = 1;
		    }

		 if(N_Particle[id].y_center > Box->y[3]){
		     N_Particle[id].cm_back_count = 1;
		    }  
		
		 if(N_Particle[id].y_center < Box->y[0]){ 
		     N_Particle[id].cm_front_count = 1;   
		    }
		
		 if(N_Particle[id].z_center > Box->z[4]){
		     N_Particle[id].cm_top_count = 1;  
		    }
		
		 if(N_Particle[id].z_center < Box->z[0]){
		     N_Particle[id].cm_bottom_count = 1; 
		    } 
	   
	     N_Particle[id].cm_out = N_Particle[id].cm_right_count + N_Particle[id].cm_left_count
						+ N_Particle[id].cm_back_count + N_Particle[id].cm_front_count
						+ N_Particle[id].cm_top_count + N_Particle[id].cm_bottom_count;
	  

        }


	void particles::Check_Periodic_BC(int id, box* Box){
		
		 N_Particle[id].left_count = 0;
		 N_Particle[id].right_count = 0;
		 N_Particle[id].front_count = 0;
		 N_Particle[id].back_count = 0;
		 N_Particle[id].top_count = 0;
		 N_Particle[id].bottom_count = 0;
		
		 N_Particle[id].sum_edge_out = 0;
  
		 for(int k=0;k<N_Particle[id].edge_N;k++){
			   N_Particle[id].edge_out[k] = 0;  
            }  
 
		 for(int k=0;k<N_Particle[id].edge_N;k++){
  
			 if(N_Particle[id].x[k] > Box->x[1]){
				 N_Particle[id].right_count = 1;
				 N_Particle[id].edge_out[k] = 1;
                }

			 if(N_Particle[id].x[k] < Box->x[0]){
				 N_Particle[id].left_count = 1;
				 N_Particle[id].edge_out[k] = 1;
				}

			 if(N_Particle[id].y[k] > Box->y[2]){
				 N_Particle[id].back_count = 1;
				 N_Particle[id].edge_out[k] = 1;
				}
    
   
			 if(N_Particle[id].y[k] < Box->y[0]){
				 N_Particle[id].front_count = 1;
				 N_Particle[id].edge_out[k] = 1;
				}
    
			 if(N_Particle[id].z[k] > Box->z[4]){
				 N_Particle[id].top_count = 1;
				 N_Particle[id].edge_out[k] = 1;
                }
    
			 if(N_Particle[id].z[k] < Box->z[0]){
				 N_Particle[id].bottom_count = 1;
				 N_Particle[id].edge_out[k] = 1;
                 }
   
			 N_Particle[id].sum_edge_out = N_Particle[id].sum_edge_out + N_Particle[id].edge_out[k];
  
    
            }  
  
 
 
        }
