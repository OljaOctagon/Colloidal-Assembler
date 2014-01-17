#include "move.h"

    void move::Translate(particles& Particles, box* Box, fileio& Fileio, int id, int mc_time){
	
		
		 rand_x = gsl_ran_gaussian(r01, sigma_trans);
		 rand_y = gsl_ran_gaussian(r01, sigma_trans);
		 rand_z = gsl_ran_gaussian(r01, sigma_trans);
		 
		 trans_vec.x = rand_x;
		 trans_vec.y = rand_y;
		 trans_vec.z = rand_z;

		 //Update position
		 
		Trans_Update_Positions(Particles, id, trans_vec);	
		 
		 // is center of mass outside of the box?
		Particles.Check_Periodic_CM(id, Box);

		Update_Periodic_Positions(Particles, Box, id);	


		Particles.N_Particle[id].Calculate_Axis();


		Particles.Update_Cell_List(id, Box);
		
				
		Particles.Collision_List[id].Calculate(Box, id, Particles.Id_Cell_List, Particles.Cell_List, Particles.Cell, Particles.N_Particle, Particles.N_Particle[0].cut_off, Particles.MAX_coll_p);	
		
		//TEST 
		/*
		int M;
		ofstream xyz_out("pos.xyz", ios::out | ios::app);
		M=Box->N*Particles.N_Particle[0].edge_N;		
		xyz_out<<M<<endl;
		xyz_out<<"Particles of frame "<<mc_time<<endl;
		
		
		for(int k=0;k<Box->N;k++){
		
				for(int j=0;j<Particles.N_Particle[k].edge_N;j++){
					
					 xyz_out<<"cb"<<"       "<<Particles.N_Particle[k].x[j]<<"   "<<Particles.N_Particle[k].y[j]<<"   "<<Particles.N_Particle[k].z[j]<<endl;
					 
					}   
			}
			
		xyz_out.close();
		*/
		//TEST	
		
	
			
		 exit_status = 0;
		 col_count = 0;
		
		 Collision_Test( Particles, Box, id, Particles.Collision_List);
	
	     if(exit_status>=1){
	
			 Reset_Positions(Particles, id);
			 
			 Particles.N_Particle[id].Calculate_Axis();
			 
			 Particles.Reset_Cell_List(Box, id, Particles.c_id, Particles.n_id, Particles.id_num);
			

	 	    }	
			
			
			
		 if (exit_status == 0){
			
			 Set_Positions(Particles, id);
	
			 accept_translate = accept_translate +1;
		 
		    }

	

		 Particles.N_Particle[id].trans_periodic[0] = 0.0;
		 Particles.N_Particle[id].trans_periodic[1] = 0.0;
		 Particles.N_Particle[id].trans_periodic[2] = 0.0;
		
		
		 N_trans_moves = N_trans_moves + 1;
 
		
    } 
		
			 
