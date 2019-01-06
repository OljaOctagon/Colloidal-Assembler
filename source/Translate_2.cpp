#include "pmove.h"

    void pmove::Translate(particles& Particles, box* Box, fileio& Fileio, int id, int mc_time){
	
		 rand_x = gsl_rng_uniform(r01);
		 rand_y = gsl_rng_uniform(r01);
		 rand_z = gsl_rng_uniform(r01);

		 trans_vec.x = dmax_t*(rand_x - 0.5);
		 trans_vec.y = dmax_t*(rand_y - 0.5);
		 
		if (is_2D==1){
			trans_vec.z = 0.0; 
		}
		 
		 else{	 
			trans_vec.z = dmax_t*(rand_z - 0.5);
		}	
		 //Update position
		 
		Trans_Update_Positions(Particles, id, trans_vec);	
		 
		 // is center of mass outside of the box?
		Particles.Check_Periodic_CM(id, Box);

		Particles.Update_Cell_List(id, Box);
		 
		Update_Periodic_Positions(Particles, Box, id);
				
		Particles.Collision_List[id].Calculate(Box, id, Particles.Id_Cell_List, Particles.Cell_List, Particles.Cell, Particles.N_Particle, Particles.N_Particle[0].cut_off);	
		 
		
			
			
		 exit_status = 0;
		 col_count = 0;
		
		 Collision_Test( Particles, Box, id, Particles.Collision_List);
	
	     if(exit_status>=1){
	
			 Reset_Positions(Particles, id);
			 
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
		
			 
