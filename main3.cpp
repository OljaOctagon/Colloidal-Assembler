#include "main3.h"

 
int main(int argc, char * argv[]){
	
	     
		r = gsl_rng_alloc(gsl_rng_ranlxd2);
		r01 = gsl_rng_alloc(gsl_rng_ranlxd2); 	
			
		strcpy(runtype, argv[1]);
		
		strcpy(new_runtype, "-n");
		strcpy(former_runtype, "-f");
		strcpy(tis_runtype, "-tis");
		strcpy(flux_runtype, "-flux");
		
		
		checkpoint_time = atoi(argv[2]);
		start_time = checkpoint_time;
	
		term_time = atoi(argv[3]);
		
		wv_0_in = 0;
		wv_i_in = 1;
		wv_ip1_in = 2;
		
		
		fileio Fileio;
		Fileio.Set_Vars();
		Fileio.Read_Parameters();
		
		
		calc_frequency = Fileio.calc_frequency_in;
		frame_frequency = Fileio.frame_frequency_in;
		checkpoint_frequency = Fileio.checkpoint_frequency_in;
	
		box* Box;
		Box = new box;

		
		strcpy(traj_type_in, "forward");
		
		particles Particles(Fileio.number_of_cells_in, Fileio.N_in, Fileio.MAX_coll_p_in, Fileio.MAX_fshell_p_in);
	
		
		// init new run
		if(strcmp(runtype, new_runtype)==0){

			//initialize box
			
			Box->Startconfig(Fileio.N_in, Fileio.P_sigma_in, Fileio.packing_fraction_in, Particles.N_Particle[0].V);
			
			start_time = 0;
			
			Box->edges_from_center();

			// initialize lattice
			Particles.Startconfig(Box);
			
			
			Particles.Make_Cell_List(Box);
			Particles.Set_Cell_List(Box);
			Particles.Make_Cell_Neighbour_List(); 
		
		
			gsl_rng_set(r01,0xf143);
			gsl_rng_set(r, 0x01a23);
			
			
			
	    }
	    
	    int time;
	    
	    time = start_time;
	   
				
		move Move(Box->N, Particles.edge_N, Fileio.delta_tmax, Fileio.delta_rmax, Fileio.delta_Vmax, Fileio.is_translation_ON, Fileio.is_rotation_ON, Fileio.is_volumemove_ON);
	
		Move.Set_Random_State(r,r01);
	
		
		int MC_cycle_time;
		time = start_time;
		 				
		Fileio.Write_Positions(time, Box, Particles);
		Fileio.Write_Box_Positions(time, Box);
				
		
	
		Box->P_sigma = Fileio.P_sigma_in;
		
		cout<<"Packing_fraction: "<<Box->packing_fraction<<endl;
		
	
		  
		do{  
			
			time = time + 1;
			cout<<"time:"<<time<<endl;
			
			MC_cycle_time = Move.mt_sum;
			
			//cout<<"MC_cycle_time:" <<MC_cycle_time<<endl;
			
		
			for(int cycle_time = 0; cycle_time<MC_cycle_time; cycle_time++){
				//cout<<"Move.Interate"<<endl; 			 	
				Move.Iterate(Particles, Box, Fileio, time);
				
					
			}
				
				
			if(time%frame_frequency==0){
				
				Fileio.Write_Positions(time, Box, Particles);
				Fileio.Write_Box_Positions(time, Box);
				
			}
			
		
			if((time%checkpoint_frequency==0)||(time == start_time +1 )){
				
				Fileio.Save_Config(Box, Particles, Move.r, Move.r01, Move.dmax_t, Move.dmax_alpha, Move.dmax_beta, Move.dmax_gamma, Move.dmax_V, Move.dmax_L, time); 
				N_check_points = N_check_points + 1;
				
				
			}
				
				
			
			if((exit_code==true)&&(time%checkpoint_frequency==0)){  
				term_code=true; 
			}	  
		

		
		}while((time<term_time)&&(exit_code==false)); 
		
		
		Fileio.Save_Config(Box, Particles, Move.r, Move.r01, Move.dmax_t, Move.dmax_alpha, Move.dmax_beta, Move.dmax_gamma, Move.dmax_V, Move.dmax_L, term_time);
		
		cout<<"End of Simulation. BYE BYE"<<endl;
		
		
	    
	    
}
