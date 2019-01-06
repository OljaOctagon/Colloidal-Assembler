
#include "main.h"

 
int main(int argc, char * argv[]){
	
		
		     
		r = gsl_rng_alloc(gsl_rng_ranlxd2);
		r01 = gsl_rng_alloc(gsl_rng_ranlxd2); 	
			
		strcpy(runtype, argv[1]);
		
		strcpy(new_runtype, "-n");
		strcpy(former_runtype, "-f");
	
		checkpoint_time = atoi(argv[2]);
		start_time = checkpoint_time;
	
		term_time = atoi(argv[3]);
		
		fileio Fileio;
		Fileio.Set_Vars();
		Fileio.Read_Parameters();
		
		//term_time = Fileio.time_in;
		calc_frequency = Fileio.calc_frequency_in;
		frame_frequency = Fileio.frame_frequency_in;
		checkpoint_frequency = Fileio.checkpoint_frequency_in;
	
		box* Box;
		
		Box = new box;

		
		particles Particles(Fileio.number_of_cells_in, Fileio.N_in);
	
		
		// init new run
		if(strcmp(runtype, new_runtype)==0){

			//initialize box
			
			Box->Startconfig(Fileio.N_in, Fileio.P_sigma_in, Fileio.packing_fraction_in);
			
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
	    
		if(strcmp(runtype, former_runtype)==0){

			//initialize former box
			Fileio.Read_Box(Box, checkpoint_time);
			Box->Startconfig_former(Fileio.N_in, Fileio.P_sigma_in);
			Box->edges_from_center();
			
			//random state initialize
			Fileio.Read_Random_State(r,r01, checkpoint_time);

			// initialize former positions
			Particles.Startconfig(Box);
			Fileio.Read_Positions(Box, Particles, checkpoint_time);
			Fileio.Read_Orientations(Box, Particles, checkpoint_time);
			//Fileio.Read_Orientations_Axis_Angle(Box, Particles, checkpoint_time);
			//Fileio.Read_Orientations_EULER(Box, Particles.N_Particle);
			
			
			
			Particles.Set_former_Config(Box);
			
			Particles.Make_Cell_List(Box);
			Particles.Set_Cell_List(Box);
			Particles.Make_Cell_Neighbour_List(); 
			
		
	    
	    }
	    
	    
		// init moves
	    order_parameter Order_Parameter(Box, Particles);
		cluster Cluster(Box);
		
		pmove Move(Box->N, Particles.edge_N, Fileio.delta_tmax, Fileio.delta_amax, Fileio.delta_thetamax, Fileio.delta_Vmax, Fileio.delta_Lmax, Fileio.delta_alpha_max, Fileio.delta_beta_max, Fileio.delta_gamma_max);
		//Move.Calculate_Cluster_List(Particles, Box);
		Move.Set_Random_State(r,r01);
	
		int time;
		int MC_cycle_time;
		time = start_time;
	
		Box->P_sigma = Fileio.P_sigma_in;
	
		//Fileio.Save_Config(Box, Particles, Move.r, Move.r01, Move.dmax_t, Move.dmax_alpha, Move.dmax_beta, Move.dmax_gamma, Move.dmax_V, Move.dmax_L, time);
		 
		cout<<"Box.phi "<<Box->packing_fraction<<endl; 
		  
		  
		do{  
			
			time = time + 1;
			
			MC_cycle_time = 2*Box->N + 1;
		
		
			for(int cycle_time = 0; cycle_time<MC_cycle_time; cycle_time++){
				 	
				Move.Iterate(Particles, Box, Fileio, time);
				//cout<<"cycle_time "<<cycle_time<<endl;	
			}
		  
			if(time%calc_frequency==0){
			  
				//cout<<"time:"<<time<<endl;
			  
				//D2_order = Order_Parameter.D2(Box, Particles.N_Particle);
				//D4_order = Order_Parameter.D4(Box, Particles.N_Particle);
				
				double Cut_Off = 1.4;
				Order_Parameter.Calculate_Local_Order_Parameters(Particles, Particles.N_Particle, Box, Cut_Off);
				Order_Parameter.Caculate_is_nPhase(Box);
				Cluster.Calculate(Order_Parameter, Particles, Box);
				
				
				Fileio.Write_Order_Parameters(time, D2_order, D4_order);
				Fileio.Write_Local_Order_Parameters(time, Particles, Box, Order_Parameter, Cluster);
				//Fileio.Write_Cluster_List(Move.Cluster_List, Move.cluster_size, Particles, Order_Parameter);
				
				Move.Calculate_Acceptances(time);
			
				Fileio.Write_Acceptances(time, Move.accept_translate_procent, Move.accept_rotate_procent, Move.accept_iso_vol_procent, Move.accept_complete_procent);
				Fileio.Write_NPT(time, Box);
				
		    } 

			
			if(time%frame_frequency==0){
				
				Fileio.Write_Positions(time, Box, Particles);
				Fileio.Write_Box_Positions(time, Box);
				
			}
		
		
			if(time%checkpoint_frequency==0){
			
			Fileio.Save_Config(Box, Particles, Move.r, Move.r01, Move.dmax_t, Move.dmax_alpha, Move.dmax_beta, Move.dmax_gamma, Move.dmax_V, Move.dmax_L, time);
				
			}
		
		
		}while((time<term_time)&&(exit_code==false)); 
	
		
		Fileio.Save_Config(Box, Particles, Move.r, Move.r01, Move.dmax_t, Move.dmax_alpha, Move.dmax_beta, Move.dmax_gamma, Move.dmax_V, Move.dmax_L, term_time);
		//Fileio.Write_g_r(term_time,Order_Parameter);
		//Order_Parameter.write_histogram_info();
		
		cout<<"End of Simulation. BYE BYE"<<endl;
		
		
	
	}  

     	 
     
     	
		
