
#include "main2.h"

 
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
		
		//term_time = Fileio.time_in;
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
	    
		if(strcmp(runtype, former_runtype)==0){

			//initialize former box
			Fileio.Read_Box(Box, checkpoint_time);
			Box->Startconfig_former(Fileio.N_in, Fileio.P_sigma_in);
			Box->edges_from_center();
			
			//random state initialize
			Fileio.Read_Random_State(r,r01, checkpoint_time);
			
			//gsl_rng_set(r01,Fileio.seed1_in);
			//gsl_rng_set(r, Fileio.seed2_in);

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
	    
	    if(strcmp(runtype, tis_runtype)==0){

			
		
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
			
			Particles.Set_former_Config(Box);
			
			// initialize Cell List
			Particles.Make_Cell_List(Box);
			Particles.Set_Cell_List(Box);
			Particles.Make_Cell_Neighbour_List(); 
			
			
			wv_0_in = atoi(argv[4]);
			wv_i_in = atoi(argv[5]);
			wv_ip1_in = atoi(argv[6]);
			
			
			seed_r01 = atoi(argv[7]);
			seed_r = atoi(argv[8]);
			
			strcpy(traj_type_in, argv[9]);
			
			gsl_rng_set(r01,seed_r01);
			gsl_rng_set(r, seed_r);
			
			N_check_points = int(double(start_time)/double(checkpoint_frequency));		 
			
 	
		    
		
	    }
	   
	    if(strcmp(runtype, flux_runtype)==0){

			
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
			
			Particles.Set_former_Config(Box);
			
			// initialize Cell List
			Particles.Make_Cell_List(Box);
			Particles.Set_Cell_List(Box);
			Particles.Make_Cell_Neighbour_List(); 
			
			
			wv_0_in = atoi(argv[4]);
			wv_i_in = atoi(argv[5]);
			//wv_ip1_in = atoi(argv[6]);
			
			strcpy(traj_type_in, argv[6]);
			
		
			N_check_points = int(double(start_time)/double(checkpoint_frequency));		 
 	
		    //os.system("./McPoly -tis " + str(t_slice) +" "+ str(Nn_max_length) + " " + str(w_limit_1)+ " " + str(w_limit_2) + " " +  str(r_seed) + " " +str(r01_seed))
			
		
	    }
	    
	    
	    cout<<"traj_type: "<<traj_type_in<<endl;
	   
		tis_window Tis_Window(wv_0_in, wv_i_in, wv_ip1_in, term_time, traj_type_in);
	   
	   
		if((strcmp(runtype, tis_runtype)==0)){
			Tis_Window.Read_Ptraj_info(checkpoint_time);
			
		}	
		
	   
		// init moves
		
	    order_parameter Order_Parameter(Box, Particles, Fileio.MAX_coll_p_in);
		cluster Cluster(Box);
		
		move Move(Box->N, Particles.edge_N, Fileio.delta_tmax, Fileio.delta_rmax, Fileio.delta_Vmax, Fileio.is_translation_ON, Fileio.is_rotation_ON, Fileio.is_volumemove_ON);
		//Move.Calculate_Cluster_List(Particles, Box);
		Move.Set_Random_State(r,r01);
	
		int time;
		int MC_cycle_time;
		time = start_time;
	
		Box->P_sigma = Fileio.P_sigma_in;
		
		cout<<"Packing_fraction: "<<Box->packing_fraction<<endl;
		
	
		//Fileio.Save_Config(Box, Particles, Move.r, Move.r01, Move.dmax_t, Move.dmax_alpha, Move.dmax_beta, Move.dmax_gamma, Move.dmax_V, Move.dmax_L, time);
		  
		do{  
			
			time = time + 1;
			//cout<<"time:"<<time<<endl;
			//MC_cycle_time = 2*Box->N + 1;
			MC_cycle_time = Move.mt_sum;
			
		
			for(int cycle_time = 0; cycle_time<MC_cycle_time; cycle_time++){
				 	
			
				 	
				Move.Iterate(Particles, Box, Fileio, time);
				//TEST
				//cout<<"time "<<(time-1)*MC_cycle_time + cycle_time<<endl;
				//TEST
					
			}
			
		  
			if(time%calc_frequency==0){
			  
				cout<<"time:"<<time<<endl;
			  
			 
				Order_Parameter.Calculate_Local_Order_Parameters(Particles, Particles.N_Particle, Box);
				Order_Parameter.Caculate_is_nPhase(Box);
				Cluster.Calculate(Order_Parameter, Particles, Box);
				cout<<"Cluster.Size_R: "<<Cluster.Size_R<<endl;
				
				Fileio.Write_Order_Parameters(time, D2_order, D4_order);
				Fileio.Write_Local_Order_Parameters(time, Particles, Box, Order_Parameter, Cluster);
				//Fileio.Write_Cluster_List(Move.Cluster_List, Move.cluster_size, Particles, Order_Parameter);
				
				Move.Calculate_Acceptances(time);
			
				Fileio.Write_Acceptances(time, Move.accept_translate_procent, Move.accept_rotate_procent, Move.accept_iso_vol_procent, Move.accept_complete_procent);
				Fileio.Write_NPT(time, Box);
				
				
				if(strcmp(runtype, tis_runtype)==0){
			
					//cout<<"exit_code "<<exit_code<<endl;
					//cout<<"Cluster.Size_R: "<<Cluster.Size_R<<endl;
					current_window_value = Cluster.Size_R;
					exit_code = Tis_Window.tis_controller(current_window_value, time, N_check_points, checkpoint_frequency);
					//time = Tis_Window.time_N;
				}	
				
				if(strcmp(runtype, flux_runtype)==0){
			
					current_window_value = Cluster.Size_R;
					Tis_Window.Get_FLUX_path_info(current_window_value);
					
				}
					
				if(strcmp(runtype, tis_runtype)==0){	
						current_window_value = Cluster.Size_R;
						Tis_Window.Get_Ptraj_info(current_window_value, time);
						//Tis_Window.Ptraj_controller(current_window_value, time);
						
				}		
				
			
				
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
	
		
		
		if(strcmp(runtype, tis_runtype)==0){	
			Tis_Window.Write_Ptraj_info();
		}	
			
			
		if(strcmp(runtype, flux_runtype)==0){	
			Tis_Window.Write_FLUX_path_info(time, Move.dmax_t);
			
		}
		
	

		Fileio.Save_Config(Box, Particles, Move.r, Move.r01, Move.dmax_t, Move.dmax_alpha, Move.dmax_beta, Move.dmax_gamma, Move.dmax_V, Move.dmax_L, term_time);
		

		
		//Fileio.Write_g_r(term_time,Order_Parameter);
		//Order_Parameter.write_histogram_info();
		
		cout<<"End of Simulation. BYE BYE"<<endl;
		
		
	
	}  

     	 
     
		
