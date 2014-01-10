# include "pp_tis_window.h"



	tis_window::tis_window( int wv_0_in, int wv_i_in, int wv_ip1_in, int term_time, char traj_type_in[128]){
		
		//w_number = wn_in;
		//wv_current = wv_current_in;
		wv_0 = 	wv_0_in;
		wv_i = wv_i_in;
		wv_ip1 = wv_ip1_in;
		
	
		w_term_time = term_time;
		
		is_proper_path = 0;
		is_lambda_0_path = 0;
		is_lambda_p1_path = 0;
		exceded_max_path_length = 0;
		comming_from_basinA = 1;
		has_crossed_wv_i = 0;
		has_crossed_wv_im1 = 0;
		has_crossed_wv_0 = 0;
		time_tau_b = 0;
		time_tau_f = 0;
		no_crossing_with_lambda_i = 0;
		no_crossing_with_lambda_0 = 0;
		comming_from_inside=0;
		comming_from_outside=0;
		has_crossed_NOW = 0;
		crossing_time = -1;
		//is_lambda_i_trajectory = 0;
		traj_info = 0;
		wv_max = 0;
		
		exit_code = false;	
		
		update_state = 0;
		set_state = 0;
		
		N_bins = 4;
		delta_lambda = 5.0;
		tis_time = 1;
		//xmin = double(wv_i - 1);
		//xmax = double(wv_ip1 + 4);

		//histogram = gsl_histogram_alloc (N_bins);
		//set_state = gsl_histogram_set_ranges_uniform (histogram, xmin, xmax);    
			
		N_positive_crossing = 0;  
		
		random_counter = 0;
		string_random1 = "RANDOM_STATE_R_";
		string_random2 = "RANDOM_STATE_R1_";
		string_bin = ".bin";
		
		strcpy(traj_type,traj_type_in);
		
		/*
		for (int i=0; i<100;i++){
			random_list[i] = -1;
			
		}
		*/
		
			
	}



	void tis_window::Get_Ptraj_info(int wv, int mc_time){
		
		if((wv>=wv_i)&&(traj_info==0)){
			traj_info = 1;
			crossing_time = mc_time;
			//has_crossed_NOW = 1;
			//cout<<"TRAJ_INFO: 1"<<endl;
			
		}
		

	}
	
	
	void tis_window::Get_Ptraj_info_BACKWARDS(int wv, int mc_time){
		if((wv<=wv_i)&&(traj_info==0)){
			traj_info = 1;
			crossing_time = mc_time;
		
		}
	
	}
	

	void tis_window::Read_Ptraj_info(int checkpoint_time){
		
	
		tis_string = "traj.info";	
		
		ifstream traj_info_in(tis_string.c_str());
		traj_info_in>>crossing_time;
		traj_info_in.close();
        
        if(crossing_time>checkpoint_time){
			traj_info = 0;
        }
        
        if(crossing_time <= checkpoint_time){
			traj_info = 1;
		}
       
		has_crossed_wv_i = traj_info;
		
		
	
	}
	
	void tis_window::Write_Ptraj_info(){
		
	
		tis_string = "traj.info";	
		
		ofstream traj_info_out(tis_string.c_str());
		traj_info_out<<crossing_time<<endl;
		traj_info_out.close();
			
			
	}
	
	void tis_window::Ptraj_controller(int wv, int mc_time){
	
	
		if(strcmp(traj_type,"forwards")==0){
			Get_Ptraj_info(wv, mc_time);
		}	
		
		if(strcmp(traj_type,"backwards")==0){
			Get_Ptraj_info_BACKWARDS(wv, mc_time);
			
		}	
			
		
		
	}	
	


	void tis_window::Get_path_info(int wv, int mc_time){
		
		
		
		
		if((wv<wv_0) && (has_crossed_wv_i==1)){
			
			cout<<"Reached lambda 0 with crossing lambda i before"<<endl;
			is_proper_path = 1;
			is_lambda_0_path = 1;
			
		}	
		
		
		
		if((wv>=wv_i)&&(has_crossed_wv_i==0)){
		
			has_crossed_wv_i = 1;
			
			//time_tau_f = 0;
			
		}	
		
		if((wv<wv_0)&&(has_crossed_wv_i==0)){
		
			no_crossing_with_lambda_i = 1;
			//time_tau_f = 0;
			//tis_time = 0;
			
			
		}	
		
	
		
		if(wv >= wv_ip1){
			cout<<"Reached lambda i+1 with crossing lambda i before"<<endl;
			is_proper_path = 1;
			is_lambda_p1_path = 1;
			
		}	
		
		if(mc_time == w_term_time){
			cout<<"Exceded maximum path length"<<endl;
			exceded_max_path_length = 1;
			
		
		}
		
		
		
		if((wv>wv_0)&&(mc_time==1)){
		
		  is_proper_path = 0;
		  no_crossing_with_lambda_0 =1;
		
		}
		
	
		
	}


	
	void tis_window::Get_path_info_BACKWARDS(int wv, int mc_time){
		
	
		if((wv>wv_0) && (has_crossed_wv_i==1)){
			
			cout<<"Reached lambda 0 with crossing lambda i before"<<endl;
			is_proper_path = 1;
			is_lambda_0_path = 1;
			
		}	
		
		
		
		if((wv<=wv_i)&&(has_crossed_wv_i==0)){
		
			has_crossed_wv_i = 1;
			
			//time_tau_f = 0;
			
		}	
		
		if((wv>wv_0)&&(has_crossed_wv_i==0)){
		
			no_crossing_with_lambda_i = 1;
			//time_tau_f = 0;
			//tis_time = 0;
			
			
		}	
		
	
		
		if(wv <= wv_ip1){
			cout<<"Reached lambda i+1 with crossing lambda i before"<<endl;
			is_proper_path = 1;
			is_lambda_p1_path = 1;
			
		}	
		
		if(mc_time == w_term_time){
			cout<<"Exceded maximum path length"<<endl;
			exceded_max_path_length = 1;
			
		
		}
		
		
		
		if((wv<wv_0)&&(mc_time==1)){
		
		  is_proper_path = 0;
		  no_crossing_with_lambda_0 =1;
		
		}
			
		
	}


	void tis_window::Get_FLUX_path_info(int wv){
		
		cout<<"wv_i: "<<wv_i<<endl;
		cout<<"wv: "<<wv<<endl;
		
		
		
		if((wv > wv_i)&&(comming_from_inside == 1)){
			cout<<"Crossing from inside "<<endl;
			N_positive_crossing = N_positive_crossing +1;
			comming_from_inside=0;
			comming_from_outside=1;
			
			/*
			ofstream cross_out("crossing_points.dat", ios::out | ios::app);
			cross_out<<mc_time<<" "<<wv<<endl;
			cross_out.close();
			*/
			
		}
		
		if((wv < wv_i)&&(comming_from_outside == 1)){
			
			N_positive_crossing = N_positive_crossing +1;
			cout<<"Crossing from outside"<<endl;
			comming_from_outside=0;
			comming_from_inside=1;
			
		}		
		
	
		
		
	}

	void tis_window::Get_FLUX_path_info_BACKWARDS(int wv){
		
		cout<<"wv_i: "<<wv_i<<endl;
		cout<<"wv: "<<wv<<endl;
		
		
		
		if((wv < wv_i)&&(comming_from_inside == 1)){
			cout<<"Crossing from inside "<<endl;
			N_positive_crossing = N_positive_crossing +1;
			comming_from_inside=0;
			comming_from_outside=1;
			
			/*
			ofstream cross_out("crossing_points.dat", ios::out | ios::app);
			cross_out<<mc_time<<" "<<wv<<endl;
			cross_out.close();
			*/
			
		}
		
		if((wv > wv_i)&&(comming_from_outside == 1)){
			
			N_positive_crossing = N_positive_crossing +1;
			cout<<"Crossing from outside"<<endl;
			comming_from_outside=0;
			comming_from_inside=1;
			
		}		
		
	
		
		
	}


	void tis_window::Write_path_info(int N_check_points, int check_point_frequency, int mc_time){
		
		N_check_points = mc_time/check_point_frequency;
		
		ofstream path_info_out("path.info");
		
		path_info_out<<is_proper_path<<endl;
		path_info_out<<is_lambda_0_path<<endl;
		path_info_out<<is_lambda_p1_path<<endl;
		path_info_out<<exceded_max_path_length<<endl;
		path_info_out<<no_crossing_with_lambda_i<<endl;
		path_info_out<<mc_time<<endl;
		path_info_out<<N_check_points<<endl;
		path_info_out<<check_point_frequency<<endl;
		//path_info_out<<wv_max<<endl;
		
		
		path_info_out.close();
		
		
		
	}			
	
	void tis_window::Write_FLUX_path_info(int mc_time, int delta_t){
		
		
		ofstream path_info_out("flux_path.info");
		
		path_info_out<<N_positive_crossing<<endl;
		path_info_out<<mc_time<<endl;
		path_info_out<<delta_t<<endl;
		
		path_info_out.close();
		
		
	}			
	

	bool tis_window::tis_controller(int wv, int mc_time, int N_check_points, int check_point_frequency){
	
		
		if(strcmp(traj_type,"forwards")==0){
			Get_path_info(wv, mc_time);
		
		}
		
		if(strcmp(traj_type,"backwards")==0){
			Get_path_info_BACKWARDS(wv, mc_time);
		
		}	
		
		
		
		//if((is_proper_path ==1)||(exceded_max_path_length==1)||(no_crossing_with_lambda_i==1)||(no_crossing_with_lambda_0==1)){
		if((is_proper_path ==1)||(exceded_max_path_length==1)||(no_crossing_with_lambda_i==1)||(no_crossing_with_lambda_0==1)){	
			
			Write_path_info(N_check_points, check_point_frequency, mc_time);
			//write_histogram_info();
		 	exit_code = true;
		}

		return exit_code;

	}
		
	
	bool tis_window::FLUX_controller(int wv, int mc_time, int delta_t){	

		Get_FLUX_path_info(wv);
		
		if(exceded_max_path_length==1){
			Write_FLUX_path_info(mc_time, delta_t);
		 	exit_code = true;
		}
	
		return exit_code;

	}


	/*

	void tis_window::set_histogram(int N_bins, int xmin, int xmax){
	
	
		histogram = gsl_histogram_alloc(N_bins);
		
		set_state = gsl_histogram_set_ranges_uniform (histogram, xmin, xmax);      
		
	}

	void tis_window::unset_histogram(){
		
		gsl_histogram_free (histogram);			
			
	}


	void tis_window::update_histogram(int wv){
		
		//cout<<"hello histo"<<endl;
		r_size = double(wv);
		cout<<r_size<<endl;
		if((r_size >=xmin)&&(r_size<=xmax)){
			update_state = gsl_histogram_increment (histogram, r_size);
			cout<<"update_histogram: "<<update_state<<endl;
		}
	
	}	

	void tis_window::write_histogram_info(){
		
		FILE *fp_h;
		fp_h = fopen("histogram.info","w");
		save_state = gsl_histogram_fprintf (fp_h, histogram, "%f", "%f");	
		close_state = fclose(fp_h);
		
			
	}	


	*/





