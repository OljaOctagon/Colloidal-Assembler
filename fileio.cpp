

#include "fileio.h"


	void fileio::Write_Order_Parameters(int mc_time, double D2_order, double D4_order){
	
		string_3 = string_D2 + string_dat;
		ofstream D2_out(string_3.c_str(), ios::out | ios::app);
		
		string_3 = string_D4 + string_txt;
		ofstream D4_out(string_3.c_str(), ios::out | ios::app);
		
	
		D2_out<<mc_time<<"     "<<D2_order<<endl;
		
		D4_out<<mc_time<<"     "<<D4_order<<endl;
		
		D2_out.close();
		D4_out.close();
		
	
	
      }


	void fileio::Write_Local_Order_Parameters(int mc_time, particles& Particles, box* Box, order_parameter& Order_Parameter, cluster& Cluster){			
		
		/*
		 string_3 = string_local_OP + string_txt;
		 ofstream op_out(string_3.c_str(), ios::out | ios::app);
		
		 for(id=0;id<Box->N;id++){
			
			 //op_out<<mc_time<<"     "<<Order_Parameter.il2[id]<<"      "<<Order_Parameter.il4[id]<<"        "<<Order_Parameter.q4.ql[id]<<"      "<<Particles.N_Particle[id].x_center<<"      "<<Particles.N_Particle[id].y_center<<"      "<<Particles.N_Particle[id].z_center<<endl;
			 op_out<<mc_time<<"    "<<Order_Parameter.il2[id]<<"    "<<Order_Parameter.il4[id]<<"    "<<Order_Parameter.q4.ql[id]<<"    "<<Order_Parameter.N_obonds[id]<<"   "<<Particles.N_Particle[id].x_center<<"      "<<Particles.N_Particle[id].y_center<<"      "<<Particles.N_Particle[id].z_center<<endl;	
			
			}	
		
		  op_out.close();
	
		*/
		
		
			
		ofstream r_out("Largest_R.dat", ios::out | ios::app);
		r_out<<mc_time<<" "<<Cluster.Size_R<<endl;
		r_out.close();	
		
		
		
		/*		
		ofstream vmd_op1_out("op1.dat", ios::out | ios::app);
		
		vmd_op1_out<<Box->N<<endl;
		vmd_op1_out<<"Cube Particles"<<endl;
		
		
		for(id=0;id<Box->N;id++){
			//vmd_op1_out<<Order_Parameter.il2[id]<<"   "<<Order_Parameter.q4.ql[id]<<"    "<<Order_Parameter.N_obonds[id]<<endl;
			vmd_op1_out<<Order_Parameter.q4.ql[id]<<endl;
		}	
		
		vmd_op1_out.close();
		
		
		ofstream vmd_op2_out("op2.dat", ios::out | ios::app);
		
		vmd_op2_out<<Box->N<<endl;
		vmd_op2_out<<"Cube Particles"<<endl;
		
		
		for(id=0;id<Box->N;id++){
			//vmd_op1_out<<Order_Parameter.il2[id]<<"   "<<Order_Parameter.q4.ql[id]<<"    "<<Order_Parameter.N_obonds[id]<<endl;
			vmd_op2_out<<Order_Parameter.N_obonds[id]<<endl;
		}	
			
		vmd_op2_out.close();
		
		
		ofstream vmd_op3_out("op3.dat", ios::out | ios::app);
		
		vmd_op3_out<<Box->N<<endl;
		vmd_op3_out<<"Cube Particles"<<endl;
		
		
		for(id=0;id<Box->N;id++){
		
			vmd_op3_out<<Order_Parameter.cubic_op3_av[id]<<endl;
		}	
			
		vmd_op3_out.close();
		
		
		ofstream vmd_op4_out("op4.dat", ios::out | ios::app);
		
		vmd_op4_out<<Box->N<<endl;
		vmd_op4_out<<"Cube Particles"<<endl;
		
		
		for(id=0;id<Box->N;id++){
			
			vmd_op4_out<<Order_Parameter.q4.N_qbonds[id]<<endl;
			
		}	
		
		vmd_op4_out.close();
		*/
		
		
	
		
	}	
		

	
	void fileio::Write_g_r(int mc_time, order_parameter& Order_Parameter){
		
		string_3 = string_g_r + string_dat;
		ofstream gr_out(string_3.c_str(), ios::out | ios::app);
		
			for(int p = 0;p<Order_Parameter.N_histo_points;p++){
		
				gr_out<<mc_time<<"    "<<Order_Parameter.g_radial_distribution[p].x<<"   "<<Order_Parameter.g_radial_distribution[p].y<<endl;
			}
			
		gr_out.close();	
	
		
	}



	void fileio::Write_Cluster_List(int* Cluster_List, int cluster_size, particles& Particles, order_parameter& Order_Parameter){
		
		 string_3 = string_cluster_list + string_dat;
		
		 ofstream cluster_out(string_3.c_str(), ios::out | ios::app);
			
		 int cluster_id;	
		
	     for(id=0;id<cluster_size; id++){
			  
			  cluster_id = Cluster_List[id];
			  cluster_out<<cluster_id<<"    "<<Order_Parameter.q4.ql[cluster_id]<<"      "<<Particles.N_Particle[cluster_id].x_center<<"      "<<Particles.N_Particle[cluster_id].y_center<<"      "<<Particles.N_Particle[cluster_id].z_center<<endl;
			
		 }	
		
		 cluster_out.close();

	    }


	void fileio::Write_Vol_Packing(int mc_time, box* Box){
		
		 string_3 = string_Vol_Packing + string_dat;
		 ofstream vol_out(string_3.c_str(), ios::out | ios::app);
			
		 vol_out<<mc_time<<"     "<<Box->V<<"     "<<Box->packing_fraction<<endl;
			
		 vol_out.close();
			
		}


	void fileio::Save_Config(box* Box, particles& Particles, gsl_rng *r, gsl_rng *r01, double dmax_t, double dmax_alpha, double dmax_beta, double dmax_gamma, double dmax_V, double dmax_L, int mc_time){

		
		 //save particle positions
		 Save_Positions(Box, Particles, mc_time);

		 //save particle orientation
		 Save_Orientations(Box, Particles, mc_time);
		//Save_Orientations_EULER(Box, Particles.N_Particle);
		 //Save_Orientations_Axis_Angle(Box, Particles, mc_time);

		 //save Box coordinates in binary format
		 Save_Box(Box, mc_time);

		 //save state of random number generator
		 Save_Random_State(r, r01, mc_time);	

	    }


	void fileio::Read_Config(box* Box, particles& Particles, gsl_rng *r, gsl_rng *r01, int mc_time){
		
		 //read Box config file	
		 Read_Box(Box, mc_time);	

		 //read positions config file
		 Read_Positions(Box, Particles, mc_time);

		 // read orientation config file
		 Read_Orientations(Box, Particles, mc_time);
		 //Read_Orientations_EULER(Box, Particles);
		 //Read_Orientations_Axis_Angle(Box, Particles, mc_time);
		 
		 // read former state of random number generator
		 Read_Random_State(r,r01, mc_time);
			
	    }


	void fileio::Read_Parameters(){
	
		cout<<"start reading parameters"<<endl;
	
		boost::property_tree::ptree pt;
		boost::property_tree::ini_parser::read_ini("para.ini", pt);
		
	
		

		N_in 				= pt.get<int>("System.Number_of_Particles");
		P_sigma_in  		= pt.get<double>("System.Pressure");
		packing_fraction_in = pt.get<double>("System.Packing_Fraction");
		number_of_cells_in  = pt.get<int>("System.Number_of_Cells");
		MAX_coll_p_in 		= pt.get<int>("Cubes.Maximum_Collision_Partners");
		MAX_fshell_p_in 	= pt.get<int>("Cubes.Maximum_First_Shell_Partners");
		
	
		
	
		calc_frequency_in       = pt.get<int>("Output.Calculation_Frequency");
		frame_frequency_in  	= pt.get<int>("Output.Frame_Frequency");
		checkpoint_frequency_in = pt.get<int>("Output.Checkpoint_Frequency"); 
		
		is_translation_ON = pt.get<string>("Monte_Carlo_parameters.Translation");
		is_rotation_ON    = pt.get<string>("Monte_Carlo_parameters.Rotation");
		is_volumemove_ON  = pt.get<string>("Monte_Carlo_parameters.Volume_Move");
		
		
	
		delta_tmax = pt.get<double>("Monte_Carlo_parameters.Sigma_Translation");
		delta_rmax = pt.get<double>("Monte_Carlo_parameters.Sigma_Rotation");
		delta_Vmax = pt.get<double>("Monte_Carlo_parameters.Sigma_Volume");
		
		seed1_in = pt.get<int>("Monte_Carlo_parameters.rng_seed1_in");
		seed2_in = pt.get<int>("Monte_Carlo_parameters.rng_seed2_in");
		
		
		//Initiialize Random State
			//gsl_rng_set(r01,Fileio.seed1_in);
			//gsl_rng_set(r, Fileio.seed2_in);
		
	
		cout<<"Finished reading Parameters"<<endl;
			
	
	   
    }
  

	void fileio::Save_Positions(box* Box, particles& Particles, int mc_time){
	
		//output in float/double format to a txt file
		/*
  	    string_3 = string_positions + string_txt;
		
		ofstream pos_out(string_3.c_str());
		
		for(int k=0;k<Box->N;k++){
			
			 pos_out<<Particles.N_Particle[k].x_center<<"     "<<Particles.N_Particle[k].y_center<<"     "<<Particles.N_Particle[k].z_center<<endl;
		
		    }
		    
		pos_out.close();
	
		*/
	
	
		stringstream mout;
		mout << mc_time;
		string_m = mout.str();	
		
		
		string_3 = string_positions + string_m + string_bin;	
			
		ofstream pos_out1(string_3.c_str(), ios::out | ios::binary);
		
		for(int id = 0; id<Box->N; id++){
			pos_out1.write(reinterpret_cast<char*>(&Particles.N_Particle[id].x_center), sizeof(Particles.N_Particle[id].x_center));
			pos_out1.write(reinterpret_cast<char*>(&Particles.N_Particle[id].y_center), sizeof(Particles.N_Particle[id].y_center));
			pos_out1.write(reinterpret_cast<char*>(&Particles.N_Particle[id].z_center), sizeof(Particles.N_Particle[id].z_center));
			
		}
		
		pos_out1.close();
			
    }
    

	void fileio::Read_Positions(box* Box, particles& Particles, int mc_time){
		/*
		 ifstream pos_in;
		
		 pos_in.open("positions.txt");
		
		 if ( ! pos_in ) {
			 cout << "Error: Can't open the file. Please use the filename 'positions.txt' for position input \n";
			 exit(1);
		    }

		 for(int id=0; id<Box->N;id++){
			
			 pos_in>>Particles.N_Particle[id].x_center;
			 pos_in>>Particles.N_Particle[id].y_center;
			 pos_in>>Particles.N_Particle[id].z_center;
			 
			}

		 pos_in.close();
		*/
		
		stringstream mout;
		mout << mc_time;
		string_m = mout.str();	
		
		
		
	
		string_3 = string_positions + string_m + string_bin;	
			
		ifstream pos_in1(string_3.c_str(), ios::in | ios::binary);
		
		for(int id = 0; id<Box->N; id++){
			pos_in1.read(reinterpret_cast<char*>(&Particles.N_Particle[id].x_center), sizeof(double));
			pos_in1.read(reinterpret_cast<char*>(&Particles.N_Particle[id].y_center), sizeof(double));
			pos_in1.read(reinterpret_cast<char*>(&Particles.N_Particle[id].z_center), sizeof(double));
			
		}
		
		pos_in1.close();
			
    } 


	void fileio::Save_Orientations(box* Box, particles& Particles, int mc_time){
	
		
		stringstream mout;
		mout << mc_time;
		string_m = mout.str();	
		
		
		string_3 = string_orientations + string_m + string_bin;
		
		ofstream orient_out(string_3.c_str(), ios::out | ios::binary);
		
		for(int id=0;id<Box->N;id++){
		
				
			orient_out.write(reinterpret_cast<char*>(&Particles.N_Particle[id].q.x), sizeof(Particles.N_Particle[id].q.x));
			orient_out.write(reinterpret_cast<char*>(&Particles.N_Particle[id].q.y), sizeof(Particles.N_Particle[id].q.y));
			orient_out.write(reinterpret_cast<char*>(&Particles.N_Particle[id].q.z), sizeof(Particles.N_Particle[id].q.z));
			orient_out.write(reinterpret_cast<char*>(&Particles.N_Particle[id].q.w), sizeof(Particles.N_Particle[id].q.w));
			
		}
		
				
		orient_out.close();
		
			
	}
	
	

	void fileio::Read_Orientations(box* Box, particles& Particles, int mc_time){
	
		stringstream mout;
		mout << mc_time;
		string_m = mout.str();	
	
	
		string_3 = string_orientations + string_m + string_bin;	
			
		ifstream orient_in(string_3.c_str(), ios::in | ios::binary);
		
		for(int id = 0; id<Box->N; id++){
			orient_in.read(reinterpret_cast<char*>(&Particles.N_Particle[id].q.x), sizeof(double));
			orient_in.read(reinterpret_cast<char*>(&Particles.N_Particle[id].q.y), sizeof(double));
			orient_in.read(reinterpret_cast<char*>(&Particles.N_Particle[id].q.z), sizeof(double));
			orient_in.read(reinterpret_cast<char*>(&Particles.N_Particle[id].q.w), sizeof(double));
			
		}
		
		orient_in.close();
			


	}



	void fileio::Save_Box(box* Box, int mc_time){
	
		
					  
		stringstream mout;
		mout << mc_time;
		string_m = mout.str();	
		
		
		string_3 = string_Box + string_m + string_bin;	
			
		ofstream box_out(string_3.c_str(), ios::out | ios::binary);
		
		box_out.write(reinterpret_cast<char*>(&Box->x_center), sizeof(double));
		box_out.write(reinterpret_cast<char*>(&Box->y_center), sizeof(double));
		box_out.write(reinterpret_cast<char*>(&Box->z_center), sizeof(double));
		
		box_out.write(reinterpret_cast<char*>(&Box->Lx), sizeof(double));
		box_out.write(reinterpret_cast<char*>(&Box->Ly), sizeof(double));
		box_out.write(reinterpret_cast<char*>(&Box->Lz), sizeof(double));
		
		box_out.write(reinterpret_cast<char*>(&Box->V), sizeof(double));
		box_out.write(reinterpret_cast<char*>(&Box->packing_fraction), sizeof(double));
		
		box_out.close();
	  
	
	}	


	void fileio::Read_Box(box* Box, int mc_time){


		
		stringstream mout;
		mout << mc_time;
		string_m = mout.str();	
		
		
		
		string_3 = string_Box + string_m + string_bin;	
			
		ifstream box_in(string_3.c_str(), ios::in | ios::binary);
		
		box_in.read(reinterpret_cast<char*>(&Box->x_center), sizeof(double));
		box_in.read(reinterpret_cast<char*>(&Box->y_center), sizeof(double));
		box_in.read(reinterpret_cast<char*>(&Box->z_center), sizeof(double));
		
		box_in.read(reinterpret_cast<char*>(&Box->Lx), sizeof(double));
		box_in.read(reinterpret_cast<char*>(&Box->Ly), sizeof(double));
		box_in.read(reinterpret_cast<char*>(&Box->Lz), sizeof(double));
		
		box_in.read(reinterpret_cast<char*>(&Box->V), sizeof(double));
		box_in.read(reinterpret_cast<char*>(&Box->packing_fraction), sizeof(double));
		
		box_in.close();
		
	
        }


	void fileio::Save_Random_State(gsl_rng *r, gsl_rng *r01, int mc_time){
		
		  
		stringstream mout;
		mout << mc_time;
		string_m = mout.str();		
	
		string_3 = string_random1 + string_m + string_bin;	
		

		FILE *fp_r;
		fp_r = fopen(string_3.c_str(), "w");
		save_state = gsl_rng_fwrite(fp_r, r);
		cout<<"save_state "<<save_state<<endl;
		close_state = fclose(fp_r);

		string_3 = string_random2 + string_m + string_bin;

		FILE *fp_r01;
		fp_r01 = fopen(string_3.c_str(), "w");
		save_state = gsl_rng_fwrite(fp_r01, r01);
		cout<<"save_state"<<save_state<<endl;
		close_state = fclose(fp_r01);
		


	}


	void fileio::Read_Random_State(gsl_rng *r, gsl_rng *r01, int mc_time){


		stringstream mout;
		mout << mc_time;
		string_m = mout.str();		
	
		string_3 = string_random1 + string_m + string_bin;	
		
		FILE *fp_in_r;
		fp_in_r = fopen(string_3.c_str(), "r");
		read_state=gsl_rng_fread (fp_in_r, r);	
		cout<<"read_state random1"<<endl;
		close_state = fclose(fp_in_r);

		
		string_3 = string_random2 + string_m + string_bin;	

		FILE *fp_in_r01;
		fp_in_r01 = fopen(string_3.c_str(), "r");
		read_state=gsl_rng_fread (fp_in_r01, r01);
		close_state = fclose(fp_in_r01);
		
	}


	void fileio::Set_Vars(){
  
		string_1 = "positions_"; 
		string_positions = "positions_";
		string_g = "g_positions_";
		string_2 = ".xyz";
		string_4 = "_";
		string_Cut_Off= "Cut_Off_";
		string_V_List = "V_List_";
		string_txt = ".txt";
		string_Cell = "Cell_";
		string_C_List = "C_List_";
		string_Box = "Box_";
		string_para_out = "Para_Out";
		string_para_in = "Para_In";
		string_orientations = "orientations_";
		string_D4 = "D4_order";
		string_D2 = "D2_order";	
		string_Vol_Packing = "Vol_Packing";
		string_acceptances = "Acceptances";
		string_NPT = "NPT_OUT";
		string_i2 = "i2_local_order";
		string_cluster_list = "Cluster_List";
		string_local_OP = "local_OP";
		string_g_r = "g_r";
		string_dat = ".dat";
		string_bin = ".bin";
		string_random1 = "RANDOM_STATE_R_";
		string_random2 = "RANDOM_STATE_R1_";
    
	   } 


	void fileio::Write_Positions(int mc_time, box* Box, particles& Particles){
	
		string_3 = string_positions + string_2;
		
		ofstream xyz_out(string_3.c_str(), ios::out | ios::app);
		
		M=Box->N*Particles.N_Particle[0].edge_N;		
		xyz_out<<M<<endl;
		xyz_out<<"Particles of frame "<<mc_time<<endl;
		
		
		for(int k=0;k<Box->N;k++){
		
				for(int j=0;j<Particles.N_Particle[k].edge_N;j++){
					
					 xyz_out<<"cb"<<"       "<<Particles.N_Particle[k].x[j]<<"   "<<Particles.N_Particle[k].y[j]<<"   "<<Particles.N_Particle[k].z[j]<<endl;
					 
					}   
			}
			
		xyz_out.close();
       
       	
		 
		ofstream center_out("pos_c.xyz", ios::out | ios::app);
	  
		center_out<<Box->N<<endl;
		center_out<<"Centers of frame "<<mc_time<<endl;
	  
		for(int k=0;k<Box->N;k++){
			center_out<<"cb"<<"       "<<Particles.N_Particle[k].x_center<<"   "<<Particles.N_Particle[k].y_center<<"   "<<Particles.N_Particle[k].z_center<<endl;
		}				
	
		center_out.close();
		
		
	}

	void fileio::Write_Positions_SPHERE(int mc_time, box* Box, particles& Particles){
		
			string_3 = string_positions + string_2;
			
			ofstream xyz_out(string_3.c_str(), ios::out | ios::app);
			
			M=Box->N*Particles.N_Particle[0].edge_N;
			
			xyz_out<<M<<endl;
			xyz_out<<"Particles of frame "<<mc_time<<endl;
			
			for(int k=0;k<Box->N;k++){
				
				 xyz_out<<"sp"<<"       "<<Particles.N_Particle[k].x_center<<"   "<<Particles.N_Particle[k].y_center<<"   "<<Particles.N_Particle[k].z_center<<endl;
				  
				}

			xyz_out.close();
		   
		   }	
		
	
	void fileio::Write_Box_Positions( int mc_time, box* Box){
		
		string_3 = string_Box + string_2;	
		
		ofstream bxyz_out(string_3.c_str(), ios::out | ios::app);
			
		M=9;
			
		bxyz_out<<M<<endl;
		bxyz_out<<"Box of frame "<<mc_time<<endl;
			
		bxyz_out<<"cb"<<"       "<<Box->x_center<<"   "<<Box->y_center<<"   "<<Box->z_center<<endl;
		
		for(int j=0;j<8;j++){
			
			bxyz_out<<"cb"<<"       "<<Box->x[j]<<"   "<<Box->y[j]<<"   "<<Box->z[j]<<endl;
		 
		 }   
		   
	
		bxyz_out.close();
			
       }
	

	void fileio::Write_Cell_Lists(int id, int mc_time, particles& Particles){
		
		stringstream mout;
		mout << mc_time;
		string_m = mout.str();
			 
		string_v_list = string_V_List + string_m + string_txt;	
		string_c_list = string_C_List + string_m + string_txt;	
			
		ofstream vlist_out(string_v_list.c_str());
			
		vlist_out<<"number of potential collision partners: "<<Particles.Collision_List[id].Nm<<endl;	
		
		
		for(int m = 0;m<Particles.Collision_List[id].Nm;m++){	
			vlist_out<<"collisoin partner"<<m<<" of "<<id<<" : "<<Particles.Collision_List[id].Elements[m].nl_id<<endl;
		}	
		
		vlist_out.close();


		ofstream clist_out(string_c_list.c_str());
		
		clist_out<<"Cell of id: "<<id<<": "<<Particles.Id_Cell_List[id]<<endl;
		
		clist_out.close();
		
	   }


	void fileio::Write_Cell_Positions(int number_of_cells, cell* Cell, int mc_time, box* Box){
		
		stringstream mout;
		mout << mc_time;
		string_m = mout.str();
		 
		
		string_cell = string_Cell + string_m + string_2;
 
		cout<<"Cell_edge"<<Cell[0].edge_N<<endl;
		M = Cell[0].edge_N +1;	
		
		ofstream cell_xyz_out(string_cell.c_str());
		
		
		for(int k=0;k<number_of_cells;k++){
		
			cell_xyz_out<<M<<endl;
			cell_xyz_out<<"Cell"<<k<<endl;
			cell_xyz_out<<k<<"     "<<Cell[k].x_center<<"     "<<Cell[k].y_center<<"     "<<Cell[k].z_center<<endl;
			
			for(int j=0;j<Cell[k].edge_N;j++){
				 cell_xyz_out<<k<<"       "<<Cell[k].x[j]<<"   "<<Cell[k].y[j]<<"   "<<Cell[k].z[j]<<endl;
				}   
	   
		
		   }
		
		 cell_xyz_out<<M<<endl;
		 cell_xyz_out<<"Box"<<endl;
		 cell_xyz_out<<number_of_cells<<"     "<<Box->x_center<<"     "<<Box->y_center<<"     "<<Box->z_center<<endl;
			
		 for(int j=0;j<Box->edge_N;j++){
			 cell_xyz_out<<number_of_cells<<"       "<<Box->x[j]<<"   "<<Box->y[j]<<"   "<<Box->z[j]<<endl;
			}   
	   
			
		 cell_xyz_out.close();

	    }	
		
					
	void fileio::Write_Acceptances(int mc_time, double accept_translate_procent, double accept_rotate_procent, double accept_iso_vol_procent, double accept_complete_procent){
		 
		 stringstream mout;
		 mout << mc_time;
		 string_m = mout.str();
		
		 string_3 = string_acceptances + string_txt;
			
		 ofstream accept_out(string_3.c_str(), ios::out | ios::app);
			
         accept_out<<mc_time<<"     "<<accept_complete_procent<<"     "<<accept_translate_procent<<"    "<<accept_rotate_procent<<"     "<<accept_iso_vol_procent<<endl;
					
		 accept_out.close();
		
			
		}	
			
				
	void fileio::Write_NPT(int mc_time, box* Box){
		
		 stringstream mout;
		 mout << mc_time;
		 string_m = mout.str();
		
		 string_3 = string_NPT + string_txt;
			
		 ofstream npt_out(string_3.c_str(), ios::out | ios::app);
			
	     npt_out<<mc_time<<"     "<<Box->P_sigma<<"     "<<Box->P<<"      "<<Box->V<<"      "<<Box->packing_fraction<<endl;
					
		 npt_out.close();
		
		}			
			
			
		
		
		
		
