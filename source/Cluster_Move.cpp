#include "move.h"
#include <gsl/gsl_math.h>

double move::Calculate_Potential(int a, int b, particles& Particles, box* Box){

	double patch_distance_squared;
	double epsilon_ij;
	int ela;
	int elb;
	double patch_energy_ab;

	epsilon_ij = 0;
	patch_energy_ab = 0;

	for (int pid1=0; pid1<Particles.N_Particle[a]->N_patches;pid1++){

		for (int pid2=0; pid2<Particles.N_Particle[b]->N_patches;pid2++){

			patch_distance_squared = Calculate_Patch_Distance(a,b , pid1, pid2, Particles, Box);
				
			//ONLY if all the radii are the same width

			if ( patch_distance_squared < Particles.N_Particle[a]->patch_cutoff_squared[pid1] && a!= b){
				
				ela = Particles.N_Particle[a]->patch_type[pid1];
				elb = Particles.N_Particle[b]->patch_type[pid2];

				patch_energy_ab = Particles.N_Particle[a]->patch_energy[ela][elb];
						
				epsilon_ij = epsilon_ij + patch_energy_ab;
					
			}	

		}	
	}

	
	return epsilon_ij;

}


void move::Reset_Pseudo_Cluster(box* Box){

	for(int i=0; i<Box->N;i++){
		is_element[i]=false;
		for (int j=0;j<Box->N;j++){		
		pseudo_cluster_info[i][j] = -1;
		}
		List[i] = -2;

	}
}


void move::Pseudocluster_Recursion(int id_j, int cn, int fl, particles& Particles, box* Box){

	int j;

	if (is_element[id_j] == false){
		List[N_List] = id_j;
		N_List = N_List+1;
		is_element[id_j]= true;

	}


	Particles.Collision_List[id_j].Calculate_OP(Box, id_j, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);

	
	for(int cp=0;cp<Particles.Collision_List[id_j].Nm;cp++){
		
		j = Particles.Collision_List[id_j].Elements[cp].nl_id;

		
		if (pseudo_cluster_info[id_j][j]==-1) {

			e_pair = Calculate_Potential(id_j,j, Particles, Box);
			 
			pseudo_cluster_info[id_j][j] = 1;
			pseudo_cluster_info[j][id_j] = 1; 
			
			b_factor = 1 - exp(beta_f*(e_pair));
			b_factor = GSL_MAX(0,b_factor);
			////cout<<"b_factor"<<b_factor<<endl;
			XI = gsl_rng_uniform(r01);

			if(XI<b_factor){

				cn = N_Bonds;
				N_Bonds = N_Bonds + 1; 
				
				Pseudocluster_Recursion(j, cn, fl, Particles, Box);
			}	
			

		}		
	}	

}


void move::Rot_Cluster_Move(particles& Particles, box* Box, fileio& Fileio, int mc_time) {
 	

	double Total_Energy_old;
	double delta_U;

	

 	N_List = 0;
 	N_Bonds= 0;
 	r_id = -1;
 	double phi_t;
 	double phi_r;

	//choose particle id

	id =  gsl_rng_uniform_int(r,Box->N);
	

	Reset_Pseudo_Cluster(Box);

	// calculate translation 
		
		
	//Rot_Update_Quarternions_VON_MISES(Particles, id);	
	double rand_phi;
	
	//rand_phi = gsl_ran_gaussian(r01, kappa);
	rand_phi = gsl_rng_uniform(r01);

	phi_t = kappa*(rand_phi - 0.5);
	
	Rot_mat[0] = cos(phi_t);
	Rot_mat[1] = -sin(phi_t);
	Rot_mat[2] = 0;
		
	Rot_mat[3] = sin(phi_t);
	Rot_mat[4] = cos(phi_t);
	Rot_mat[5] = 0;
		
	Rot_mat[6] = 0;
	Rot_mat[7] = 0;
	Rot_mat[8] = 1;
		
	int links, failed_links;
	links = 0;
	failed_links = 0;

	
	Pseudocluster_Recursion(id, links, failed_links, Particles, Box);

	// Collision Test
	
	
	double theta_x, theta_y, theta_av_x, theta_av_y;
	double xp1, xp2, yp1, yp2;

	xp1=0;
	xp2=0;
	yp1=0;
	yp2=0;

	for (int k=0; k<N_List;k++){
		int j;
 	 	 j= List[k];
 	 	 theta_x = (Particles.N_Particle[j]->x_center/Box->Lx)*2*M_PI;
 	 	 theta_y = (Particles.N_Particle[j]->y_center/Box->Ly)*2*M_PI;

 	 	 xp1 = xp1 + (Box->Lx/(2*M_PI))*cos(theta_x);
 	 	 xp2 = xp2 + (Box->Lx/(2*M_PI))*sin(theta_x);

 	 	 yp1 = yp1 + (Box->Ly/(2*M_PI))*cos(theta_y);
 	 	 yp2 = yp2 + (Box->Ly/(2*M_PI))*sin(theta_y);
	}

	xp1=xp1/double(N_List);
	xp2=xp2/double(N_List);
	yp1=yp1/double(N_List);
	yp2=yp2/double(N_List);
	N_List=int(N_List);

	theta_av_x =atan2(-xp2, -xp1) + M_PI;
	theta_av_y =atan2(-yp2, -yp1) + M_PI;

	center_mass.x = (Box->Lx/(2*M_PI))*(theta_av_x);
	center_mass.y = (Box->Ly/(2*M_PI))*(theta_av_y);
	center_mass.z = Box->Lz/2.0;


	//Check_Periodic_Center_of_Mass(center_mass, Box); 	 	 

	//cout<<"N_List"<<N_List<<endl;
	//cout<<"N_Bonds"<<N_Bonds<<endl;

	for (int k=0;k<N_List;k++){
 	 	int j;
 	 	j= List[k];
		Rot_Move_Map(Particles, j, Box, center_mass, Rot_mat);
		Particles.Check_Periodic_CM(j, Box);		
		Update_Periodic_Positions(Particles, Box, j);	
		Particles.N_Particle[j]->Calculate_Axis();
		Particles.N_Particle[j]->Calculate_Patch_Position();	
		
	}
	
	/*
	int M;	
	ofstream xyz_out("pos_trial_move.xyz", ios::out | ios::app);
	M=Particles.Res_Size*(Particles.N_Particle[0]->edge_N);	
	xyz_out<<M<<endl;
	xyz_out<<"Particles of frame x"<<endl;
		
		
	for(int k=0;k<Particles.Res_Size;k++){
		for(int j=0;j<Particles.N_Particle[k]->edge_N;j++){
	 		xyz_out<<"cb"<<"       "<<Particles.N_Particle[k]->x[j]<<"   "<<Particles.N_Particle[k]->y[j]<<"   "<<Particles.N_Particle[k]->z[j]<<endl;
					 
					}   
			}
			
	xyz_out.close();
    */

	int k=0;
	////cout<<"N_List   "<<N_List<<endl;
	////cout<<"N_Bonds  "<<N_Bonds<<endl;

	exit_status = 0;
	col_count = 0;

	int Max_Length;
	Max_Length = rint(1./gsl_rng_uniform(r01));
	//cout<<"Max_Length"<<Max_Length<<endl;
	//cout<<"N_List"<<N_List<<endl;

	if (N_List > Max_Length){
			exit_status = 1;
	}

	
	do {
		 Particles.Collision_List[List[k]].Calculate_OP(Box, List[k], Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);		
		 //Collision Test 
		 Collision_Test( Particles, Box, List[k], Particles.Collision_List);
		 k++;

	}while((exit_status==0)&&(k<N_List));

	//cout<<"exit status "<<exit_status<<endl;

	if (exit_status > 0){

		//Reset all positions!
		for (int k=0;k<N_List;k++){
			int j;
			j= List[k];
			Reset_Positions(Particles, j);
			Particles.N_Particle[j]->Calculate_Axis();
			Particles.N_Particle[j]->Calculate_Patch_Position();
			//Set_Positions(Particles, j);	
	
		}	
	}

	
 	if (exit_status == 0){

 		//cout<<"N_List "<<N_List<<" accepted"<<endl;
 		////cout<<"exit_status zero"<<endl;

 		Total_Energy_old = Particles.Total_Energy;
 		Calculate_Pair_Potential( Particles, Box);		
 		delta_U = Total_Energy - Total_Energy_old;
 		//cout<<"delta_U"<<delta_U<<endl;

 		
 		b_factor = GSL_MIN(1, exp((beta_f-beta)*delta_U));
 		XI = gsl_rng_uniform(r01);
 	
 		//b_factor=b_factor/double(N_List);
 		

 		// Cluster Move Acceptance criterium
 		//b_factor=1;

 		if (b_factor<XI){

 			for (int k=0;k<N_List;k++){
 				int j;
 				j= List[k];
 				Reset_Positions(Particles, j);
				Particles.N_Particle[j]->Calculate_Axis();
				Particles.N_Particle[j]->Calculate_Patch_Position();

 			}	
 			
 			Particles.Total_Energy = Total_Energy_old;
 			Total_Energy = Total_Energy_old;


 		}

 		////cout<<"b_factor energy "<<b_factor<<endl; 
 	

 		if (b_factor>=XI){

 			//cout<<"accept"<<endl;
 			//cout<<"N_List "<<N_List<<" accepted"<<endl;

 			for (int k=0;k<N_List;k++){
 				int j;
 				j = List[k];


 				Particles.N_Particle[j]->phi = Particles.N_Particle[j]->phi + phi_t;

 			
			 if (Particles.N_Particle[j]->phi > 2*M_PI){
				Particles.N_Particle[j]->phi = Particles.N_Particle[j]->phi - 2.0*M_PI;
				}
				
			 if (Particles.N_Particle[j]->phi < 0){	
				Particles.N_Particle[j]->phi = Particles.N_Particle[j]->phi + 2.0*M_PI;
				}

			

 				Set_Positions(Particles, j);	
 			}


 			Particles.Total_Energy = Total_Energy;


 		}
			

 	}

 	Reset_Pseudo_Cluster(Box);
 	
	
}


void move::Trans_Pseudocluster_Recursion(int id_j, int cn, int fl, particles& Particles, box* Box){

	
	int j;

	if (is_element[id_j] == false){
		List[N_List] = id_j;
		N_List = N_List+1;
		is_element[id_j]= true;

	}

	Particles.Collision_List[id_j].Calculate_OP(Box, id_j, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);


		for(int cp=0;cp<Particles.Collision_List[id_j].Nm;cp++){ 
			Particles.Collision_List_old[id_j].Elements[cp] = Particles.Collision_List[id_j].Elements[cp];
		}

	Particles.Collision_List_old[id_j].Nm = Particles.Collision_List[id_j].Nm;


	Trans_Update_Positions(Particles, id_j, trans_vec);
	Particles.Check_Periodic_CM(id_j, Box);
	Update_Periodic_Positions(Particles, Box, id_j);
	Particles.N_Particle[id_j]->Calculate_Axis();
	Particles.N_Particle[id_j]->Calculate_Patch_Position();	

	Particles.Collision_List[id_j].Calculate_OP(Box, id_j, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);

	Reset_Positions(Particles, id_j);
	Particles.N_Particle[id_j]->Calculate_Axis();
	Particles.N_Particle[id_j]->Calculate_Patch_Position();

  

	for(int cp=0;cp<Particles.Collision_List[id_j].Nm;cp++){
		
		Trans_Update_Positions(Particles, id_j, trans_vec);
		Particles.Check_Periodic_CM(id_j, Box);
		Update_Periodic_Positions(Particles, Box, id_j);
		Particles.N_Particle[id_j]->Calculate_Axis();
		Particles.N_Particle[id_j]->Calculate_Patch_Position();	

		j = Particles.Collision_List[id_j].Elements[cp].nl_id;

		
		if (pseudo_cluster_info[id_j][j]==-1) {

			e_single = Calculate_Potential(id_j,j, Particles, Box);
			 
			pseudo_cluster_info[id_j][j] = 1;
			pseudo_cluster_info[j][id_j] = 1; 
			//if (is_interacting=true){
				
			Particles.Collision_List[j].Calculate_OP(Box, j, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);
			for(int cp=0;cp<Particles.Collision_List[j].Nm;cp++){ 
				Particles.Collision_List_old[j].Elements[cp] = Particles.Collision_List[j].Elements[cp];
			}

			Particles.Collision_List_old[j].Nm = Particles.Collision_List[j].Nm;
				
			
			Trans_Update_Positions(Particles, j, trans_vec);
			Particles.Check_Periodic_CM(j, Box);
			Update_Periodic_Positions(Particles, Box, j);	
			Particles.N_Particle[j]->Calculate_Axis();
			Particles.N_Particle[j]->Calculate_Patch_Position();	
				
			e_pair = Calculate_Potential(id_j,j, Particles, Box);
				
			q_factor = exp(beta*(e_pair - e_single));	
			b_factor = 1 - exp(beta*(e_pair - e_single));
			
			b_factor = GSL_MAX(0, b_factor);
			q_factor = GSL_MIN(1, q_factor);
			
	

			XI = gsl_rng_uniform(r01);

			Reset_Positions(Particles, id_j);
			Particles.N_Particle[id_j]->Calculate_Axis();
			Particles.N_Particle[id_j]->Calculate_Patch_Position();	
				


			if (XI>=b_factor){

				Reset_Positions(Particles, j);
				Particles.N_Particle[j]->Calculate_Axis();
				Particles.N_Particle[j]->Calculate_Patch_Position();

				if (is_element[j] == true){
					
					fl = N_failed_links;

					q_f[fl] = q_factor;
					
					Trans_Update_Positions(Particles, id_j, r_trans_vec);
					Particles.Check_Periodic_CM(id_j, Box);
					Update_Periodic_Positions(Particles, Box, id_j);
					Particles.N_Particle[id_j]->Calculate_Axis();
					Particles.N_Particle[id_j]->Calculate_Patch_Position();	


					e_single = Calculate_Potential(id_j,j, Particles, Box);
				
					q_factor = exp(beta*(e_pair - e_single));
					q_factor = GSL_MIN(1,q_factor);
					q_r[fl] = q_factor;

				
					Reset_Positions(Particles, id_j);
					Particles.N_Particle[id_j]->Calculate_Axis();
					Particles.N_Particle[id_j]->Calculate_Patch_Position();

					N_failed_links = N_failed_links + 1;
					
					
				}



			}

			if(XI<b_factor){

					// Need to distinguish number of bonds and number of particles that go into list 	
				//////cout<<"hello"<<endl;
				
				//N_Bonds = cn;
				cn=N_Bonds;
				p_f[cn] = b_factor;
				//////////cout<<"b_factor "<<b_factor<<endl; 	
				// VIRTUAL opposite move

				// go with linker particle in oppositie direction

				Trans_Update_Positions(Particles, id_j, r_trans_vec);
				Particles.Check_Periodic_CM(id_j, Box);
				Update_Periodic_Positions(Particles, Box, id_j);
				Particles.N_Particle[id_j]->Calculate_Axis();
				Particles.N_Particle[id_j]->Calculate_Patch_Position();	
					
				// reset positions of linkee

				Reset_Positions(Particles, j);
				Particles.N_Particle[j]->Calculate_Axis();
				Particles.N_Particle[j]->Calculate_Patch_Position();

				// calculate single potential again and reverse linking probablility
				e_single = Calculate_Potential(id_j,j, Particles, Box);
				
				//////////cout<<" e_pair "<<e_pair<<endl;
				//////////cout<<" e_single "<<e_single<<endl;
				
				b_factor = 1 - exp(beta*(e_pair - e_single));
				b_factor = GSL_MAX(0,b_factor);
				p_r[cn] = b_factor;

				//////////cout<<" b_factor reverse "<<b_factor<<endl;
				// move linker back to new move


				Reset_Positions(Particles, id_j);
				Particles.N_Particle[id_j]->Calculate_Axis();
				Particles.N_Particle[id_j]->Calculate_Patch_Position();	
				


				//cn = cn+1;
				N_Bonds=N_Bonds+1;

				// move linkee back to new move 

				Trans_Pseudocluster_Recursion(j, cn, fl, Particles, Box);
			}	
			

		}		
	}	

}

/*
void move::Trans_Cluster_Move(particles& Particles, box* Box, fileio& Fileio, int mc_time) {


 	N_List = 0;
 	N_Bonds= 0;
 	N_failed_links=0;

	//choose particle id

	id =  gsl_rng_uniform_int(r,Box->N);
	//////////cout<<"id "<<id<<endl;

	Reset_Pseudo_Cluster(Box);

	// calculate translation 

	rand_x = gsl_ran_gaussian(r01, sigma_trans);
	rand_y = gsl_ran_gaussian(r01, sigma_trans);
	rand_z = gsl_ran_gaussian(r01, sigma_trans);
			 
	trans_vec.x = rand_x;
	trans_vec.y = rand_y;
			 
	if (is_2D==1){
		trans_vec.z = 0.0; 
	}
			 
	else{	 
		trans_vec.z = rand_z;
	}	

	
	r_trans_vec.x = -trans_vec.x;
	r_trans_vec.y = -trans_vec.y;
	r_trans_vec.z = -trans_vec.z;

	
	int links, failed_links;
	links=0;
	failed_links=0;

	
	Trans_Pseudocluster_Recursion(id, links, failed_links, Particles, Box);

	// end up with List[k] list of N_List elements in Pseudocluster
	// linking and reverse linking probabilities 

	

	// Collision Test

	

	for (int k=0;k<N_List;k++){
 	 
		 Reset_Positions(Particles, List[k]);
		 Trans_Update_Positions(Particles, List[k], trans_vec);
		 Particles.Check_Periodic_CM(List[k], Box);
		 Update_Periodic_Positions(Particles, Box, List[k]);
		 Particles.N_Particle[List[k]]->Calculate_Axis();
		 Particles.N_Particle[List[k]]->Calculate_Patch_Position();	
		
	}
	
	


	int k=0;
	//////cout<<"N_List   "<<N_List<<endl;
	//////cout<<"N_Bonds  "<<N_Bonds<<endl;

	exit_status = 0;
	col_count = 0;

	int Max_Length;
	Max_Length = rint(1./gsl_rng_uniform(r01));
	
	if (N_List > Max_Length){
		
		exit_status = 1;

	}

	
	do {
		 Particles.Collision_List[List[k]].Calculate_OP(Box, List[k], Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);		
		 //Collision Test 
		 Collision_Test( Particles, Box, List[k], Particles.Collision_List);
		 k++;

	}while((exit_status==0 )&& (k<N_List));


	if (exit_status > 0){

		//////////cout<<"collsion!"<<endl;
		//Reset all positions!
		for (int k=0;k<N_List;k++){
			Reset_Positions(Particles, List[k]);
			Particles.N_Particle[List[k]]->Calculate_Axis();
			Particles.N_Particle[List[k]]->Calculate_Patch_Position();
			Set_Positions(Particles, List[k]);	
	
		}	
	}

	
 	
 	if (exit_status == 0){

 		int list_j;
 		int list_k;
 		int m;
 		b_factor = 1;
 		bool overlapped;

 		//////cout<<"no collision"<<endl;
 		//Calculate new factors that are interacting in the new step	
 	
 	
 	

 		for (int k=0;k<N_List;k++){
 			
 			m = List[k];
 		

 			for(int cp1=0;cp1<Particles.Collision_List[m].Nm;cp1++){
 				overlapped=false;
 				list_j = Particles.Collision_List[m].Elements[cp1].nl_id;
 					
 				for(int cp2=0; cp2<Particles.Collision_List_old[m].Nm; cp2++){
					list_k = Particles.Collision_List_old[m].Elements[cp2].nl_id;
					
					if(list_j == list_k){
						overlapped=true;
					}

 				}

 				if (overlapped==false){
 					//Calculate boltzmann factor

 					//only new potential plays role
 				
 					
					e_new = Calculate_Potential(m,list_j, Particles, Box);
					b_factor = b_factor*exp(-beta*e_new);

 				}	


 			}	

 		}

 		//////cout<<"b_factor after new add "<<b_factor<<endl;

 		// Calculate factors for particels that were interacting in the last time step 

 		for (int k=0;k<N_List;k++){
 			
 			m = List[k];

 			for(int cp2=0;cp2<Particles.Collision_List_old[m].Nm;cp2++){
 				overlapped = false;
 				list_k = Particles.Collision_List_old[m].Elements[cp2].nl_id;
 				for(int cp1=0; cp1<Particles.Collision_List[m].Nm; cp1++){
					list_j = Particles.Collision_List[m].Elements[cp1].nl_id;
					
					if(list_j==list_k){
						overlapped=true;
					}

 				}
 				if (overlapped==false){
 					//Calculate boltzmann factor

 					//only old potential plays a role
 					Reset_Positions(Particles, m);
				 	Particles.N_Particle[m]->Calculate_Axis();
					Particles.N_Particle[m]->Calculate_Patch_Position();

					e_old = Calculate_Potential(m,list_k, Particles, Box);
					b_factor = b_factor*exp(beta*e_old);


					//move back 

					Trans_Update_Positions(Particles, m, trans_vec);	
					Particles.Check_Periodic_CM(m, Box);
					Update_Periodic_Positions(Particles, Box, m);
					Particles.N_Particle[m]->Calculate_Axis();
					Particles.N_Particle[m]->Calculate_Patch_Position();



 				}	


 			}	

 		}

 		

 		// Add up all the factors
 		double p_factor = 1;
 		double s_factor = 1;

 		double p_fr;
 		double q_fr;
 		int k = 0;
 		int l = 0;

 		////////cout<<"Boltzmann factor for newly and oldy bonded particles "<<b_factor<<endl;


 		do {
 			////////cout<<" p_r[k] "<<p_r[k]<<endl;
 			////////cout<<" p_f[k]" <<p_f[k]<<endl;

 			if (p_r[k]<1e-10){
 				p_fr = 0;
 				
 			}
 			
 			else{
 				p_fr = double(p_r[k])/double(p_f[k]); 

 			}	
 				
 			p_factor = p_factor*p_fr;
 			
 			k++;

 		}while(k<N_Bonds);

 		//////cout<<"Boltzmann factor for Pseudocluster "<<p_factor<<endl;

 		do {

 			if (q_r[l]<1e-10){
 				q_fr = 0; 
 			}	

 			else{
 				q_fr = double(q_r[l])/double(q_f[l]);
 			}

 			s_factor = s_factor*q_fr;
 			l++;


 		}while(l<N_failed_links);

 		//////cout<<"Boltzmann factor failed internal links "<<s_factor<<endl;

 		//////cout<<"p_factor "<<p_factor<<endl;

 		b_factor = GSL_MIN(1, b_factor*s_factor*p_factor);

 		//////cout<<" Overall Boltzmann factor "<<b_factor<<endl;

 		XI = gsl_rng_uniform(r01);
 	
 		// Cluster Move Acceptance criterium

 		if (b_factor<XI){

 			// Reset all positions
 			//////////cout<<"reset cluster move"<<endl;
 			//////////cout<<"Number of not moved particles "<<N_List<<endl;
 			for (int k=0;k<N_List;k++){
 				Reset_Positions(Particles, List[k]);
				Particles.N_Particle[List[k]]->Calculate_Axis();
				Particles.N_Particle[List[k]]->Calculate_Patch_Position();

 			}	
 		

 		}

 	
 		if (b_factor>=XI){


 			for (int k=0;k<N_List;k++){
 				Set_Positions(Particles, List[k]);	
 			}


 			Calculate_Pair_Potential(Particles, Box);
 			Particles.Total_Energy = Total_Energy;


 		}
			

 	}

 	Reset_Pseudo_Cluster(Box);
 	
	
}
*/


void move::Trans_Cluster_Move(particles& Particles, box* Box, fileio& Fileio, int mc_time) {

 	N_List = 0;
 	N_Bonds= 0;
 	N_failed_links=0;
 	double Total_Energy_old;
 	double delta_U;

	//choose particle id

	id =  gsl_rng_uniform_int(r,Box->N);
	//////////cout<<"id "<<id<<endl;

	Reset_Pseudo_Cluster(Box);

	// calculate translation 

	rand_x = gsl_ran_gaussian(r01, sigma_trans);
	rand_y = gsl_ran_gaussian(r01, sigma_trans);
	rand_z = gsl_ran_gaussian(r01, sigma_trans);
			 
	trans_vec.x = rand_x;
	trans_vec.y = rand_y;
			 
	if (is_2D==1){
		trans_vec.z = 0.0; 
	}
			 
	else{	 
		trans_vec.z = rand_z;
	}	
	
	int links, failed_links;
	links=0;
	failed_links=0;
	
	Pseudocluster_Recursion(id, links, failed_links, Particles, Box);

	// Collision Test


	for (int k=0;k<N_List;k++){
 	 
		 
		 Trans_Update_Positions(Particles, List[k], trans_vec);
		 Particles.Check_Periodic_CM(List[k], Box);
		 Update_Periodic_Positions(Particles, Box, List[k]);
		 Particles.N_Particle[List[k]]->Calculate_Axis();
		 Particles.N_Particle[List[k]]->Calculate_Patch_Position();	
		
	}
	
	int k=0;

	exit_status = 0;
	col_count = 0;

	int Max_Length;
	Max_Length = rint(1./gsl_rng_uniform(r01));
	
	if (N_List > Max_Length){
		
		exit_status = 1;

	}

	
	do {
		 Particles.Collision_List[List[k]].Calculate_OP(Box, List[k], Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);		
		 //Collision Test 
		 Collision_Test( Particles, Box, List[k], Particles.Collision_List);
		 k++;

	}while((exit_status==0 )&& (k<N_List));


	if (exit_status > 0){

		for (int k=0;k<N_List;k++){
			Reset_Positions(Particles, List[k]);
			Particles.N_Particle[List[k]]->Calculate_Axis();
			Particles.N_Particle[List[k]]->Calculate_Patch_Position();
			//Set_Positions(Particles, List[k]);	
	
		}	
	}

	
 	
 	if (exit_status == 0){

 	
 		Total_Energy_old = Particles.Total_Energy;
 		Calculate_Pair_Potential( Particles, Box);		
 		delta_U = Total_Energy - Total_Energy_old;
 		
 		
 		b_factor = GSL_MIN(1, exp((beta_f-beta)*delta_U));
 		XI = gsl_rng_uniform(r01);
 	
 	
 		if (b_factor<XI){

 			for (int k=0;k<N_List;k++){
 				int j;
 				j = List[k];
 				Reset_Positions(Particles, j);
				Particles.N_Particle[j]->Calculate_Axis();
				Particles.N_Particle[j]->Calculate_Patch_Position();

 			}	
 			
 			Particles.Total_Energy = Total_Energy_old;
 			Total_Energy = Total_Energy_old;


 		}
 	
 		if (b_factor>=XI){
 			//cout<<"accept"<<endl;
 			for (int k=0;k<N_List;k++){
 				int j;
 				j = List[k];
 				Set_Positions(Particles, j);	
 			}

 			Particles.Total_Energy = Total_Energy;

 		}
			
 	}

 	Reset_Pseudo_Cluster(Box);
 	
	
}



		
	
