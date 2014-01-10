
#include "move.h"


    move::move(){}
	
    double move::minimum(double a,double b){
			

			if(a<=b){
			return_value = a;
			}		

			if(b<a){
			return_value = b;
			}
				
			return return_value;	
					
	}	

	move::move(int size, int edge_length, double delta_tmax, double delta_rmax, double delta_Vmax, string is_translation_ON, string is_rotation_ON, string is_volumemove_ON){
	 
		L_axis = new m_vector[44];
		Cluster_List = new int[size];
			 
		for(int j=0;j<44;j++){
			  
			L_axis[j].x =0;
			L_axis[j].y =0;
			L_axis[j].z =0;
				 
		}	 
				 
			 
			 
		accept_translate = 0;
		accept_rotate = 0;
		accept_iso_vol = 0;
		accept_aniso_vol = 0;
		

		N_trans_moves = 0;
		N_rotate_moves = 0;
		N_iso_vol_moves = 0;
		N_aniso_vol_moves = 0;
			 
			 
			 
		sigma_trans = delta_tmax;
		dmax_V = delta_Vmax;
		kappa= delta_rmax; 
			
		cout<<"kappa"<<endl;	
			
		Unity[0] = 1.0;
		Unity[1] = 0.0;
		Unity[2] = 0.0;
			 
		Unity[3] = 0.0;
		Unity[4] = 1.0;
		Unity[5] = 0.0;
			 
		Unity[6] = 0.0;
		Unity[7] = 0.0;
		Unity[8] = 1.0;
			 
		zero_vec.x = 0.0;
		zero_vec.y = 0.0;
		zero_vec.z = 0.0;
			 
			  
		col_count = 0;
			 	 	  
		r = gsl_rng_alloc(gsl_rng_ranlxd2);
		r01 = gsl_rng_alloc(gsl_rng_ranlxd2); 

		for(int j=0;j<9;j++){
			Rot_mat[j] = Unity[j];
		}  
			  
		q_0 = 1.0;
		q_1 = 0.0;
		q_2 = 0.0,
		q_3 = 0.0;
			 
				
		T = 100.0;
		kB = 1.3806488*pow(10.0,-23.0);
		beta = 1./(kB*T);
			  
		
		//init discrete distribution for move_type 	 
		P_mt = new double[3];
		
		
		
		mt_w[0]=0;
		mt_w[1]=0;
		mt_w[2]=0;
		
		mt_sum = 0;
		
		cout<<"is_translation_ON "<<is_translation_ON<<endl;
		cout<<"is_rotation_ON "<<is_rotation_ON<<endl;
		cout<<"is_volumemove_ON "<<is_volumemove_ON<<endl; 
		
		
		if(is_translation_ON.compare("on")==0){
			mt_w[0]=size;
		}		
		
		if(is_rotation_ON.compare("on")==0){
			mt_w[1]=size;
		}	
		
		if(is_volumemove_ON.compare("on")==0){
			mt_w[2]=1;
		}	
		
		mt_sum = mt_w[0] + mt_w[1] + mt_w[2];
		cout<<"mt_sum: "<<mt_sum<<endl;
		
	
		P_mt[0] = double(mt_w[0])/double(mt_sum);
		P_mt[1] = double(mt_w[1])/double(mt_sum);
		P_mt[2] = double(mt_w[2])/double(mt_sum);
		
		sizep=3;
		 
		mt = gsl_ran_discrete_preproc(sizep, P_mt);	  
		
	
	}

	move::~move(){
		delete [] L_axis;
	}  

	void move::Set_Random_State(gsl_rng *r_start, gsl_rng *r01_start){		
		random_init_r = gsl_rng_memcpy (r, r_start);	
		random_init_r01 = gsl_rng_memcpy(r01, r01_start);
			 
	}  
    
	void move::Trans_Update_Positions(particles& Particles, int id, m_vector &trans_vec){
		
		Particles.N_Particle[id].x_center = Particles.N_Particle[id].x_center + trans_vec.x; 
		Particles.N_Particle[id].y_center = Particles.N_Particle[id].y_center + trans_vec.y; 
		Particles.N_Particle[id].z_center = Particles.N_Particle[id].z_center + trans_vec.z; 
			 
		// update positions  
			  
		for(int j=0;j<Particles.N_Particle[id].edge_N;j++){
					
			Particles.N_Particle[id].x[j] = Particles.N_Particle[id].x[j] + trans_vec.x; 
			Particles.N_Particle[id].y[j] = Particles.N_Particle[id].y[j] + trans_vec.y; 
			Particles.N_Particle[id].z[j] = Particles.N_Particle[id].z[j] + trans_vec.z; 
			  
		} 
		  
	}  


	void move::Rot_Update_Positions(particles& Particles, int id, double (&Rot_mat)[9]){
	   
		Particles.N_Particle[id].edges_from_center();
		Particles.N_Particle[id].distance_from_center();
	   
		for(int j=0;j<Particles.N_Particle[id].edge_N;j++){
		  	  
			Particles.N_Particle[id].new_dist_x[j] = Rot_mat[0]*Particles.N_Particle[id].dist_x[j] +Rot_mat[1]*Particles.N_Particle[id].dist_y[j] + Rot_mat[2]*Particles.N_Particle[id].dist_z[j]; 
			Particles.N_Particle[id].new_dist_y[j] = Rot_mat[3]*Particles.N_Particle[id].dist_x[j] +Rot_mat[4]*Particles.N_Particle[id].dist_y[j] + Rot_mat[5]*Particles.N_Particle[id].dist_z[j]; 
			Particles.N_Particle[id].new_dist_z[j] = Rot_mat[6]*Particles.N_Particle[id].dist_x[j] +Rot_mat[7]*Particles.N_Particle[id].dist_y[j] + Rot_mat[8]*Particles.N_Particle[id].dist_z[j];
				
				
			Particles.N_Particle[id].x[j] = Particles.N_Particle[id].x_center + Particles.N_Particle[id].new_dist_x[j]; 
			Particles.N_Particle[id].y[j] = Particles.N_Particle[id].y_center + Particles.N_Particle[id].new_dist_y[j];
			Particles.N_Particle[id].z[j] = Particles.N_Particle[id].z_center + Particles.N_Particle[id].new_dist_z[j]; 
				
		
				   
		}  

	}
	 	 
	int move::Random(int N){
	    id = gsl_rng_uniform_int(r,N);  
	    return id;
	}  

	void move::Iterate(particles& Particles, box* Box, fileio& Fileio, int mc_time){

		value=gsl_ran_discrete(r,mt);
	
	
		if(value==0){
							
			//cluster_counter=gsl_rng_uniform_int(r,cluster_size); 
			//id= Cluster_List[cluster_counter];
			id =  gsl_rng_uniform_int(r,Box->N);  
			Translate(Particles, Box, Fileio, id, mc_time);
					 
		}	
				
			
		if(value==1){	
	
			//cluster_counter=gsl_rng_uniform_int(r,cluster_size); 
			//id= Cluster_List[cluster_counter];
			id =  gsl_rng_uniform_int(r,Box->N);
			Rotate(Particles, Box, Fileio, id, mc_time);

		}		

		if(value==2){
				   
			Iso_Vol_Change (Particles, Box, Fileio, mc_time);
			//Aniso_Vol_Change (Box, Particles.N_Particle, Particles.N_Particle_old, Particles, Fileio, mc_time);
				
		}  	
	
		

    } 
		
	void move::Calculate_Acceptances(int mc_time){
		
		accept_translate_procent = double(accept_translate)/double(N_trans_moves);
		accept_rotate_procent = double(accept_rotate)/double(N_rotate_moves);
		accept_iso_vol_procent = double(accept_iso_vol)/double(N_iso_vol_moves);
		accept_complete_procent = double(accept_translate + accept_rotate + accept_iso_vol)/double(mc_time);
	}	

	void move::Calculate_Cluster_List( particles& Particles, box* Box){
		
		cluster_radius = Box->Lx/3.0;
		cluster_radius_square = (Box->Lx*Box->Lx)/9.0;
			
		cluster_counter=0;
			
		for(int id=0;id<Box->N;id++){			
					 Cluster_List[id] = -100;
		
	    }
			
			
		for(int id=0; id<Box->N;id++){
				
			diffc_x = Particles.N_Particle[id].x_center - Box->x_center;
			diffc_y = Particles.N_Particle[id].y_center - Box->y_center;
			diffc_z = Particles.N_Particle[id].z_center - Box->z_center;
				 	
			diffc_square = diffc_x*diffc_x + diffc_y*diffc_y + diffc_z*diffc_z;
					
			if(diffc_square<cluster_radius_square){
				
				Cluster_List[cluster_counter] = id;
				cluster_counter=cluster_counter+1;
							 				
			}
					
			cluster_size=cluster_counter;	
					
		}	
						
	}
		
	void move::Update_Periodic_Positions(particles& Particles, box* Box, int id){
		
		if(Particles.N_Particle[id].cm_out>=1){
				
			Particles.N_Particle[id].trans_periodic[0] =  double(Particles.N_Particle[id].cm_left_count)*Box->Lx  -  double(Particles.N_Particle[id].cm_right_count)*Box->Lx; 
			Particles.N_Particle[id].trans_periodic[1] =  double(Particles.N_Particle[id].cm_front_count)*Box->Ly -  double(Particles.N_Particle[id].cm_back_count)*Box->Ly;
			Particles.N_Particle[id].trans_periodic[2] =  double(Particles.N_Particle[id].cm_bottom_count)*Box->Lz - double(Particles.N_Particle[id].cm_top_count)*Box->Lz; 
				
		    	
			Particles.N_Particle[id].x_center = Particles.N_Particle[id].x_center + Particles.N_Particle[id].trans_periodic[0];
			Particles.N_Particle[id].y_center = Particles.N_Particle[id].y_center + Particles.N_Particle[id].trans_periodic[1];
			Particles.N_Particle[id].z_center = Particles.N_Particle[id].z_center + Particles.N_Particle[id].trans_periodic[2];
				 	   
			//update positions  
		
			for(int j=0;j<Particles.N_Particle[id].edge_N;j++){
						
				Particles.N_Particle[id].x[j] = Particles.N_Particle[id].x[j] + Particles.N_Particle[id].trans_periodic[0];
				Particles.N_Particle[id].y[j] = Particles.N_Particle[id].y[j] + Particles.N_Particle[id].trans_periodic[1];
				Particles.N_Particle[id].z[j] = Particles.N_Particle[id].z[j] + Particles.N_Particle[id].trans_periodic[2];  
				  
			} 
				    
		}
				 
		
	}	
	
	void move::Reset_Positions(particles& Particles, int id){
	
		Particles.N_Particle[id].x_center = Particles.N_Particle_old[id].x_center;
		Particles.N_Particle[id].y_center = Particles.N_Particle_old[id].y_center;
		Particles.N_Particle[id].z_center = Particles.N_Particle_old[id].z_center;
			
		for(int j=0;j<Particles.N_Particle[id].edge_N;j++){
			
			Particles.N_Particle[id].x[j] = Particles.N_Particle_old[id].x[j];
			Particles.N_Particle[id].y[j] = Particles.N_Particle_old[id].y[j];
			Particles.N_Particle[id].z[j] = Particles.N_Particle_old[id].z[j];
		
		}	
    
	}
	
	
	void move::Set_Positions(particles& Particles, int id){
	
		Particles.N_Particle_old[id].x_center = Particles.N_Particle[id].x_center;
		Particles.N_Particle_old[id].y_center = Particles.N_Particle[id].y_center;
	    Particles.N_Particle_old[id].z_center = Particles.N_Particle[id].z_center;
			
		for(int j=0;j<Particles.N_Particle[id].edge_N;j++){
				
			Particles.N_Particle_old[id].x[j] = Particles.N_Particle[id].x[j];
			Particles.N_Particle_old[id].y[j] = Particles.N_Particle[id].y[j];
			Particles.N_Particle_old[id].z[j] = Particles.N_Particle[id].z[j];
				
			
		}	

	
	}
	
	

	void move::Rot_Update_Quarternions_VON_MISES(particles& Particles, int id){
			
		
		mu.w= Particles.N_Particle[id].q.w;
		mu.x= Particles.N_Particle[id].q.x;
		mu.y= Particles.N_Particle[id].q.y;
		mu.z= Particles.N_Particle[id].q.z;
		
			
		//for max displacement: 0.2: kappa = 25
		//for max displacement: 0.07 kappa = 204
		//for max displacement: 2.0, kappa = 0.25
		//for max displacement: 1.0 kappa = 1.0
		
	
		b = -kappa + sqrt(kappa*kappa + 1);
		x0 = (1.0 - b)/(1.0 + b);
		c = kappa*x0 + 2*log(1 - x0*x0);	
		
		
		do{
			
			// Generate a random value with the beta(a, a) distribution with a = 1.5.
			//(1.5 because the parameters for beta are (m - 1)/2 and m = 4.

			do{
			  u = gsl_ran_flat(r01,-1.0,1.0);
			  v = gsl_rng_uniform(r01);

			  s = u*u + v*v;

			  }while(s>1.0);
			
			z = 0.5 + u*v*sqrt(1 - s)/double(s);
			
							
			u = gsl_rng_uniform(r01);
			w = (1.0 - (1.0 + b)*z)/(1.0 - (1.0 - b)*z);
			t = kappa*w + 2.0*log(1.0 - x0 * w) - c;
			
		}while(t<log(u));
				
			
		do {
			u0 = gsl_ran_flat(r01,-1.0,1.0);
			u1 = gsl_ran_flat(r01,-1.0,1.0);
						
			sq_u01 = u0*u0 + u1*u1;
				
		} while(sq_u01>1.0);
					
		sp_x = 2.*u0*sqrt(1.0-sq_u01);
		sp_y = 2.*u1*sqrt(1.0-sq_u01);
		sp_z = 1.0 - 2.*sq_u01;
		
		
		q_t.w = w;
		q_t.x = sp_x*sqrt(1.0 - w*w);
		q_t.y = sp_y*sqrt(1.0 - w*w); 
		q_t.z = sp_z*sqrt(1.0 - w*w);
			
		
		q_s.w = q_t.w*mu.w - q_t.x*mu.x - q_t.y*mu.y - q_t.z*mu.z;
		q_s.x = q_t.w*mu.x + q_t.x*mu.w + q_t.y*mu.z - q_t.z*mu.y; 
		q_s.y = q_t.w*mu.y - q_t.x*mu.z + q_t.y*mu.w + q_t.z*mu.x;
		q_s.z = q_t.w*mu.z + q_t.x*mu.y - q_t.y*mu.x + q_t.z*mu.w;
		
		Particles.N_Particle[id].q.w = q_s.w;
		Particles.N_Particle[id].q.x = q_s.x;
		Particles.N_Particle[id].q.y = q_s.y;
		Particles.N_Particle[id].q.z = q_s.z;
		
		double quat_norm;
		
		quat_norm = q_s.w*q_s.w + q_s.x*q_s.x + q_s.y*q_s.y + q_s.z*q_s.z;
		quat_norm = sqrt(quat_norm);
		
		/*
		ofstream orient_out("quarternion_von_mises.dat",  ios::out | ios::app);
		orient_out<<Particles.N_Particle[id].q.w<<"  "<<Particles.N_Particle[id].q.x<<"  "<<Particles.N_Particle[id].q.y<<"  "<<Particles.N_Particle[id].q.z<<"  "<<quat_norm<<endl;
		orient_out.close();
		*/
		
	
		
	}	
	
	

	
	void move::Rot_Update_Quarternions_RANDOM(particles& Particles, int id){
		
		
		 
		
		double scalar_q;
		 
		 
		//do {	
		
			do {
				u0 = gsl_ran_flat(r01, -1.0,1.0);
				u1 = gsl_ran_flat(r01,-1.0,1.0);
				
				 
				sq_u01 = u0*u0 + u1*u1;
				 
			 
			}while(sq_u01>1);
			 
			do {
				u2 = gsl_ran_flat(r01, -1.0,1.0);
				u3 = gsl_ran_flat(r01,-1.0,1.0);
			
				 
				sq_u23 = u2*u2 + u3*u3;
			 
			}while(sq_u23>1);
			 
			
			square_term = sqrt((1-sq_u01)/sq_u23);
		
			q_t.w = u0;
			q_t.x = u1;
			q_t.y= u2*square_term;
			q_t.z = u3*square_term;
			
			//scalar_q = q_t.x*Particles.N_Particle[id].q.x + q_t.y * Particles.N_Particle[id].q.y + q_t.z*Particles.N_Particle[id].q.z + q_t.w * Particles.N_Particle[id].q.w;
		
			
		//}while(scalar_q < 0.95);
			
			
		Particles.N_Particle[id].q.w = q_t.w;
		Particles.N_Particle[id].q.x = q_t.x;
		Particles.N_Particle[id].q.y = q_t.y;
		Particles.N_Particle[id].q.z = q_t.z;
			
	
    }	
	

	
	
