#include "order_parameter_2.h"


order_parameter::~order_parameter(){
	
	
}

order_parameter::order_parameter(box* Box, particles& Particles, int MAX_coll_p_in){
		
		
		//global
		
		MAX_coll_p = MAX_coll_p_in;
		cout<<"MAX_coll_p: order param "<<MAX_coll_p<<endl;
		
		
		//boost read parameters;
		
	
		boost::property_tree::ptree pt;
		boost::property_tree::ini_parser::read_ini("para.ini", pt);

		cutoff_coll_list = pt.get<double>("OP_general.Cutoff_collisionlist");
		cutoff_op= pt.get<double>("OP_general.Cutoff_op");
		
		cubic_op3_T = pt.get<double>("OP_cubic_op3.Cutoff");
		N_bins= pt.get<int>("OP_cubic_op3.Histogram_Nbins");
		xmin= pt.get<double>("OP_cubic_op3.Histogram_xmin");
		xmax= pt.get<double>("OP_cubic_op3.Histogram_xmax");
		cubic_op3_order= pt.get<int>("OP_cubic_op3.Order");
		
		qcorr_T_in =pt.get<double>("OP_q.Cutoff_q_corr");
		
		Cutoff_Type=pt.get<string>("OP_Cluster.Cutoff_Type");
		//Cutoff_Type="N_obonds";
		//strcpy(Cutoff_Type, "N_obonds");
		cout<<"Cutoff_Type: "<<Cutoff_Type<<endl;
		
		
		phase_T=pt.get<double>("OP_Cluster.Cutoff");
		
		max_r = pt.get<double>("g_radial.max_r");
		delta_r = pt.get<double>("g_radial.delta_r");
				
		
	
		// set variables for g_radial
		r_start = Particles.N_Particle[0].Lx;
		r_stop  = Particles.N_Particle[0].Lx*max_r;
		
		
		N_histo_points = int( ((r_stop - r_start)/(delta_r))); 
		g_radial_histo = new int[N_histo_points];
		g_radial_distribution = new data_point_2D[N_histo_points];
		
		gx_histo = new int[N_histo_points];
		gy_histo = new int[N_histo_points];
		gz_histo = new int[N_histo_points];
		
		gx_distribution = new data_point_2D[N_histo_points];
		gy_distribution = new data_point_2D[N_histo_points];
		gz_distribution = new data_point_2D[N_histo_points];
		
	
		//for cubic orderparameters
		N_obonds = new int[Box->N];
		N_obonds_cubic4 = new int[Box->N];
		cubic_op4_av = new double[Box->N];
		
		
		
		//spherical harmonics						
		q4.Set(4, Box->N, MAX_coll_p, qcorr_T_in);
		q6.Set(6, Box->N, MAX_coll_p, qcorr_T_in);
		
		is_nPhase = new int[Box->N];
		
	
		//q3op inital
		
		cubic_op3 = 0.0;
		cubic_op3_av = new double[Box->N];	

		

		histogram = gsl_histogram_alloc (N_bins);
		set_state = gsl_histogram_set_ranges_uniform (histogram, xmin, xmax);    



	}
	

	void order_parameter::Calculate_Local_Order_Parameters(particles& Particles, cube* N_Particle, box* Box){
	
		
			
			for( int id=0;id<Box->N; id++){
				Particles.Collision_List[id].Calculate(Box, id, Particles.Id_Cell_List, Particles.Cell_List, Particles.Cell, Particles.N_Particle, cutoff_coll_list, Particles.MAX_coll_p);	
			}	
			
			
			
			for( int id=0;id<Box->N; id++){
				//Particles.Collision_List[id].Calculate_Neighbours(cutoff_op);
				Particles.Collision_List[id].Calculate_Neighbours_Cubic(cutoff_op,  Particles.N_Particle, id);
			}
		
		/*		
		cout<<"start calculating q4"<<endl;
		q4.Calculate_local_op(Particles, Box);
		cout<<"finish calculating q4"<<endl;
		q4.Calculate_Corr_Funct( Particles, Box);
		*/
	
		Calculate_cubic_op3(Particles, Box);
		/*
		double nqbonds;
		
		for(int id=0;id<Box->N; id++){
			nqbonds= N_obonds[id];
			update_histogram(nqbonds);
		}
		*/
		
		
		
		
		/*
		double g_cutoff;
		g_cutoff=0;
		g_cutoff= Calculate_g_r(Particles, Box);
	    Calculate_g_xyz(Particles, Box);
		*/

	}	


/*
void order_param::g_r(){
		
	}	
	
void order_param::g_o(){
	}
void order_param::g_xyz(){
	}
  
*/
    

	void spher_harm_op::Set(int l, int N, int MAX_coll_p_in, double qcorr_T_in){
	
		l_op = l;
		
		m_total=2*l_op +1;

		MAX_coll_p = MAX_coll_p_in;
		cout<<"MAX_coll_p spher harm"<<MAX_coll_p<<endl;
		
		qcorr_T = qcorr_T_in;
		
		
			
		 Yij_m = new complex<double>**[N];
			for(int i=0;i<N;i++){
				Yij_m[i] = new complex<double>*[MAX_coll_p];	
					for(int j=0;j<MAX_coll_p;j++){
						Yij_m[i][j]= new complex<double>[m_total];	
					
				    }	 
				 
			    }
				
		q_id_m = new complex<double>*[N];
			for(int i=0;i<N;i++){
				q_id_m[i] = new complex<double>[m_total];	
					
			    }
		
		q_id_m_av= new complex<double>*[N];
			for(int i=0;i<N;i++){
				q_id_m_av[i] = new complex<double>[m_total];	
					
			    }
		

		ql = new double[N];
		
		N_qbonds = new int[N];
		q_id_norm = new double[N];
		
		
		ql_ij_scalar = new double*[N];
		for(int i=0;i<N;i++){
			ql_ij_scalar[i] = new double[MAX_coll_p];
		}
		
		qq_id_m = new complex<double>*[N];
		for(int i=0;i<N;i++){
			qq_id_m[i] = new complex<double>[m_total];
		}
		
		
	}
	
	
	void spher_harm_op::Calculate_local_op(particles& Particles, box* Box){
		       
		int N_count;       
		       
        m_total = 2*l_op+1;
    
		for(int id1=0; id1<Box->N;id1++){
			
			//cout<<MAX_coll_p<<endl;
			
			for(int j=0;j<MAX_coll_p;j++){
				
				for(int m=0;m<m_total;m++){					
					Yij_m[id1][j][m] = complex<double>(0.0,0.0);	
				}	
					
			}		
		}		
		
		
		for(int id=0; id<Box->N;id++){
					
			for(int j = 0;j<Particles.Collision_List[id].Nm;j++){
			
			   
				if((Particles.Collision_List[id].Elements[j].is_neighbour==1)&&(id!=Particles.Collision_List[id].Elements[j].nl_id)){	
					   
					 list_j = Particles.Collision_List[id].Elements[j].nl_id;
				
					 dist_cos_theta = (Particles.Collision_List[id].Elements[j].distance.z)/Particles.Collision_List[id].Elements[j].distance_norm;
					 dist_phi= atan((Particles.Collision_List[id].Elements[j].distance.y)/(Particles.Collision_List[id].Elements[j].distance.x + 1.0e-14));
						
					 for(int s=l_op; s<m_total;s++){					
						 int m;
						 m=s-l_op;
						 legendre_spher_harm = gsl_sf_legendre_sphPlm (l_op,m, dist_cos_theta);
						
						 Yij_m[id][j][s] = exp(complex<double>(0.0,m*dist_phi))*legendre_spher_harm;
											
						}
								
					  
								
					  for(int s=0;s<l_op;s++){	
						 int m;
						 m=s-l_op;
						 Yij_m[id][j][s] = gsl_pow_int(-1.0,m)*conj(Yij_m[id][j][s-2*(s-l_op)]);
						}
								
			    }
	        }
        }
    
		
		
    
		for(int id=0;id<Box->N; id++){
				  
			N_count = 0;
			
				  
			for(int m=0; m<m_total; m++){  	  
				q_id_m[id][m] = complex<double>(0.0,0.0);
				q_id_m_av[id][m] = complex<double>(0.0,0.0);
				
			
			}
				  
			for(int j=0;j<Particles.Collision_List[id].Nm; j++){
					
					if(Particles.Collision_List[id].Elements[j].is_neighbour==1){
						
						
					    N_count = N_count +1;
						
						for(int m=0; m<m_total; m++){  
								
							q_id_m[id][m] = q_id_m[id][m] +  Yij_m[id][j][m]; 
					
					} 
				}
			}	
			 
		 	
			for(int m=0; m<m_total; m++){  
					
				//q_id_m[id][m] = q_id_m[id][m]/((double(Particles.Collision_List[id].NNm-1)));
				q_id_m[id][m] = q_id_m[id][m]/double(N_count);
				
			}
				
			
					
			/*
			ql[id]=0.0;
					
			for(int m=0;m<m_total;m++){
						
				ql[id] = ql[id] + norm(q_id_m[id][m]);
		
					   
				}
						
			ql[id]= sqrt(((4.0*M_PI)/(double(m_total)))*ql[id]);	
					
			*/		
	  
	    }
				
				
		//for averaged local bond order parameters
					
		for(int id=0;id<Box->N; id++){
			N_count=0;
			
			for(int j=0;j<Particles.Collision_List[id].Nm; j++){
				if((Particles.Collision_List[id].Elements[j].is_neighbour==1)){	
					N_count = N_count + 1;
						
					for(int m=0; m<m_total; m++){  
				
						list_j = Particles.Collision_List[id].Elements[j].nl_id;				
						q_id_m_av[id][m] = q_id_m_av[id][m] + q_id_m[list_j][m];
						
					}
				}
			}		
				
			for(int m=0; m<m_total; m++){  
		    		//q_id_m_av[id][m] = q_id_m_av[id][m]/(double(Particles.Collision_List[id].NNm));
					q_id_m_av[id][m] = q_id_m_av[id][m]/double(N_count);
					
				
			}				
		}	
		
		
		
		for(int id=0;id<Box->N; id++){				
				
			ql[id]=0.0;
					
			for(int m=0;m<m_total;m++){
						
				ql[id] = ql[id] + norm(q_id_m_av[id][m]);
			   
			}
						
			ql[id]= sqrt(((4.0*M_PI)/(double(m_total)))*ql[id]);	
							
		}
			//calcualte global Q
			
			Q_global = 0.0;
			
		for(int id=0;id<Box->N; id++){
					
				Q_global = Q_global + ql[id];
				
			}		
				
	
    
    }	
    
    
    void spher_harm_op::Calculate_Corr_Funct(particles& Particles, box* Box){
		
		double N_count; 
		
		m_total = 2*l_op+1;
    
		for(int id1=0; id1<Box->N;id1++){
			
			
			for(int j=0;j<MAX_coll_p;j++){
				
				for(int m=0;m<m_total;m++){					
					Yij_m[id1][j][m] = complex<double>(0.0,0.0);	
				}	
					
			}		
		}		
		
		
		for(int id=0; id<Box->N;id++){
					
			for(int j = 0;j<Particles.Collision_List[id].Nm;j++){
			
			   
				if((Particles.Collision_List[id].Elements[j].is_neighbour==1)&&(id!=Particles.Collision_List[id].Elements[j].nl_id)){	
					   
					 list_j = Particles.Collision_List[id].Elements[j].nl_id;
				
					 dist_cos_theta = (Particles.Collision_List[id].Elements[j].distance.z)/Particles.Collision_List[id].Elements[j].distance_norm;
					 dist_phi= atan((Particles.Collision_List[id].Elements[j].distance.y)/(Particles.Collision_List[id].Elements[j].distance.x + 1.0e-14));
						
					 for(int s=l_op; s<m_total;s++){					
						 int m;
						 m=s-l_op;
						 legendre_spher_harm = gsl_sf_legendre_sphPlm (l_op,m, dist_cos_theta);
						
						 Yij_m[id][j][s] = exp(complex<double>(0.0,m*dist_phi))*legendre_spher_harm;
											
						}
								
					  
								
					  for(int s=0;s<l_op;s++){	
						 int m;
						 m=s-l_op;
						 Yij_m[id][j][s] = gsl_pow_int(-1.0,m)*conj(Yij_m[id][j][s-2*(s-l_op)]);
						 //cout<<"Yij_m "<<id<<"  "<<" "<<j<<" "<<s<<"  "<<Yij_m[id][j][s]<<endl;
						}
								
			    }
	        }
        }
    
	
		for(int id=0;id<Box->N; id++){
			
			N_count= 0;	  
				  
			for(int m=0; m<m_total; m++){  	  
				q_id_m[id][m] = complex<double>(0.0,0.0);
			 
				
			
			}
				  
			for(int j=0;j<Particles.Collision_List[id].Nm; j++){
					
					if(Particles.Collision_List[id].Elements[j].is_neighbour==1){
						
						N_count=N_count+1;
						
						for(int m=0; m<m_total; m++){  
								
							q_id_m[id][m] = q_id_m[id][m] +  Yij_m[id][j][m]; 
					
					} 
				}
			}	
			 
		 	
			for(int m=0; m<m_total; m++){  
					
				//q_id_m[id][m] = q_id_m[id][m]/((double(Particles.Collision_List[id].NNm-1)));
				q_id_m[id][m] = q_id_m[id][m]/double(N_count);
				
			}
	
		}
	
		for(int id = 0; id<Box->N;id++){
			
			q_id_norm[id] = 0;
			
			for( int m=0;m<m_total;m++){
				q_id_norm[id] = q_id_norm[id] + norm(q_id_m[id][m]);
			}
			
			q_id_norm[id] = sqrt(q_id_norm[id]);
			
			for(int m=0;m<m_total;m++){
				qq_id_m[id][m] = q_id_m[id][m]/q_id_norm[id];
				//cout<<"qq_id_m: "<<id<<" "<<m<<"  "<<qq_id_m[id][m]<<endl;
			}
			
		}
		
		for(int id=0; id<Box->N;id++){
			
			N_qbonds[id] = 0;
			
			for(int j = 0;j<Particles.Collision_List[id].Nm;j++){
				if((Particles.Collision_List[id].Elements[j].is_neighbour==1)&&(id!=Particles.Collision_List[id].Elements[j].nl_id)){	
					   
					list_j = Particles.Collision_List[id].Elements[j].nl_id;
				
					ql_ij_scalar[id][j] = 0.0;
					qscalar_id_j = complex<double>(0.0,0.0);
					
					for(int m = 0; m<m_total; m++){		
						qscalar_id_j = qscalar_id_j + qq_id_m[id][m]*conj(qq_id_m[list_j][m]);
						
						
					}
				
					ql_ij_scalar[id][j] = fabs(real(qscalar_id_j)); 
					//cout<<ql_ij_scalar[id][j]<<endl;
				
					if(ql_ij_scalar[id][j]>qcorr_T){
						
						N_qbonds[id] = N_qbonds[id] + 1;	
						
					}
		
				}
			}
			
		}
	}
    
    
	void order_parameter::Caculate_is_nPhase(box* Box){
	
	
		for(int id=0; id<Box->N;id++){ 
		 is_nPhase[id] = 0;
		}
		
		int l;
		l=0;
		
		if (Cutoff_Type.compare("N_obonds")==0){
			
			cout<<"hello N_obonds"<<endl;
			
			for(int id=0; id<Box->N;id++){ 
				if(N_obonds[id]>=phase_T){		
					l = l+1;
				 
					is_nPhase[id] = 1;
				}
			}
		}	
			
		if (Cutoff_Type.compare("N_qbonds")==0){
			
			for(int id=0; id<Box->N;id++){ 
				if(q4.N_qbonds[id]>=phase_T){		
					l = l+1;
				 
					is_nPhase[id] = 1;
				}
			}
		}		
			
			
		
		if (Cutoff_Type.compare("q6")==0){
			for(int id=0; id<Box->N;id++){ 
				if(q6.ql[id]>=phase_T){		
					l = l+1;
				 
					is_nPhase[id] = 1;
				}
			}
		}		
		
		if (Cutoff_Type.compare("q4")==0){
			for(int id=0; id<Box->N;id++){ 
				if(q4.ql[id]>=phase_T){		
					l = l+1;
				 
					is_nPhase[id] = 1;
				}
			}
		}		
		
		
		if (Cutoff_Type.compare("cubic_op4_av")==0){
			for(int id=0; id<Box->N;id++){ 
				if(cubic_op4_av[id]>=phase_T){		
					l = l+1;
				 
					is_nPhase[id] = 1;
				}
			}
		}		
		
		
		
		
				
		cout<<"Number of particles in new phase: "<<l<<endl;
	}	


	void order_parameter::Calculate_cubic_op3(particles& Particles, box* Box){
		
		
		for( int id =0; id<Box->N;id++ ){
			Particles.N_Particle[id].Calculate_Axis();	
		}


		for(int id=0; id<Box->N;id++){
			
			cubic_op3_av[id] = 0.0;
			
			a1x = Particles.N_Particle[id].ax_1.x;
			a1y = Particles.N_Particle[id].ax_1.y;
			a1z = Particles.N_Particle[id].ax_1.z;
					
			a2x = Particles.N_Particle[id].ax_2.x;
			a2y = Particles.N_Particle[id].ax_2.y;
			a2z = Particles.N_Particle[id].ax_2.z;
					
			a3x = Particles.N_Particle[id].ax_3.x;
			a3y = Particles.N_Particle[id].ax_3.y;
			a3z = Particles.N_Particle[id].ax_3.z;
			
			N_obonds[id] = 0;
			
			int N_count = 0;
			
			for(int j = 0;j<Particles.Collision_List[id].Nm;j++){
				if((Particles.Collision_List[id].Elements[j].is_neighbour==1)&&(id!=Particles.Collision_List[id].Elements[j].nl_id)){	

					N_count = N_count +1;
		
					cubic_op3 = 0;	

					list_j = Particles.Collision_List[id].Elements[j].nl_id;	
					
					b1x = Particles.N_Particle[list_j].ax_1.x;
					b1y = Particles.N_Particle[list_j].ax_1.y;
					b1z = Particles.N_Particle[list_j].ax_1.z;
							
					b2x = Particles.N_Particle[list_j].ax_2.x;
					b2y = Particles.N_Particle[list_j].ax_2.y;
					b2z = Particles.N_Particle[list_j].ax_2.z;
							
					b3x = Particles.N_Particle[list_j].ax_3.x;
					b3y = Particles.N_Particle[list_j].ax_3.y;
					b3z = Particles.N_Particle[list_j].ax_3.z;
					
					scp[0] = a1x*b1x + a1y*b1y + a1z*b1z;
					scp[1] = a1x*b2x + a1y*b2y + a1z*b2z;
					scp[2] = a1x*b3x + a1y*b3y + a1z*b3z;
					
					scp[3] = a2x*b1x + a2y*b1y + a2z*b1z;
					scp[4] = a2x*b2x + a2y*b2y + a2z*b2z;
					scp[5] = a2x*b3x + a2y*b3y + a2z*b3z;
					
					scp[6] = a3x*b1x + a3y*b1y + a3z*b1z;
					scp[7] = a3x*b2x + a3y*b2y + a3z*b2z;
					scp[8] = a3x*b3x + a3y*b3y + a3z*b3z;
					
					
					if(cubic_op3_order==6){
					
						for(int p=0; p<9; p++){
							cubic_op3 = cubic_op3 + scp[p]*scp[p]*scp[p]*scp[p]*scp[p]*scp[p];
						}
					
					}
					
					if(cubic_op3_order==4){
							for(int p=0; p<9; p++){
							cubic_op3 = cubic_op3 + scp[p]*scp[p]*scp[p]*scp[p];
						}
				
					}
					
					
					
					if(cubic_op3 > cubic_op3_T){
						N_obonds[id] = N_obonds[id] + 1;	
					}
					
				
					cubic_op3_av[id] = cubic_op3_av[id] + cubic_op3;
					
					
				}
			}	
		
			//cubic_op3_av[id] = cubic_op3_av[id]/double(Particles.Collision_List[id].NNm-1);
			cubic_op3_av[id] = cubic_op3_av[id]/double(N_count);
			//update_histogram(cubic_op3_av[id]);
		
		}
		

	}

	/*
	int Calculate_Dominant_Axis(int id, m_vector b){
		
		int which_axes;
		
		m_vector diff_vec[8];
		double diff_Nv[8];
		double surface_sum[6];
		double surface_min;
		int ax_N;
		
		
		for(int i=0;i<8;i++){
		
			diff_vec[i].x = Particles.N_Particle[id].x[i] - b.x;
			diff_vec[i].y = Particles.N_Particle[id].y[i] - b.y;
			diff_vec[i].z = Particles.N_Particle[id].z[i] - b.z;
		
			diff_Nv[i] = diff_vec[i].norm();
				
		}	
		
		surface_sum[0] = diff_Nv[0] + diff_Nv[1] + diff_Nv[2] + diff_Nv[3];
		surface_sum[1] = diff_Nv[4] + diff_Nv[5] + diff_Nv[6] + diff_Nv[7];
		surface_sum[2] = diff_Nv[0] + diff_Nv[1] + diff_Nv[4] + diff_Nv[5];
		surface_sum[3] = diff_Nv[2] + diff_Nv[3] + diff_Nv[6] + diff_Nv[7];
		surface_sum[4] = diff_Nv[0] + diff_Nv[3] + diff_Nv[4] + diff_Nv[7];
		surface_sum[5] = diff_Nv[1] + diff_Nv[2] + diff_Nv[5] + diff_Nv[6];
		
		surface_min = 1e10;
		
		for(int k=0;k<6;k++){
			if(surface_sum[k]<surface_min){
				surface_min=surface_sum[k];
				ax_N = k;
				
			}
		}	
		
		if ((ax_N == 0)||(ax_N == 1)){	
			which_axes=3;
		}	
		
		if ((ax_N == 2)||(ax_N == 3)){	
			which_axes=2;
		}	
		
		if ((ax_N == 4)||(ax_N == 5)){	
			which_axes=1;
		}	
		
		
		return which_axes;
	}





	void order_parameter::Calculate_cubic_op4(int Np){
		
		
		m_vector a1;
		m_vector b1;
		m_vector c1;
		
		m_vector a2;
		m_vector b2;
		m_vector c2;
	
		m_vector d;
		
		m_vector dom_ax1;
		m_vector dom_ax2;
		
		double d_norm;
		
		double scpd;
		
		
		for( int id =0; id<Np;id++ ){
			Particles.N_Particle[id].Calculate_Axis();	
		}


		for(int id=0; id<Np;id++){
			
			cubic_op4_av[id] = 0.0;
			
			a1.x = Particles.N_Particle[id].ax_1.x;
			a1.y = Particles.N_Particle[id].ax_1.y;
			a1.z = Particles.N_Particle[id].ax_1.z;
					
			b1.x = Particles.N_Particle[id].ax_2.x;
			b1.y = Particles.N_Particle[id].ax_2.y;
			b1.z = Particles.N_Particle[id].ax_2.z;
					
			c1.x = Particles.N_Particle[id].ax_3.x;
			c1.y = Particles.N_Particle[id].ax_3.y;
			c1.z = Particles.N_Particle[id].ax_3.z;
			
			N_obonds_cubic_4[id] = 0;
			
			int N_count= 0;
			
			for(int j = 0;j<Collision_List[id].Nm;j++){
				if((Collision_List[id].Elements[j].is_neighbour==1)&&(id!=Collision_List[id].Elements[j].nl_id)){	

				
					N_count=N_count+1;

					cubic_op4 = 0;	

					list_j = Collision_List[id].Elements[j].nl_id;	
					
					a2.x = Particles.N_Particle[list_j].ax_1.x;
					a2.y = Particles.N_Particle[list_j].ax_1.y;
					a2.z = Particles.N_Particle[list_j].ax_1.z;
							
					b2.x = Particles.N_Particle[list_j].ax_2.x;
					b2.y = Particles.N_Particle[list_j].ax_2.y;
					b2.z = Particles.N_Particle[list_j].ax_2.z;
							
					c2.x = Particles.N_Particle[list_j].ax_3.x;
					c2.y = Particles.N_Particle[list_j].ax_3.y;
					c2.z = Particles.N_Particle[list_j].ax_3.z;
					
					
					d.x = Particles.N_Particle[list_j].x_center - Particles.N_Particle[id].x_center;
                    d.y = Particles.N_Particle[list_j].y_center - Particles.N_Particle[id].y_center;
                    d.z = Particles.N_Particle[list_j].z_center - Particles.N_Particle[id].z_center;
                    
                    d.x = d.x -  Box.Lx*rint(d.x/Box.Lx);
                    d.y = d.y -  Box.Ly*rint(d.y/Box.Ly);
                    d.z = d.z -  Box.Lz*rint(d.z/Box.Lz);
                    
                
					d_norm = d.norm();
					d.x = d.x/d_norm;
					d.y = d.y/d_norm;
					d.z = d.z/d_norm;
				
					dom_ax1.x = 0;
					dom_ax1.y = 0;
					dom_ax1.z = 0;
					
					dom_ax2.x = 0;
					dom_ax2.y = 0;
					dom_ax2.z = 0;	
					
					
					int which_axes;
					
					which_axes=Calculate_Dominant_Axis(id,d);
					
					if(which_axes == 1){		
					  dom_ax1.x = a1.x;
					  dom_ax1.y = a1.y;
					  dom_ax1.z = a1.z;
					  //cout<<"Dominant Axes is a"<<endl
					}
				
					if(which_axes == 2){	
					  dom_ax1.x = b1.x;
					  dom_ax1.y = b1.y;
					  dom_ax1.z = b1.z;
					  //cout<<"Dominant Axes is b"<<endl
					}
				
					if(which_axes == 3){	
					  dom_ax1.x = c1.x;
					  dom_ax1.y = c1.y;
					  dom_ax1.z = c1.z;
					  //cout<<"Dominant Axes is c"<<endl;
					}
				
				
					which_axes=Calculate_Dominant_Axis(list_j,d);
					
					if(which_axes == 1){		
					  dom_ax2.x = a2.x;
					  dom_ax2.y = a2.y;
					  dom_ax2.z = a2.z;
					  //cout<<"Dominant Axes is a"<<endl
					}
				
					if(which_axes == 2){	
					  dom_ax2.x = b2.x;
					  dom_ax2.y = b2.y;
					  dom_ax2.z = b2.z;
					  //cout<<"Dominant Axes is b"<<endl
					}
				
					if(which_axes == 3){	
					  dom_ax2.x = c2.x;
					  dom_ax2.y = c2.y;
					  dom_ax2.z = c2.z;
					  //cout<<"Dominant Axes is c"<<endl;
					}
				
						
						scpd = scalar_product(dom_ax1, dom_ax2);
						scpd = fabs(2*scpd*scpd -1.0);
						
						cubic_op4 = scpd;
				
					
					//cubic_op4 = scp_max*scp_max*scp_max*scp_max*scp_max*scp_max;
						
					if(cubic_op4 > cubic_op4_T){
						N_obonds_cubic_4[id] = N_obonds_cubic_4[id] + 1;	
					}
					
				
					cubic_op4_av[id] = cubic_op4_av[id] + cubic_op4;
					
					
				}
			}	
		
			cubic_op4_av[id] = cubic_op4_av[id]/double(N_count);
			
			N_obonds_cubic_4[id] = N_obonds_cubic_4[id]/double(N_count);
			
		    //update_histogram(cubic_op4_av[id]);
		   
		   
		}
		

	}


	*/



	double order_parameter::Calculate_g_r(particles& Particles, box* Box){
	
		for(int p=0;p<N_histo_points;p++){
			
			g_radial_histo[p]=0; 
		
		}
	
		for(int id=0;id<Box->N;id++){
			
			
				for(int j=0;j<Particles.Collision_List[id].Nm;j++){
						
					num_r1 = r_start;	
						
					for(int p=0;p<N_histo_points;p++){
					 
						num_r2 = r_start + double(delta_r*(p+1));
						
						r_dist = Particles.Collision_List[id].Elements[j].distance_norm;
						 
					
						if((r_dist>=num_r1)&&(r_dist<=num_r2)&&(id!=Particles.Collision_List[id].Elements[j].nl_id)){
							
						   g_radial_histo[p] = g_radial_histo[p] + 1;	
							
						}
					 
						num_r1 = num_r2;
					 
					}
					
			    
			    }
	
	    } 

		rn = r_start;

		for(int p = 0;p<N_histo_points;p++){
		
			rn = r_start + double(delta_r*p);
			

			g_norm = 4.0*M_PI*rn*rn*delta_r*double(Box->N)*(double(Box->N)/Box->V);
			
			g_radial_distribution[p].x = rn;
			g_radial_distribution[p].y = double(g_radial_histo[p])/g_norm;
			
			cout<<"rn: "<<rn<<endl;
			cout<<"g_radial_distribution[p].y:"<<g_radial_distribution[p].y<<endl;
			
			
		}	
	
		g_min = 10.0;
	
		for(int p=0;p<N_histo_points;p++){
			
			if(g_radial_distribution[p].y<g_min){
				
				p_min = p;
				g_min = g_radial_distribution[p].y;			
		    }
		    
					
		}
	
	
        g_Cut_Off = g_radial_distribution[p_min].x;
	
		
		
		ofstream gr_out("g_r.dat", ios::out | ios::app);
		
			for(int p = 0;p<N_histo_points;p++){
		
				gr_out<<g_radial_distribution[p].x<<"   "<<g_radial_distribution[p].y<<endl;
			}
			
		gr_out.close();	
	
	
	
	    
		return g_Cut_Off;
	
	
    }

	
	void order_parameter::Calculate_g_xyz(particles& Particles, box* Box){
	
		for(int p=0;p<N_histo_points;p++){
			
			gx_histo[p]=0; 
			gy_histo[p]=0;
			gz_histo[p]=0; 
			
			
		
		}
	
		for(int id=0;id<Box->N;id++){
			
			
				for(int j=0;j<Particles.Collision_List[id].Nm;j++){
						
					num_r1 = r_start;	
					
					Particles.N_Particle[id].Calculate_Axis();
						
					nbd.x = Particles.Collision_List[id].Elements[j].distance.x; 
					nbd.y = Particles.Collision_List[id].Elements[j].distance.y; 
					nbd.z = Particles.Collision_List[id].Elements[j].distance.z; 
						
					//calculate distance.x, distance.y, distance.z in body centered coordinate system of particle id	
					newcoord_distance.x = Particles.N_Particle[id].ax_1.x*nbd.x + Particles.N_Particle[id].ax_2.x*nbd.y + Particles.N_Particle[id].ax_3.x*nbd.z;
					newcoord_distance.y = Particles.N_Particle[id].ax_1.y*nbd.x + Particles.N_Particle[id].ax_2.y*nbd.y + Particles.N_Particle[id].ax_3.y*nbd.z;
					newcoord_distance.z = Particles.N_Particle[id].ax_1.z*nbd.x + Particles.N_Particle[id].ax_2.z*nbd.y + Particles.N_Particle[id].ax_3.z*nbd.z;
						
					rv_dist.x = fabs(newcoord_distance.x);
					rv_dist.y = fabs(newcoord_distance.y);
					rv_dist.z = fabs(newcoord_distance.z);
					
					
					for(int p=0;p<N_histo_points;p++){
					 
						num_r2 = r_start + double(delta_r*(p+1));
					
						
						if((rv_dist.x>=num_r1)&&(rv_dist.x<=num_r2)&&
						   //(rv_dist.y>=0)&&(rv_dist.y<=0.25)&&
						   //(rv_dist.z>=0)&&(rv_dist.z<=0.25)&&
						   (id!=Particles.Collision_List[id].Elements[j].nl_id)){		
							  gx_histo[p] = gx_histo[p] + 1;	
						}	
					
					
					
					
					
						if((rv_dist.y>=num_r1)&&(rv_dist.y<=num_r2)&&(id!=Particles.Collision_List[id].Elements[j].nl_id)){		
							  gy_histo[p] = gy_histo[p] + 1;	
						}	
						
					
						if((rv_dist.z>=num_r1)&&(rv_dist.z<=num_r2)&&(id!=Particles.Collision_List[id].Elements[j].nl_id)){		
							  gz_histo[p] = gz_histo[p] + 1;	
						}	
						
					 
						num_r1 = num_r2;
					 
					}
					
			    
			    }
	
	    } 

		rn = r_start;

		for(int p = 0;p<N_histo_points;p++){
		
			rn = r_start + double(delta_r*p);
			g_norm = double(Box->V)/double(Box->N*Box->N);
			cout<<"g_norm"<<endl;
			
			gx_distribution[p].x = rn;
			gx_distribution[p].y = double(gx_histo[p])*g_norm;
			
			gy_distribution[p].x = rn;
			gy_distribution[p].y = double(gy_histo[p])*g_norm;
			
			gz_distribution[p].x = rn;
			gz_distribution[p].y = double(gz_histo[p])*g_norm;
			
			
			
		}	
	

			
		ofstream gx_out("gx.dat", ios::out | ios::app);
		
			for(int p = 0;p<N_histo_points;p++){
		
				gx_out<<gx_distribution[p].x<<"   "<<gx_distribution[p].y<<endl;
			}
			
		gx_out.close();	
		
		
		ofstream gy_out("gy.dat", ios::out | ios::app);
		
			for(int p = 0;p<N_histo_points;p++){
		
				gy_out<<gy_distribution[p].x<<"   "<<gy_distribution[p].y<<endl;
			}
			
		gy_out.close();	
		
		
		ofstream gz_out("gz.dat", ios::out | ios::app);
		
			for(int p = 0;p<N_histo_points;p++){
		
				gz_out<<gz_distribution[p].x<<"   "<<gz_distribution[p].y<<endl;
			}
			
		gz_out.close();	
	
	



    }



	void order_parameter::unset_histogram(){
		
		gsl_histogram_free (histogram);			
			
	}


	void order_parameter::update_histogram(double hv){
	
		if((hv >=xmin)&&(hv<=xmax)){
			update_state = gsl_histogram_increment (histogram, hv);
			
		}
	
	}	

	void order_parameter::write_histogram_info(){
		
		FILE *fp_h;
		fp_h = fopen("histogram.info","w");
		save_state = gsl_histogram_fprintf (fp_h, histogram, "%f", "%f");	
		close_state = fclose(fp_h);
		
			
	}	


