
#include "gear.h"
	
	
	gear::gear(){

		Ng=8;
		edge_N=8;
		N_all=Ng*8;
		
		N_independent_faces=12;
		
		edge_N_vec = 3;

		gear_wheels = new cube[Ng];
	
		 x = new double[edge_N];
		 y = new double[edge_N];
		 z = new double[edge_N];
		   
		 dist_x = new double[edge_N];
		 dist_y = new double[edge_N];
		 dist_z = new double[edge_N];
		 
		 new_dist_x = new double[edge_N];
		 new_dist_y = new double[edge_N];
		 new_dist_z = new double[edge_N];
		 
		 edges= new m_vector[edge_N_vec];
		 facenormal = new m_vector[N_independent_faces];
		 
		 edge_out = new int[edge_N];
		 
		 copy_count = 0;
		
		 trans_periodic = new double[3];
		
		 trans_old = new double[3];
		 Rot_old = new double[9];

		 rel_pos = new m_vector[edge_N];
			


	}

	 gear::~ gear(){
  

		 delete [] x; 
		 delete [] y; 
		 delete [] z; 
		 
		 delete [] dist_x; 
		 delete [] dist_y; 
		 delete [] dist_z; 
		 
		 delete [] edge_out;
		 
		 delete [] new_dist_x;
		 delete [] new_dist_y;
		 delete [] new_dist_z;
		 
		 delete[] trans_periodic;
		
		 
		 delete[] trans_old;
		 delete[] Rot_old;
		
  
        }  


	void gear::Set_Lengths(){

		phi=0;
		phi_old=phi;
		rel_pos = new m_vector[Ng];
		base_v = new m_vector[Ng];

		r=1.0;
		u=2*M_PI*r;
		Lx_max=u/double(Ng*2);
		p_Lx=1.0;

		Lx=Lx_max*p_Lx;
		d=(u-Ng*Lx)/double(Ng);
        
		for (int i=0;i<Ng;i++){
			gear_wheels[i].Lx = Lx;
		}

		theta=2.0*asin(Lx/r);

		As= (r*r*theta)/2.0;
		At= (h*Lx)/2.0;
		Ap= As-At;

		V=2.0*M_PI*r*Lx + Ng*(Lx*Lx*Lx - Ap*Lx);

		for (int i=0;i<Ng;i++){
			gear_wheels[i].Set_Lengths(Lx);
		}
		get_initial_relative_wheel_positions(Ng);
		get_relative_orientation();
		
	
		//R = r + Lx - h_t;
		R=r+Lx*sqrt(2.0);
		cut_off = 2.*R;
		cut_off_squared = cut_off*cut_off;
		cut_off_inner=2.*r;
		//cut_off_small=Lx*sqrt(2.0);
		cut_off_small=Lx*sqrt(3);
		cut_off_small_squared=cut_off_small*cut_off_small;



	}

	
	void gear::get_initial_relative_wheel_positions(int Ng){

		get_circle_partition(Ng);

		h   = sqrt(r*r - (Lx*Lx)/4.);
		h_t = r-h;
		ds  = r + Lx/2. - h_t;

		for(int i=0;i<Ng;i++){
			rel_pos[i].x = ds*base_v[i].x;
			rel_pos[i].y = ds*base_v[i].y;
			rel_pos[i].z = ds*base_v[i].z;


		}

	}

	void gear::get_circle_partition(int Ng){

		d_phi=(2.0*M_PI)/double(Ng);
		
		phi_n = 0;

		for (int i=0;i<Ng;i++){
			
			phi_n = phi_n + i*d_phi; 
			gear_wheels[i].phi= phi_n;
			
			base_v[i].x= cos(phi_n);
			base_v[i].y= sin(phi_n);
			base_v[i].z= 0.0;

		}

	}

	void gear::Update_Rotation(){

		phi_n = phi;

		for (int i=0;i<Ng;i++){
			
			phi_n = phi_n + i*d_phi; 

			if (phi_n > 2*M_PI){
				phi_n= phi_n - 2*M_PI;
			}	
			
			
		
			gear_wheels[i].phi= phi_n;
			
			rel_pos[i].x= ds*cos(phi_n);
			rel_pos[i].y= ds*sin(phi_n);
			rel_pos[i].z= ds*0.0;
			

		}

		//edges_from_center();

	}


	void gear::get_relative_orientation(){

		m_vector d_pos;

		for(int i=0;i<Ng;i++){

			gear_wheels[i].relative_distances();

			Rot_mat[0] = cos(gear_wheels[i].phi);
			Rot_mat[1] = -sin(gear_wheels[i].phi);
			Rot_mat[2] = 0;
		
			Rot_mat[3] =  sin(gear_wheels[i].phi);
			Rot_mat[4] =  cos(gear_wheels[i].phi);
			Rot_mat[5] =  0;
		
			Rot_mat[6] = 0;
			Rot_mat[7] = 0;
			Rot_mat[8] = 1;
		
			for(int k=0;k<gear_wheels[0].edge_N;k++){
			
				d_pos.x = Rot_mat[0]*gear_wheels[i].rel_pos[k].x +Rot_mat[1]*gear_wheels[i].rel_pos[k].y+ Rot_mat[2]*gear_wheels[i].rel_pos[k].z;
				d_pos.y = Rot_mat[3]*gear_wheels[i].rel_pos[k].x +Rot_mat[4]*gear_wheels[i].rel_pos[k].y+ Rot_mat[5]*gear_wheels[i].rel_pos[k].z;
				d_pos.z = Rot_mat[6]*gear_wheels[i].rel_pos[k].x +Rot_mat[7]*gear_wheels[i].rel_pos[k].y+ Rot_mat[8]*gear_wheels[i].rel_pos[k].z;

				gear_wheels[i].rel_pos[k].x = d_pos.x;
				gear_wheels[i].rel_pos[k].y = d_pos.y;
				gear_wheels[i].rel_pos[k].z = d_pos.z; 
				//gear_wheels[i].rel_pos[k].z = 0.0;



			}


		}


	}


	void gear::edges_from_center(){

		if(phi!=phi_old){
				Update_Rotation();
			
		}	
			
		for(int i=0;i<Ng;i++){
			gear_wheels[i].x_center = x_center + rel_pos[i].x;	
			gear_wheels[i].y_center = y_center + rel_pos[i].y;
			gear_wheels[i].z_center = z_center + rel_pos[i].z;
			
			x[i] = gear_wheels[i].x_center;
			y[i] = gear_wheels[i].y_center;
			z[i] = gear_wheels[i].z_center;

			get_relative_orientation();

			for(int k=0;k<gear_wheels[i].edge_N;k++){

				gear_wheels[i].x[k] = gear_wheels[i].x_center + gear_wheels[i].rel_pos[k].x;
				gear_wheels[i].y[k] = gear_wheels[i].y_center + gear_wheels[i].rel_pos[k].y;
				gear_wheels[i].z[k] = gear_wheels[i].z_center + gear_wheels[i].rel_pos[k].z;



			}


		}
		
	}


	
	void gear::Set_Start_Lattice_Position(int id, double box_Lx, double box_Ly, double box_Lz, int N_box){
		
		//for cubic lattice in 2D or anisotropic box shape
		
		int N_sitespx, N_sitespy;
		double N_sitespx_float, N_sitespy_float;
		double l_distpx, l_distpy;
		
		N_sitespx = rint(pow(double(N_box),1./2.));
		N_sitespx_float = double(N_sitespx);  
		N_sitespy = N_sitespx;
		N_sitespy_float = double(N_sitespy);  
		
		cout<<"Nsitespxy= "<<N_sitespx<<" "<<N_sitespy<<endl;
		
	
		//N_sitespz = box_Lz/Lx;
		//N_sitespz_float = double(N_sitespz);  
		
		  
		l_distpx = (box_Lx - N_sitespx_float)/N_sitespx_float;  
		l_distpy = (box_Ly - N_sitespy_float)/N_sitespy_float;   
		//l_distpz = (box_Lz - N_sitespz_float)/N_sitespz_float;   	 
		
		
		x_center = R/2.0 + l_distpx/2.0 + double(id %N_sitespx + l_distpx*(id % N_sitespx));
		y_center = R/2.0 + l_distpy/2.0 + double((id/N_sitespx)%N_sitespy + l_distpx*((id/N_sitespy)%N_sitespx)); 
		z_center = Lx/2.0;
		//z_center = 0.5 + l_distpz/2.0 + double(id/int(pow((double)N_sitespx,2)) + l_distpz*(id/int((double)pow(N_sitespx,2))));
		  	 
		
		edges_from_center();
	
		
	}

	void gear::write_positions(ofstream& fout){
			for(int i=0;i<Ng;i++){
				for(int k=0;k<gear_wheels[i].edge_N;k++){
					fout<<"zr"<<"      "<<gear_wheels[i].x[k]<<"   "<<gear_wheels[i].y[k]<<"   "<<gear_wheels[i].z[k]<<endl;	
			    }
			}					
	}	



	void gear::distance_from_center(){

	}
	void gear::Calculate_Axis(){


	} 
	void gear::Set_Axis(){

	}	

	void gear::Set_Start_Lattice_Position(int id, double box_Lx, int N_box){

	}

	

	double gear::Calculate_Projection_to_Separating_Axis(m_vector laxis){

		int Rp;
		Rp=0;

		return Rp;

	}


	void gear::Calculate_Face_Normals(){
		
	
		int k=0;
		double fnorm;
		
		for(int i=0;i<Ng;i++){
			
			facenormal[k].x= gear_wheels[i].x[1] - gear_wheels[i].x[0];
			facenormal[k].y= gear_wheels[i].y[1] - gear_wheels[i].y[0];
			facenormal[k].z= gear_wheels[i].z[1] - gear_wheels[i].z[0];
			
			fnorm=facenormal[k].norm();
			
			facenormal[k].x=facenormal[k].x/fnorm;
			facenormal[k].y=facenormal[k].y/fnorm;
			facenormal[k].z=facenormal[k].z/fnorm;
			
		  
			k++; 
		  
			facenormal[k].x = gear_wheels[i].x[3] - gear_wheels[i].x[0];
			facenormal[k].y = gear_wheels[i].y[3] - gear_wheels[i].y[0];
			facenormal[k].z = gear_wheels[i].z[3] - gear_wheels[i].z[0];
			
			fnorm=facenormal[k].norm();
			
			facenormal[k].x=facenormal[k].x/fnorm;
			facenormal[k].y=facenormal[k].y/fnorm;
			facenormal[k].z=facenormal[k].z/fnorm;
			
		 
			k++;
		
		
		}	
			
		
			
			
	}	

	/*
    void gear::Calculate_Face_Normals(){
		
		// face normals are the same as axis in the cube:
		
		facenormal[0].x = ax_1.x;
		facenormal[0].y = ax_1.y;
		facenormal[0].z = ax_1.z;
		
		facenormal[1].x = ax_2.x;
		facenormal[1].y = ax_2.y;
		facenormal[1].z = ax_2.z;
		
		facenormal[2].x = ax_3.x;
		facenormal[2].y = ax_3.y;
		facenormal[2].z = ax_3.z;
		
		
		edges[0].x = ax_1.x;
		edges[0].y = ax_1.y;
		edges[0].z = ax_1.z;
		
		edges[1].x = ax_2.x;
		edges[1].y = ax_2.y;
		edges[1].z = ax_2.z;
		
		edges[2].x = ax_3.x;
		edges[2].y = ax_3.y;
		edges[2].z = ax_3.z;
		
		//cout<<"Calculate Face Normals"<<endl;
		
    }		


	double gear::Calculate_Projection_to_Separating_Axis(m_vector laxis){
		
		double Rp;
	
		
		Rp = (fabs(ax_1.x*laxis.x + ax_1.y*laxis.y + ax_1.z*laxis.z) + 
			  fabs(ax_2.x*laxis.x + ax_2.y*laxis.y + ax_2.z*laxis.z) + 
			  fabs(ax_3.x*laxis.x + ax_3.y*laxis.y + ax_3.z*laxis.z) )*(Lx/2.0);
		
	
		
		return Rp;
		
		//cout<<"Calculate Projection to Separating Axis"<<endl;
		
		
	}	
		*/
