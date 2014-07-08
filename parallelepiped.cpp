#include "parallelepiped.h"



	parallelepiped::parallelepiped(){
 
		 edge_N = 8; 
		 N_independent_faces=3;
		 N_cross_edges = 3;
		 edge_N_vec = 3;
 
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
		
		
		
  
        }  


    parallelepiped::~ parallelepiped(){
  

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



    void parallelepiped::edges_from_center(){ // only valid if cubes are aligned with coordinate axes!
  


		 x[0] = x_center + (-Lx -2*a_x)/2.0;
		 y[0] = y_center + (-hb - a_x)/2.0;
		 z[0] = z_center + (-hz)/2.0; 

		 x[1] = x_center + (Lx-2*a_x)/2.0;
		 y[1] = y_center + (-hb - a_x)/2.0;
		 z[1] = z_center + (-hz)/2.0;

		 x[2] = x_center + (Lx)/2.0;
		 y[2] = y_center + (hb - a_x)/2.0;
		 z[2] = z_center + (-hz)/2.0;

		 x[3] = x_center + (-Lx)/2.0;
		 y[3] = y_center + (hb - a_x)/2.0;
		 z[3] = z_center + (-hz)/2.0;

		 x[4] = x_center + (-Lx)/2.0;
		 y[4] = y_center + (-hb + a_x)/2.0;
		 z[4] = z_center + (hz)/2.0;

		 x[5] = x_center + (Lx)/2.0;
		 y[5] = y_center + (-hb + a_x)/2.0;
		 z[5] = z_center + (hz)/2.0;

		 x[6] = x_center + (Lx+2*a_x)/2.0;
		 y[6] = y_center + (hb + a_x)/2.0;
		 z[6] = z_center + (hz)/2.0;

		 x[7] = x_center + (-Lx+2*a_x)/2.0;
		 y[7] = y_center + (hb + a_x)/2.0;
		 z[7] = z_center + (hz)/2.0;



        }  


	void parallelepiped::distance_from_center(){
 
	    for(int j=0;j<edge_N;j++){
	   
		     dist_x[j] =  x[j] - x_center;
		     dist_y[j] =  y[j] - y_center;
		     dist_z[j] =  z[j] - z_center;
	  
	    }
	    
	    //cout<<"Distance from center"<<endl;
	    
    }
    
    


	void parallelepiped::Calculate_Axis(){
		
		double norm_ax;
		
		 ax_1.x =  x[1] - x[0];
		 ax_1.y =  y[1] - y[0];
	     ax_1.z =  z[1] - z[0];
			
		norm_ax = ax_1.norm();
				 
	     ax_1.x = ax_1.x/norm_ax;
	     ax_1.y = ax_1.y/norm_ax; 
		 ax_1.z = ax_1.z/norm_ax;  
				 
		 ax_2.x =  x[3] - x[0];
		 ax_2.y =  y[3] - y[0];
		 ax_2.z =  z[3] - z[0];
				 
		norm_ax = ax_2.norm();		 
				 
		 ax_2.x = ax_2.x/norm_ax; 
		 ax_2.y = ax_2.y/norm_ax; 
		 ax_2.z = ax_2.z/norm_ax;
					 
		 ax_3.x =  x[4] - x[0];
	     ax_3.y =  y[4] - y[0];
	     ax_3.z =  z[4] - z[0];
		
		norm_ax = ax_3.norm();		 
					
		 ax_3.x = ax_3.x/norm_ax;
		 ax_3.y = ax_3.y/norm_ax; 
		 ax_3.z = ax_3.z/norm_ax;  
		
	
	}

	
	void parallelepiped::Set_Axis(){

	     ax_1_old.x =  ax_1.x;
		 ax_1_old.y =  ax_1.y;
		 ax_1_old.z =  ax_1.z;
			 	 
		 ax_2_old.x =  ax_2.x;
		 ax_2_old.y =  ax_2.y;
		 ax_2_old.z =  ax_2.z;
			
	 	 ax_3_old.x =  ax_3.x;
		 ax_3_old.y =  ax_3.y;
		 ax_3_old.z =  ax_3.z;

		//cout<<"Set Axis"<<endl;			
			
    } 	



    void parallelepiped::Set_Lengths(){
		
		Lx = 1.0;
		Ly = 1.0;
		Lz = 1.0;
		
		alpha=80.0*M_PI/180.0;
		gamma = 40.0*M_PI/180.0;
		
		alpha_2 = M_PI - alpha;
		gamma_2 = M_PI - gamma;
		
		
		hb=Ly*sin(alpha);
		hb_2= double(hb)/2.0;
		
		hz=Lz*sin(gamma);
		hz_2=double(hz)/2.0;
		
		
		a_x = sqrt(Ly*Ly - hb*hb);
		
		diag2_short = sqrt( (Lx-a_x)*(Lx-a_x) + hb*hb);
		diag2_long = sqrt( (Lx+a_x)*(Lx+a_x) + hb*hb);
		
		diag3_short = sqrt( diag2_short*diag2_short + Lz*Lz);
		diag3_long = sqrt( diag2_long*diag2_long + Lz*Lz);
		
		cut_off = diag3_long;
		cut_off_squared = diag3_long*diag3_long;
		

		V = Lx*hb*hz;
	
	

	}

	void parallelepiped::Set_Lengths(double Lx, double Ly, double Lz, double alpha, double gamma){
		
		
	
		alpha_2 = M_PI - alpha;
		gamma_2 = M_PI - gamma;
		
		
		hb=Ly*sin(alpha);
		hb_2= double(hb)/2.0;
		
		hz=Lz*sin(alpha);
		hz_2=double(hb)/2.0;
		
		
		a_x = sqrt(Ly*Ly - hb*hb);
		
		diag2_short = sqrt( (Lx-a_x)*(Lx-a_x) + hb*hb);
		diag2_long = sqrt( (Lx+a_x)*(Lx+a_x) + hb*hb);
		
		diag3_short = sqrt( diag2_short*diag2_short + Lz*Lz);
		diag3_long = sqrt( diag2_long*diag2_long + Lz*Lz);
		
		cut_off = diag3_long;
		cut_off_squared = diag3_long*diag3_long;
		

		V = Lx*hb*hz;
		
	
	
		
	}


   void parallelepiped::Set_Start_Lattice_Position(int id, double box_Lx, int N_box){
		
		//for cubic lattice
		
		int N_sitesp1, N_sitesp2, N_sitesp3;
		double N_sitesp1_float, N_sitesp2_float, N_sitesp3_float;
		double l_distpx, l_distpy, l_distpz;
		
		double V_q, a_q, b_q, c_q;
		double a,b;
		
		a_q=(Lx+a_x);
		b_q=hb;
		c_q=hz;
		
		V_q=a_q*b_q*c_q;
		
		
		a=b_q/a_q;
		b=c_q/a_q;
		
		
		N_sitesp1=floor(pow(double(N_box)/double(a*b), 1./3.));
		N_sitesp2=floor(N_sitesp1*a);
		N_sitesp3=ceil(double(N_box)/double(N_sitesp1*N_sitesp2));
	
	
		N_sitesp1_float = double(N_sitesp1);
		N_sitesp2_float = double(N_sitesp2);
		N_sitesp3_float = double(N_sitesp3);
		
		
		l_distpx = (box_Lx - N_sitesp1_float*a_q)/N_sitesp1_float;  
		l_distpy = (box_Lx - N_sitesp2_float*b_q)/N_sitesp2_float;
		l_distpz = (box_Lx - N_sitesp3_float*c_q)/N_sitesp3_float;
		 
		cout<<"l1 "<<l_distpx<<endl;
		cout<<"l2 "<<l_distpy<<endl;
		cout<<"l3 "<<l_distpz<<endl; 
		 
		
		x_center = Lx + l_distpx/2.0       + double(id%N_sitesp1)*(2.0*Lx + l_distpx );
		y_center = h_Lx + l_distpy/2.0     + double((id/N_sitesp1)%N_sitesp2)*(2.0*h_Lx + l_distpy );
		z_center = height_t + l_distpz/2.0 + double(id/int(N_sitesp1*N_sitesp2))*(2.0*height_t + l_distpz );
		
		
		edges_from_center();
		
		
		
	}

    void parallelepiped::Calculate_Face_Normals(){
		
		
		
		double face_ax;
		
		// face normals are the same as axis in the cube:
		
		
		
		facenormal[0].x = ax_1.y*ax_2.z - ax_1.z*ax_2.y;
		facenormal[0].y = ax_1.z*ax_2.x - ax_1.x*ax_2.z;
		facenormal[0].z = ax_1.x*ax_2.y - ax_1.y*ax_2.x;
		
		face_ax = facenormal[0].norm();		 
						
		facenormal[0].x = facenormal[0].x/face_ax;
		facenormal[0].y = facenormal[0].y/face_ax;
		facenormal[0].z = facenormal[0].z/face_ax;
					
	
		
		facenormal[1].x = ax_1.y*ax_3.z - ax_1.z*ax_3.y;
		facenormal[1].y = ax_1.z*ax_3.x - ax_1.x*ax_3.z;
		facenormal[1].z = ax_1.x*ax_3.y - ax_1.y*ax_3.x;
		
		face_ax = facenormal[1].norm();		 
						
		facenormal[1].x = facenormal[1].x/face_ax;
		facenormal[1].y = facenormal[1].y/face_ax;
		facenormal[1].z = facenormal[1].z/face_ax;
		
	
		facenormal[2].x = ax_2.y*ax_3.z - ax_2.z*ax_3.y;
		facenormal[2].y = ax_2.z*ax_3.x - ax_2.x*ax_3.z;
		facenormal[2].z = ax_2.x*ax_3.y - ax_2.y*ax_3.x;
		
		face_ax = facenormal[2].norm();		 
						
		facenormal[2].x = facenormal[2].x/face_ax;
		facenormal[2].y = facenormal[2].y/face_ax;
		facenormal[2].z = facenormal[2].z/face_ax;
		
		
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




	double parallelepiped::Calculate_Projection_to_Separating_Axis(m_vector laxis){
	
		
		double Rp;
		double rmax;
		double rmin;
		double scp_oc;
		
		double norm_ax;
	
	    distance_from_center();
	    
	    rmin=  dist_x[0]*laxis.x + dist_y[0]*laxis.y + dist_z[0]*laxis.z;
		rmax= rmin;
 
		for (int j=1;j<edge_N;j++){
			
			scp_oc = dist_x[j]*laxis.x + dist_y[j]*laxis.y + dist_z[j]*laxis.z;
			
			if (scp_oc < rmin) {
				rmin = scp_oc;
			} 
			else if (scp_oc > rmax) {
				rmax = scp_oc;
			}
			
			
		}	
		
	    Rp = fabs((rmax-rmin)/2.0);
		//Rp=rmax-rmin;	
		
		
		return Rp;
		
		
	}	
		
		
	
		
		
