#include "rhombohedron.h"



	rhombohedron::rhombohedron(){
 
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
		
		
         //cout<<"Cube Constructor"<<endl;		
		
  
        }  


    rhombohedron::~ rhombohedron(){
  

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



    void rhombohedron::edges_from_center(){ // only valid if cubes are aligned with coordinate axes!
  
		/*
		 x[0] = x_center + (-3*a_x + Lx/2.0);
		 y[0] = y_center  - h_2;
		 z[0] = z_center  - h_2;

		 x[1] = x_center + (-3*a_x + 3.0*Lx/2.0);
		 y[1] = y_center - h_2;
		 z[1] = z_center - h_2;

		 x[2] = x_center + (-2*a_x + 3.0*Lx/2.0);
		 y[2] = y_center + h_2;
		 z[2] = z_center - h_2;

		 x[3] = x_center + (-2*a_x + Lx/2.0);
		 y[3] = y_center + h_2;
		 z[3] = z_center - h_2;

		 x[4] = x_center + (-2*a_x + Lx/2.0);
		 y[4] = y_center - h_2;
		 z[4] = z_center + h_2;

		 x[5] = x_center + (-2*a_x + 3.0*Lx/2.0);
		 y[5] = y_center - h_2;
		 z[5] = z_center + h_2;

		 x[6] = x_center + (-a_x + 3.0*Lx/2.0);
		 y[6] = y_center + h_2;
		 z[6] = z_center + h_2;

		 x[7] = x_center - (-a_x + Lx/2.0);
		 y[7] = y_center + h_2;
		 z[7] = z_center + h_2;

		//cout<<"Edges from center"<<endl;
		*/  

		 x[0] = x_center + (-Lx -2*a_x)/2.0;
		 y[0] = y_center + (-h - a_x)/2.0;
		 z[0] = z_center + (-h)/2.0; 

		 x[1] = x_center + (Lx-2*a_x)/2.0;
		 y[1] = y_center + (-h - a_x)/2.0;
		 z[1] = z_center + (-h)/2.0;

		 x[2] = x_center + (Lx)/2.0;
		 y[2] = y_center + (h - a_x)/2.0;
		 z[2] = z_center + (-h)/2.0;

		 x[3] = x_center + (-Lx)/2.0;
		 y[3] = y_center + (h - a_x)/2.0;
		 z[3] = z_center + (-h)/2.0;

		 x[4] = x_center + (-Lx)/2.0;
		 y[4] = y_center + (-h + a_x)/2.0;
		 z[4] = z_center + (h)/2.0;

		 x[5] = x_center + (Lx)/2.0;
		 y[5] = y_center + (-h + a_x)/2.0;
		 z[5] = z_center + (h)/2.0;

		 x[6] = x_center + (Lx+2*a_x)/2.0;
		 y[6] = y_center + (h + a_x)/2.0;
		 z[6] = z_center + (h)/2.0;

		 x[7] = x_center + (-Lx+2*a_x)/2.0;
		 y[7] = y_center + (h + a_x)/2.0;
		 z[7] = z_center + (h)/2.0;



        }  


	void rhombohedron::distance_from_center(){
 
	    for(int j=0;j<edge_N;j++){
	   
		     dist_x[j] =  x[j] - x_center;
		     dist_y[j] =  y[j] - y_center;
		     dist_z[j] =  z[j] - z_center;
	  
	    }
	    
	    //cout<<"Distance from center"<<endl;
	    
    }
    
    


	void rhombohedron::Calculate_Axis(){
		
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

	
	void rhombohedron::Set_Axis(){

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



    void rhombohedron::Set_Lengths(){
		
		Lx = 1.0;
		alpha=85.0*M_PI/180.0;
		beta = M_PI - alpha;
		
		h=Lx*sin(alpha);
		
		h_2= double(h)/2.0;
		
		a_x = sqrt(Lx*Lx - h*h);
		
		diag2_short = sqrt( (Lx-a_x)*(Lx-a_x) + h*h);
		diag2_long = sqrt( (Lx+a_x)*(Lx+a_x) + h*h);
		
		diag3_short = sqrt( diag2_short*diag2_short + Lx*Lx);
		diag3_long = sqrt( diag2_long*diag2_long + Lx*Lx);
		
		cut_off = diag3_long;
		cut_off_squared = diag3_long*diag3_long;
		
	
		
		
		V = Lx*h*h;
		
	
		//cout<<"Set lengths"<<endl;
		
	

	}

	void rhombohedron::Set_Lengths(double Lx, double alpha){
		
		beta = M_PI - alpha;
		
		h=Lx*sin(alpha);
		
		h_2= double(h)/2.0;
		
		a_x = sqrt(Lx*Lx - h*h);
		
		diag2_short = sqrt( (Lx-a_x)*(Lx-a_x) + h*h);
		diag2_long = sqrt( (Lx+a_x)*(Lx+a_x) + h*h);
		
		diag3_short = sqrt( diag2_short*diag2_short + Lx*Lx);
		diag3_long = sqrt( diag2_long*diag2_long + Lx*Lx);
		
		cut_off = diag3_long;
		cut_off_squared = diag3_long*diag3_long;
		

		V = Lx*h*h;
		
	
			
		
	}


   void rhombohedron::Set_Start_Lattice_Position(int id, double box_Lx, int N_box){
		
		//for cubic lattice
		
		int N_sitesp_x;
		int N_sitesp_yz;
		
		double N_sitesp_x_float;
		double N_sitesp_yz_float;
		
		double l_distp_x;
		double l_distp_yz;
		
		
		N_sitesp_x = rint(pow(double(N_box*sin(alpha)*sin(alpha)),1./3.));
		N_sitesp_x_float = double(N_sitesp_x); 
		
	
		
		N_sitesp_yz = ceil (N_sitesp_x_float/sin(alpha));
		N_sitesp_yz_float = double(N_sitesp_yz);
	
		 
		  
		l_distp_x = (box_Lx - N_sitesp_x_float*Lx)/N_sitesp_x_float;   
		l_distp_yz = (box_Lx - N_sitesp_yz_float*h)/N_sitesp_yz_float;   
		
		
		/*
		x_center = p +  l_distp/2.0 + double(id %N_sitesp + l_distp*(id % N_sitesp));
		y_center = 0.5 +  l_distp/2.0 + double((id/N_sitesp)%N_sitesp + l_distp*((id/N_sitesp)%N_sitesp)); 
		z_center = 0.5 +  l_distp/2.0 + double(id/int(pow((double)N_sitesp,2)) + l_distp*(id/int((double)pow(N_sitesp,2))));
		*/
		
		
		x_center = a_x/2.0 +   l_distp_x/2.0  + double(id %N_sitesp_x)*Lx + double(id %N_sitesp_x)*l_distp_x;
		y_center = h_2 + l_distp_yz/2.0 + double((id/N_sitesp_x)%N_sitesp_yz)*h + double((id/N_sitesp_x)%N_sitesp_yz)*l_distp_yz;
		z_center = h_2 + l_distp_yz/2.0 + double(id/int(N_sitesp_x*N_sitesp_yz))*h + double((id/N_sitesp_x)%N_sitesp_yz)*l_distp_yz;
		 
		edges_from_center();
		
		//cout<<"Set start lattice"<<endl;
		
		
		
	}

    void rhombohedron::Calculate_Face_Normals(){
		
		
		
		double face_ax;
		
		// face normals are the same as axis in the cube:
		
		/*
		facenormal[0].x = ax_1.x;
		facenormal[0].y = ax_1.y;
		facenormal[0].z = ax_1.z;
		
		facenormal[1].x = ax_2.x;
		facenormal[1].y = ax_2.y;
		facenormal[1].z = ax_2.z;
		
		facenormal[2].x = ax_3.x;
		facenormal[2].y = ax_3.y;
		facenormal[2].z = ax_3.z;
		*/
		
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




	double rhombohedron::Calculate_Projection_to_Separating_Axis(m_vector laxis){
	
		
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
		
		
	
		
		
		
