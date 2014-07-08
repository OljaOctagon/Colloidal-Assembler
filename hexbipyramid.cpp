

#include "hexbipyramid.h"



	hexbipyramid::hexbipyramid(){
 
		 edge_N = 18;
		 N_independent_faces = 7;
		 N_cross_edges=9;
		 edge_N_vec = 14;
		
		 x = new double[edge_N];
		 y = new double[edge_N];
		 z = new double[edge_N];
		   
		 dist_x = new double[edge_N];
		 dist_y = new double[edge_N];
		 dist_z = new double[edge_N];
		 
		 new_dist_x = new double[edge_N];
		 new_dist_y = new double[edge_N];
		 new_dist_z = new double[edge_N];
		
		 facenormal = new m_vector[N_independent_faces];	
		 edges = new m_vector[edge_N_vec];
		 
		 //trans_init = new double[3];
		 
		 edge_out = new int[edge_N];
		 
		 trans_periodic = new double[3]; 
		 trans_old = new double[3];
		 Rot_old = new double[9];
		
		
		 copy_count = 0;
		
	
  
        }  


    hexbipyramid::~ hexbipyramid(){
  

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
		 
		 
		 delete[] trans_old;
		 delete[] Rot_old;
		
  
        }  



    void hexbipyramid::edges_from_center(){ //for hexbipyramid aligned with coordinate ax
    
		
		 x[0] = x_center + Lx;
		 y[0] = y_center;
		 z[0] = z_center; 
		
		 x[1] = x_center + Lx/2.0;
		 y[1] = y_center + h_Lx;	
		 z[1] = z_center;	

		 x[2] = x_center - Lx/2.0;
		 y[2] = y_center + h_Lx;
		 z[2] = z_center;
		
		 x[3] =  x_center - Lx;
		 y[3] =  y_center;
		 z[3] =  z_center;
	
		 x[4] = x_center - Lx/2.0;
		 y[4] = y_center - h_Lx;
		 z[4] = z_center;
		 
		 x[5] = x_center + Lx/2.0;
		 y[5] = y_center - h_Lx;
		 z[5] = z_center;
		
		 x[6] = x_center  + Lx_t;
		 y[6] = y_center;
		 z[6] = z_center + height_t;

         x[7] = x_center + Lx_t/2.0;
		 y[7] = y_center + h_Lx_t;
		 z[7] = z_center + height_t;
    
		 x[8] = x_center - Lx_t/2.0;
		 y[8] = y_center + h_Lx_t;
		 z[8] = z_center + height_t;
		 
		 x[9] = x_center - Lx_t;
		 y[9] = y_center;
		 z[9] = z_center + height_t;

		 x[10] = x_center - Lx_t/2.0;
		 y[10] = y_center - h_Lx_t;
		 z[10] = z_center + height_t;
		 
		 x[11] = x_center + Lx_t/2.0;
		 y[11] = y_center - h_Lx_t;
		 z[11] = z_center + height_t;
		
		 x[12] = x_center + Lx_t;
		 y[12] = y_center; 
		 z[12] = z_center - height_t;

         x[13] = x_center + Lx_t/2.0;
		 y[13] = y_center + h_Lx_t;
		 z[13] = z_center - height_t;

		 x[14] = x_center - Lx_t/2.0;
		 y[14] = y_center + h_Lx_t;
		 z[14] = z_center - height_t;

		 x[15] = x_center - Lx_t;
		 y[15] = y_center;
		 z[15] = z_center - height_t;


		 x[16] = x_center - Lx_t/2.0;
		 y[16] = y_center - h_Lx_t;
		 z[16] = z_center - height_t;
		 
		 x[17] = x_center + Lx_t/2.0;
		 y[17] = y_center - h_Lx_t;
		 z[17] = z_center - height_t;
		
	
		 
	
	}


	void hexbipyramid::distance_from_center(){
 
	    for(int j=0;j<edge_N;j++){
	   
		     dist_x[j] =  x[j] - x_center;
		     dist_y[j] =  y[j] - y_center;
		     dist_z[j] =  z[j] - z_center;
	  
	    }
    }


	void hexbipyramid::Calculate_Axis(){
		
		int norm_ax;
		
		 ax_1.x =  x[1] - x[0];
		 ax_1.y =  y[1] - y[0];
	     ax_1.z =  z[1] - z[0];
	     
	     norm_ax=ax_1.norm();
				 
	     ax_1.x = ax_1.x/norm_ax;
	     ax_1.y = ax_1.y/norm_ax;
		 ax_1.z = ax_1.z/norm_ax;
				 
		 ax_2.x =  x[5] - x[0];
		 ax_2.y =  y[5] - y[0];
		 ax_2.z =  z[5] - z[0];
		
		norm_ax=ax_2.norm(); 
		 			
		 ax_2.x = ax_2.x/norm_ax;
		 ax_2.y = ax_2.y/norm_ax;
		 ax_2.z = ax_2.z/norm_ax;
				 
		 ax_3.x =  x[6] - x[0];
	     ax_3.y =  y[6] - y[0];
	     ax_3.z =  z[6] - z[0];
					
		norm_ax=ax_3.norm();
					
		 ax_3.x = ax_3.x/norm_ax;
		 ax_3.y = ax_3.y/norm_ax;
		 ax_3.z = ax_3.z/norm_ax;
		
	
	}


	void hexbipyramid::Set_Axis(){

	     ax_1_old.x =  ax_1.x;
		 ax_1_old.y =  ax_1.y;
		 ax_1_old.z =  ax_1.z;
			 	 
		 ax_2_old.x =  ax_2.x;
		 ax_2_old.y =  ax_2.y;
		 ax_2_old.z =  ax_2.z;
			
	 	 ax_3_old.x =  ax_3.x;
		 ax_3_old.y =  ax_3.y;
		 ax_3_old.z =  ax_3.z;
			
			
    } 

	
	void hexbipyramid::Set_Lengths(){		
			
	   
		Lx = 0.5;
		height = 1.0;
	    height_t=0.4*height;
	    height_o=0.6*height;
	    
	    Lx_t=0.6*Lx; 
	    
	 	cut_off = sqrt(Lx_t*Lx_t + height_t*height_t)*2.0;
	 	cout <<"cut_off"<<cut_off<<endl;
		
		cut_off_squared = cut_off*cut_off;
		//cf= cos(M_PI/6.0)/sin(M_PI/6.0);
		
		Vt = height_o*Lx_t*Lx_t*sqrt(3.0);
		V= height*Lx*Lx*sqrt(3.0);
		V= V-Vt;
		
		
		h_Lx= (sqrt(3)/2.0)*Lx;
		h_Lx_t= (sqrt(3)/2.0)*Lx_t;
		
		
	}


	
	void hexbipyramid::Set_Start_Lattice_Position(int id, double box_Lx, int N_box){
		
		//for cubic lattice
		
		int N_sitesp1, N_sitesp2, N_sitesp3;
		double N_sitesp1_float, N_sitesp2_float, N_sitesp3_float;
		double l_distpx, l_distpy, l_distpz;
		
		double V_q, a_q, b_q, c_q;
		double phi_q;
		double a,b;
		
		a_q=2*Lx;
		b_q=2*h_Lx;
		c_q=2*height_t;
		
		V_q=a_q*b_q*c_q;
		
		phi_q=double(N_box*V_q)/double(box_Lx*box_Lx*box_Lx);
		
		a=b_q/a_q;
		b=c_q/a_q;
		
		
		
		//N_sitesp1 = floor(box_Lx/(2.0*Lx));
		//N_sitesp2 = floor(box_Lx/(2.0*h_Lx));
		//N_sitesp3 = floor(box_Lx/(2.0*height_t));
		
		//N_sitesp1=ceil(phi_q*N_sitesp1);
		//N_sitesp2=ceil(phi_q*N_sitesp2);
		//N_sitesp3=ceil(phi_q*N_sitesp3);
		
		N_sitesp1=floor(pow(double(N_box)/double(a*b), 1./3.));
		N_sitesp2=floor(N_sitesp1*a);
		N_sitesp3=ceil(double(N_box)/double(N_sitesp1*N_sitesp2));
		
		
		cout<<"N_sitesp1"<<N_sitesp1<<endl;
		cout<<"N_sitesp2"<<N_sitesp2<<endl;
		cout<<"N_sitesp3"<<N_sitesp3<<endl;
		
	
		
		N_sitesp1_float = double(N_sitesp1);
		N_sitesp2_float = double(N_sitesp2);
		N_sitesp3_float = double(N_sitesp3);
		
		
		
		l_distpx = (box_Lx - N_sitesp1_float*2.0*Lx)/N_sitesp1_float;  
		l_distpy = (box_Lx - N_sitesp2_float*2.0*h_Lx)/N_sitesp2_float;
		l_distpz = (box_Lx - N_sitesp3_float*2.0*height_t)/N_sitesp3_float;
		 
		cout<<"l1 "<<l_distpx<<endl;
		cout<<"l2 "<<l_distpy<<endl;
		cout<<"l3 "<<l_distpz<<endl; 
		 
		
		x_center = Lx + l_distpx/2.0       + double(id%N_sitesp1)*(2.0*Lx + l_distpx );
		y_center = h_Lx + l_distpy/2.0     + double((id/N_sitesp1)%N_sitesp2)*(2.0*h_Lx + l_distpy );
		z_center = height_t + l_distpz/2.0 + double(id/int(N_sitesp1*N_sitesp2))*(2.0*height_t + l_distpz );
		
		
		edges_from_center();
		
		
	}
	
	

	void hexbipyramid::Calculate_Face_Normals(){
		
		
		edges[0].x = x[1] - x[0];
		edges[0].y = y[1] - y[0];
		edges[0].z = z[1] - z[0];
	
		edges[1].x = x[6] - x[0];
		edges[1].y = y[6] - y[0];
		edges[1].z = z[6] - z[0];
		
		edges[2].x = x[2] - x[1];
		edges[2].y = y[2] - y[1];
		edges[2].z = z[2] - z[1];
		
		edges[3].x = x[7] - x[1];
		edges[3].y = y[7] - y[1];
		edges[3].z = z[7] - z[1];
		
		edges[4].x = x[3] - x[2];
		edges[4].y = y[3] - y[2];
		edges[4].z = z[3] - z[2];
		
		edges[5].x = x[8] - x[2];
		edges[5].y = y[8] - y[2];
		edges[5].z = z[8] - z[2];
				
		edges[6].x = x[9] - x[3];
		edges[6].y = y[9] - y[3];
		edges[6].z = z[9] - z[3];
	
	    edges[7].x = x[11] - x[5];
		edges[7].y = y[11] - y[5];
		edges[7].z = z[11] - z[5];
		
		edges[8].x = x[10] - x[4];
		edges[8].y = y[10] - y[4];
		edges[8].z = z[10] - z[4];
				
		edges[9].x = x[5] - x[4];
		edges[9].y = y[5] - y[4];
		edges[9].z = z[5] - z[4];
	
		edges[10].x = x[4] - x[3];
		edges[10].y = y[4] - y[3];
		edges[10].z = z[4] - z[3];
	
		edges[11].x = x[0] - x[5];
		edges[11].y = y[0] - y[5];
		edges[11].z = z[0] - z[5];
			
		edges[12].x = x[7] - x[6];
		edges[12].y = y[7] - y[6];
		edges[12].z = z[7] - z[6];
		
		edges[13].x = x[11] - x[6];
		edges[13].y = y[11] - y[6];
		edges[13].z = z[11] - z[6];
		
		
		//calculate facenormal with cross products of edges
		
		facenormal[0].x = edges[0].y*edges[1].z - edges[0].z*edges[1].y; 
		facenormal[0].y = edges[0].z*edges[1].x - edges[0].x*edges[1].z;
		facenormal[0].z = edges[0].x*edges[1].y - edges[0].y*edges[1].x;
		
		facenormal[1].x = edges[2].y*edges[3].z - edges[2].z*edges[3].y; 
		facenormal[1].y = edges[2].z*edges[3].x - edges[2].x*edges[3].z;
		facenormal[1].z = edges[2].x*edges[3].y - edges[2].y*edges[3].x;
		
		facenormal[2].x = edges[4].y*edges[5].z - edges[4].z*edges[5].y; 
		facenormal[2].y = edges[4].z*edges[5].x - edges[4].x*edges[5].z;
		facenormal[2].z = edges[4].x*edges[5].y - edges[4].y*edges[5].x;
		
		facenormal[3].x = edges[10].y*edges[6].z - edges[10].z*edges[6].y; 
		facenormal[3].y = edges[10].z*edges[6].x - edges[10].x*edges[6].z;
		facenormal[3].z = edges[10].x*edges[6].y - edges[10].y*edges[6].x;
		
		facenormal[4].x = edges[9].y*edges[8].z - edges[9].z*edges[8].y; 
		facenormal[4].y = edges[9].z*edges[8].x - edges[9].x*edges[8].z;
		facenormal[4].z = edges[9].x*edges[8].y - edges[9].y*edges[8].x;
		
		facenormal[5].x = edges[11].y*edges[7].z - edges[11].z*edges[7].y; 
		facenormal[5].y = edges[11].z*edges[7].x - edges[11].x*edges[7].z;
		facenormal[5].z = edges[11].x*edges[7].y - edges[11].y*edges[7].x;
		
		facenormal[6].x = edges[12].y*edges[13].z - edges[12].z*edges[13].y; 
		facenormal[6].y = edges[12].z*edges[13].x - edges[12].x*edges[13].z;
		facenormal[6].z = edges[12].x*edges[13].y - edges[12].y*edges[13].x;
		
	
		double f_norm;
		
	    for(int k=0;k<N_independent_faces; k++){
		
			f_norm = facenormal[k].norm();
			
			facenormal[k].x = facenormal[k].x/f_norm;
			facenormal[k].y = facenormal[k].y/f_norm;
			facenormal[k].z = facenormal[k].z/f_norm;
			
		}
		
		
    }		




	double hexbipyramid::Calculate_Projection_to_Separating_Axis(m_vector laxis){
		
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
		//Rp=fabs(rmax-rmin);	
		
		//Rp = (fabs(ax_1.x*laxis.x + ax_1.y*laxis.y + ax_1.z*laxis.z) + fabs(ax_2.x*laxis.x + ax_2.y*laxis.y + ax_2.z*laxis.z))*(Lx/2.0) 
		//					+ fabs(ax_3.x*laxis.x + ax_3.y*laxis.y + ax_3.z*laxis.z)*(height/2.0);
		
		
		
		
		return Rp;
		
		
	}	
		
