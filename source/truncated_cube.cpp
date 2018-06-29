

#include "truncated_cube.h"



	truncated_cube::truncated_cube(){
 
		 
		 edge_N = 24;
		 N_independent_faces = 7;
		 N_cross_edges= 15;
		 edge_N_vec = 8;
		
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


    truncated_cube::~ truncated_cube(){
  

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

	void truncated_cube::write_positions(ofstream& fout){
				
	}	


    void truncated_cube::edges_from_center(){ //for truncated_cube aligned with coordinate ax
    
	
	
  
		lx = x_center - 0.5*Lx;
		ly = y_center - 0.5*Lx;
		lz = z_center - 0.5*Lx;
		
		x[0] = lx + edge_R*Lx; 
		y[0] = ly;
		z[0] = lz;
		
		x[1] = lx;
		y[1] = ly;
		z[1] = lz + edge_R*Lx;
		
		x[2] = lx;
		y[2] = ly + edge_R*Lx;
		z[2] = lz;
		
		lx = x_center + 0.5*Lx;
		ly = y_center - 0.5*Lx;
		lz = z_center - 0.5*Lx;
		
		x[3] = lx;
		y[3] = ly + edge_R*Lx;
		z[3] = lz;
		
		x[4] = lx;
		y[4] = ly;
		z[4] = lz + edge_R*Lx;
		
		x[5] = lx - edge_R*Lx;
		y[5] = ly;
		z[5] = lz;
		
		lx = x_center + 0.5*Lx;
		ly = y_center + 0.5*Lx;
		lz = z_center - 0.5*Lx;
		
		x[6] = lx - edge_R*Lx;
		y[6] = ly;
		z[6] = lz;
		
		x[7] = lx;
		y[7] = ly;
		z[7] = lz + edge_R*Lx;
		
		x[8] = lx;
		y[8] = ly - edge_R*Lx;
		z[8] = lz;
		
		lx = x_center - 0.5*Lx;
		ly = y_center + 0.5*Lx;
		lz = z_center - 0.5*Lx;
		
		x[9] = lx;
		y[9] = ly - edge_R*Lx;
		z[9] = lz;
		
		x[10] = lx;
		y[10] = ly;
		z[10] = lz + edge_R*Lx;
		
		x[11] = lx + edge_R*Lx;
		y[11] = ly;
		z[11] = lz;
		
		
		lx = x_center - 0.5*Lx;
		ly = y_center - 0.5*Lx;
		lz = z_center + 0.5*Lx;

		x[12] = lx + edge_R*Lx; 
		y[12] = ly;
		z[12] = lz;
		
		x[13] = lx;
		y[13] = ly + edge_R*Lx;
		z[13] = lz;
			
		x[14] = lx;
		y[14] = ly;
		z[14] = lz - edge_R*Lx;


		lx = x_center + 0.5*Lx;
		ly = y_center - 0.5*Lx;
		lz = z_center + 0.5*Lx;

		x[15] = lx;
		y[15] = ly + edge_R*Lx;
		z[15] = lz;
		
		x[16] = lx - edge_R*Lx;
		y[16] = ly;
		z[16] = lz;
			
		x[17] = lx;
		y[17] = ly;
		z[17] = lz - edge_R*Lx;

	
		lx = x_center + 0.5*Lx;
		ly = y_center + 0.5*Lx;
		lz = z_center + 0.5*Lx;

		
		x[18] = lx - edge_R*Lx;
		y[18] = ly;
		z[18] = lz;
		
		x[19] = lx;
		y[19] = ly - edge_R*Lx;
		z[19] = lz;

		x[20] = lx;
		y[20] = ly;
		z[20] = lz - edge_R*Lx;
		
		lx = x_center - 0.5*Lx;
		ly = y_center + 0.5*Lx;
		lz = z_center + 0.5*Lx;
		
		x[21] = lx + edge_R*Lx;
		y[21] = ly;
		z[21] = lz;
		
		x[22] = lx;
		y[22] = ly - edge_R*Lx;
		z[22] = lz;
		
		x[23] = lx;
		y[23] = ly;
		z[23] = lz - edge_R*Lx;
		
		  

	}  

	

	void truncated_cube::distance_from_center(){
 
	    for(int j=0;j<edge_N;j++){
	   
		     dist_x[j] =  x[j] - x_center;
		     dist_y[j] =  y[j] - y_center;
		     dist_z[j] =  z[j] - z_center;
	  
	    }
    }


	void truncated_cube::Calculate_Axis(){
		
		double norm_ax;
		
		 ax_1.x =  x[5] - x[0];
		 ax_1.y =  y[5] - y[0];
	     ax_1.z =  z[5] - z[0];
	     
	     norm_ax=ax_1.norm();
				 
	     ax_1.x = ax_1.x/norm_ax;
	     ax_1.y = ax_1.y/norm_ax;
		 ax_1.z = ax_1.z/norm_ax;
				 
		 ax_2.x =  x[9] - x[2];
		 ax_2.y =  y[9] - y[2];
		 ax_2.z =  z[9] - z[2];
		
		norm_ax=ax_2.norm(); 
		 			
		 ax_2.x = ax_2.x/norm_ax;
		 ax_2.y = ax_2.y/norm_ax;
		 ax_2.z = ax_2.z/norm_ax;
				 
		 ax_3.x =  x[14] - x[1];
	     ax_3.y =  y[14] - y[1];
	     ax_3.z =  z[14] - z[1];
					
		norm_ax=ax_3.norm();
					
		 ax_3.x = ax_3.x/norm_ax;
		 ax_3.y = ax_3.y/norm_ax;
		 ax_3.z = ax_3.z/norm_ax;
		
	
	}


	void truncated_cube::Set_Axis(){

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

	
	void truncated_cube::Set_Lengths(){		
		
		Lx = 1.0;
		edge_R=0.3;
	    
		cut_off = Lx*sqrt(3.0);
		
		cut_off_squared = cut_off*cut_off;
		
		V= Lx*Lx*Lx*(1.0 - (4./3.)*edge_R*edge_R*edge_R);
		
		
		
	}


	
	void truncated_cube::Set_Start_Lattice_Position(int id, double box_Lx, int N_box){
		
		//for cubic lattice
		
		int N_sitesp;
		double N_sitesp_float;
		double l_distp;
		
		
		N_sitesp = rint(pow(double(N_box),1./3.));
		N_sitesp_float = double(N_sitesp);  
		  
		l_distp = (box_Lx - N_sitesp_float)/N_sitesp_float;   
		
		x_center = 0.5 +  l_distp/2.0 + double(id %N_sitesp + l_distp*(id % N_sitesp));
		y_center = 0.5 +  l_distp/2.0 + double((id/N_sitesp)%N_sitesp + l_distp*((id/N_sitesp)%N_sitesp)); 
		z_center = 0.5 +  l_distp/2.0 + double(id/int(pow((double)N_sitesp,2)) + l_distp*(id/int((double)pow(N_sitesp,2))));
		  	 
		edges_from_center();
		
		
		
	}
	
	

	void truncated_cube::Calculate_Face_Normals(){
		
		
		//calculating facenormals for the four independent triangles
		
		
		//calculate edges for independant triangle # 1
		
		edges[0].x = x[1] - x[0];
		edges[0].y = y[1] - y[0];
		edges[0].z = z[1] - z[0];
		
		edges[1].x = x[2] - x[0];
		edges[1].y = y[2] - y[0];
		edges[1].z = z[2] - z[0];
		
		//calculate edges for independant triangle # 2
		
		edges[2].x = x[3] - x[5];
		edges[2].y = y[3] - y[5];
		edges[2].z = z[3] - z[5];
		
		edges[3].x = x[4] - x[5];
		edges[3].y = y[4] - y[5];
		edges[3].z = z[4] - z[5];
		
	    //calculate edges for independant triangle # 3
		
		edges[4].x = x[12] - x[14];
		edges[4].y = y[12] - y[14];
		edges[4].z = z[12] - z[14];
		
		edges[5].x = x[13] - x[14];
		edges[5].y = y[13] - y[14];
		edges[5].z = z[13] - z[14];
		
		
		//calculate edges for independant triangle # 4
		
		edges[6].x = x[15] - x[17];
		edges[6].y = y[15] - y[17];
		edges[6].z = z[15] - z[17];
		
		edges[7].x = x[16] - x[17];
		edges[7].y = y[16] - y[17];
		edges[7].z = z[16] - z[17];
		
		
		//for octagons facenormals are axes
		
		facenormal[0].x = ax_1.x;
		facenormal[0].y = ax_1.y;
		facenormal[0].z = ax_1.z;
		
		facenormal[1].x = ax_2.x;
		facenormal[1].y = ax_2.y;
		facenormal[1].z = ax_2.z;
		
		facenormal[2].x = ax_3.x;
		facenormal[2].y = ax_3.y;
		facenormal[2].z = ax_3.z;
		
		
		// triangle 1 -> facnormal # 4 
		
		facenormal[3].x = edges[1].y*edges[0].z - edges[1].z*edges[0].y; 
		facenormal[3].y = edges[1].z*edges[0].x - edges[1].x*edges[0].z;
		facenormal[3].z = edges[1].x*edges[0].y - edges[1].y*edges[0].x;
		
		
		// triangle 2 -> facenormal # 5
		
		facenormal[4].x = edges[3].y*edges[2].z - edges[3].z*edges[2].y; 
		facenormal[4].y = edges[3].z*edges[2].x - edges[3].x*edges[2].z;
		facenormal[4].z = edges[3].x*edges[2].y - edges[3].y*edges[2].x;
		
		
		// triangle 3 -> facenormal # 6
		
		facenormal[5].x = edges[5].y*edges[4].z - edges[5].z*edges[4].y; 
		facenormal[5].y = edges[5].z*edges[4].x - edges[5].x*edges[4].z;
		facenormal[5].z = edges[5].x*edges[4].y - edges[5].y*edges[4].x;
		
		
		// triangle 4 -> facenormal # 7
		
		facenormal[6].x = edges[7].y*edges[6].z - edges[7].z*edges[6].y; 
		facenormal[6].y = edges[7].z*edges[6].x - edges[7].x*edges[6].z;
		facenormal[6].z = edges[7].x*edges[6].y - edges[7].y*edges[6].x;
		
		
	
		double f_norm;
		
	    for(int k=0;k<N_independent_faces; k++){
		
			f_norm = facenormal[k].norm();
			
			facenormal[k].x = facenormal[k].x/f_norm;
			facenormal[k].y = facenormal[k].y/f_norm;
			facenormal[k].z = facenormal[k].z/f_norm;
			
		}
		
		
    }		




	double truncated_cube::Calculate_Projection_to_Separating_Axis(m_vector laxis){
		
		
		double scp_oc;
	    
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
		
		
		return Rp;
		
		
	}	
		
		
	




