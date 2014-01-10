


#include "octahedron.h"



	octahedron::octahedron(){
 
		 edge_N = 6
 
		 x = new double[edge_N];
		 y = new double[edge_N];
		 z = new double[edge_N];
		   
		 dist_x = new double[edge_N];
		 dist_y = new double[edge_N];
		 dist_z = new double[edge_N];
		 
		 new_dist_x = new double[edge_N];
		 new_dist_y = new double[edge_N];
		 new_dist_z = new double[edge_N];
		 
		 //trans_init = new double[3];
		 
		 edge_out = new int[edge_N];
		  
		 trans_old = new double[3];
		 Rot_old = new double[9];
		 
		 
		 
		
  
  
        }  


    octahedron::~ octahedron(){
  

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



    void octahedron::edges_from_center(){ //for octahedron aligned with coordinate ax
    
		
		 x[0] = x_center - 0.5*Lx;
		 y[0] = y_center - 0.5*Lx;
		 z[0] = z_center - 0.5*Lx;

		 x[1] = x_center + 0.5*Lx;
		 y[1] = y_center - 0.5*Lx;
		 z[1] = z_center - 0.5*Lx;

		 x[2] = x_center + 0.5*Lx;
		 y[2] = y_center + 0.5*Lx;
		 z[2] = z_center - 0.5*Lx;

		 x[3] = x_center - 0.5*Lx;
		 y[3] = y_center + 0.5*Lx;
		 z[3] = z_center - 0.5*Lx;
	
	
		 x[4] = x_center;	
		 y[4] = y_center;
		 z[4] = z_center + height;
		 
		 x[5] = x_center;	
		 y[5] = y_center;
		 z[5] = z_center - height;
		 
	
	}


	void octahedron::distance_from_center(){
 
	    for(int j=0;j<6;j++){
	   
		     dist_x[j] =  x[j] - x_center;
		     dist_y[j] =  y[j] - y_center;
		     dist_z[j] =  z[j] - z_center;
	  
	    }
    }


	void octahedron::Calculate_Axis(){
		
		int norm_ax;
		
		 ax_1.x =  x[1] - x[0];
		 ax_1.y =  y[1] - y[0];
	     ax_1.z =  z[1] - z[0];
	     
	     norm_ax=ax_1.norm();
				 
	     ax_1.x = ax_1.x/norm_ax;
	     ax_1.y = ax_1.y/norm_ax;
		 ax_1.z = ax_1.z/norm_ax;
				 
		 ax_2.x =  x[3] - x[0];
		 ax_2.y =  y[3] - y[0];
		 ax_2.z =  z[3] - z[0];
		
		norm_ax=ax_2.norm(); 
		 			
		 ax_2.x = ax_2.x/norm_ax;
		 ax_2.y = ax_2.y/norm_ax;
		 ax_2.z = ax_2.z/norm_ax;
				 
		 ax_3.x =  x[4] - x_center;
	     ax_3.y =  y[4] - y_center;
	     ax_3.z =  z[4] - z_center;
					
		norm_ax=ax_3.norm();
					
		 ax_3.x = ax_3.x/norm_ax;
		 ax_3.y = ax_3.y/norm_ax;
		 ax_3.z = ax_3.z/norm_ax;
		
	
		}


	void octahedron::Set_Axis(){

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

	
	void octahedron::Set_Lengths(){		
			
	    Lx = 1.0;
	    Ly = 1.0;
	    Lz = 1.0;
	   
	    height = (sqrt(3.0)/2.)*Lx;
	    
		cut_off = 2.0*height;
		cut_off_squared = cut_off*cutoff;
		V_p = (sqrt(2.0)/3.0)*Lx*Lx*Lx;
		
		

	}

	void octahedron::Set_Start_Lattice(int id, int packing_fraction, double box_Lx, double box_Ly, double box_Lz){
		
		//for cubic lattice
		
		int N_sitesp;
		double N_sitesp_float;
		double l_distp;
		
		
		N_sitesp1 = rint(box_Lx/Lx);
		N_sitesp2 = rint(box_Ly/Ly);
		N_sitesp3 = rint(box_Lz/Lz);
		
		N_sitesp1_float = double(N_sitesp1);  
		N_sitesp2_float = double(N_sitesp2);  
		N_sitesp3_float = double(N_sitesp3);
		
		x_center = Lx/2.0 + double(id%N_sitesp1)*Lx;
		y_center = Lx/2.0 + double((id/N_sitesp1)%N_sitesp2)*Lx;
		z_center = Lx/2.0 + double(id/int(N_sitesp1*N_sitesp2))*height;
		
		edges_from_center();
		
	}
	
	void octahedron::Calculate_Face_Normals(){
		
		// face normals are the same as axis in the cube:
		
		facenormal[0].x = ax_1.x;
		facenormal[0].y = ax_1.y;
		facenormal[0].z = ax_1.z;
		
		facenormal[1].x = ax_2.x;
		facenormal[1].y = ax_2.x;
		facenormal[1].z = ax_2.x;
		
		facenormal[2].x = ax_3.x;
		facenormal[2].y = ax_3.x;
		facenormal[2].z = ax_3.x;
		
	
		
    }		





