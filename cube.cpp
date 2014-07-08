
#include "cube.h"



	cube::cube(){
 
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


    cube::~ cube(){
  

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



    void cube::edges_from_center(){ // only valid if cubes are aligned with coordinate axes!
  
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

		 x[4] = x_center - 0.5*Lx;
		 y[4] = y_center - 0.5*Lx;
		 z[4] = z_center + 0.5*Lx;

		 x[5] = x_center + 0.5*Lx;
		 y[5] = y_center - 0.5*Lx;
		 z[5] = z_center + 0.5*Lx;

		 x[6] = x_center + 0.5*Lx;
		 y[6] = y_center + 0.5*Lx;
		 z[6] = z_center + 0.5*Lx;

		 x[7] = x_center - 0.5*Lx;
		 y[7] = y_center + 0.5*Lx;
		 z[7] = z_center + 0.5*Lx;

		//cout<<"Edges from center"<<endl;
		  

        }  


	void cube::distance_from_center(){
 
	    for(int j=0;j<edge_N;j++){
	   
		     dist_x[j] =  x[j] - x_center;
		     dist_y[j] =  y[j] - y_center;
		     dist_z[j] =  z[j] - z_center;
	  
	    }
	    
	    //cout<<"Distance from center"<<endl;
	    
    }
    
    


	void cube::Calculate_Axis(){
		
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
		
		/*
		cout<<"Calculate Axis"<<endl;
		cout<<"ax_1 "<<ax_1.x<<" "<<ax_1.y<<" "<<ax_1.z<<endl;
	    cout<<"ax_2 "<<ax_2.x<<" "<<ax_2.y<<" "<<ax_2.z<<endl;
		cout<<"ax_3 "<<ax_3.x<<" "<<ax_3.y<<" "<<ax_3.z<<endl;
		*/
		
	
	}

	
	void cube::Set_Axis(){

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



    void cube::Set_Lengths(){
		
		Lx = 1.0;
		cut_off = sqrt(3.0)*Lx;
		cut_off_squared = 3.0*Lx*Lx;
		
		V = Lx*Lx*Lx;
		
		//cout<<"Set lengths"<<endl;
		
	

	}

   void cube::Set_Start_Lattice_Position(int id, double box_Lx, int N_box){
		
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
		
		//cout<<"Set start lattice"<<endl;
		
		
		
	}

    void cube::Calculate_Face_Normals(){
		
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




	double cube::Calculate_Projection_to_Separating_Axis(m_vector laxis){
		
		double Rp;
	
		
		Rp = (fabs(ax_1.x*laxis.x + ax_1.y*laxis.y + ax_1.z*laxis.z) + 
			  fabs(ax_2.x*laxis.x + ax_2.y*laxis.y + ax_2.z*laxis.z) + 
			  fabs(ax_3.x*laxis.x + ax_3.y*laxis.y + ax_3.z*laxis.z) )*(Lx);
		
	
		
		return Rp;
		
		//cout<<"Calculate Projection to Separating Axis"<<endl;
		
		
	}	
		
		
		





