

#include "box.h"



box::box(){
  
	   packing_fraction = 0.0;
	   P_sigma = 0.0;
	   N = 0;
	   T = 0.0;
	   P = 0.0;
	   V = 0.0;
	   g_factor = 0.0; 
	   g_trans_vec = 0.0;
	   
	   Lx = 0.0;
	   Ly = 0.0;
	   Lz = 0.0; 
	   
	   
	   g_factor = 0;
	   g_trans_vec = 0;
		  
	   //V_p = 1.0;
		  
	   x_center = 0.0;
	   y_center = 0.0;
	   z_center = 0.0;
		  
	   edge_N = 8; 
	   
	   x=NULL;
	   y=NULL;
	   z=NULL;
	   
	   dist_x=NULL;
	   dist_y=NULL;
	   dist_z=NULL;
	   
	   new_dist_x=NULL;
	   new_dist_y=NULL;
	   new_dist_z=NULL;
	   
			  
		x = new double[edge_N];
	    y = new double[edge_N];
		z = new double[edge_N];
			   
		dist_x = new double[edge_N];
	    dist_y = new double[edge_N];
	    dist_z = new double[edge_N];
			 
		new_dist_x = new double[edge_N];
		new_dist_y = new double[edge_N];
		new_dist_z = new double[edge_N];

    }

void box::Startconfig(int N_in, double P_sigma_in, double packing_fraction_in, double V_p){

      packing_fraction = packing_fraction_in;
      P_sigma = P_sigma_in;
      N = N_in;
      //V_p = 1.0;
      
      T = 1.0;
      //P = P_sigma/(3.*sqrt(3.)/8.);
      P=P_sigma;
      //V = double(N)*V_p/packing_fraction;
      
	  V = double(N)*V_p/packing_fraction;
	  cout<<"V: "<<V<<" "<<"V_p:"<<V_p<<endl;	
     
      
      //for cubes
      Lx = pow(V,(1./3.));
      Ly = Lx;
      Lz = Lx;
      
      //for octahedron
      //double Vb;
      //Vb = V/sqrt(3.0);
      //Lx= pow(Vb, (1./3.));
      //Ly = Lx;
      //Lz = Lx*sqrt(3.);
	
	 
  
      g_factor = 100.0;
      g_trans_vec = g_factor*Lx;
      //V_p = 1.0;
      
      x_center = Lx/2.0;  
      y_center = Ly/2.0;
      z_center = Lz/2.0;
      
     
    }



void box::Startconfig_former(int N_in, double P_sigma_in){

     
      P_sigma = P_sigma_in;
      N = N_in;
	  //V_p = 1.0;
      T = 1.0;
      //P = P_sigma/(3.*sqrt(3.)/8.);
      P = P_sigma;
   
      g_factor = 100.0;
      g_trans_vec = g_factor*Lx;
      
     
    }



box::~box(){
	
packing_fraction = 0.0;
   P_sigma = 0.0;
   N = 0;
   T = 0.0;
   P = 0.0;
   V = 0.0;
   g_factor = 0.0; 
   g_trans_vec = 0.0;
   
   Lx = 0.0;
   Ly = 0.0;
   Lz = 0.0; 
   
	

}  


void box::edges_from_center(){ // only valid if cubes are aligned with coordinate axes!
  
  x[0] = x_center - 0.5*Lx;
  y[0] = y_center - 0.5*Ly;
  z[0] = z_center - 0.5*Lz;

  x[1] = x_center + 0.5*Lx;
  y[1] = y_center - 0.5*Ly;
  z[1] = z_center - 0.5*Lz;

  x[2] = x_center + 0.5*Lx;
  y[2] = y_center + 0.5*Ly;
  z[2] = z_center - 0.5*Lz;

  x[3] = x_center - 0.5*Lx;
  y[3] = y_center + 0.5*Ly;
  z[3] = z_center - 0.5*Lz;

  x[4] = x_center - 0.5*Lx;
  y[4] = y_center - 0.5*Ly;
  z[4] = z_center + 0.5*Lz;

  x[5] = x_center + 0.5*Lx;
  y[5] = y_center - 0.5*Ly;
  z[5] = z_center + 0.5*Lz;

  x[6] = x_center + 0.5*Lx;
  y[6] = y_center + 0.5*Ly;
  z[6] = z_center + 0.5*Lz;

  x[7] = x_center - 0.5*Lx;
  y[7] = y_center + 0.5*Ly;
  z[7] = z_center + 0.5*Lz;

  

}


void box::distance_from_center(){
 
 for(int j=0;j<8;j++){
   
 dist_x[j] =  x[j] - x_center;
 dist_y[j] =  y[j] - y_center;
 dist_z[j] =  z[j] - z_center;
  
 }
}


void box::Set_Lengths(double Vp){
	
	//for isotropic box
	
	V = double(N*Vp)/packing_fraction;
	Lx = pow(V,(1./3.));
    Ly = Lx;
    Lz = Lx;
       

}	

  

  
  
  
