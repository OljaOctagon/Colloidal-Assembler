
#include "order_parameter.h"

void order_parameter::g_radial(double r_start, double r_stop, double delta_r){
	
	
		for(int id=0;id<Box->N;id++){
			
			
				for(int j=0;j<Collision_List[id].Nm;j++){
						
					num_r1 = r_start;	
						
					for(int p=0;p<N_histo_points;p++){
					 
						num_r2 = r_start + delta_r*(p+1);
						
						r_dist = Collision_List[id].Elements[j].distance_norm;
						
						if((r_dist>=num_r1)&&(r_dist<num_r2)){
							
						   g_radial_histo[p-1] = g_radial_histo[p-1] + 1;	
							
						}
					 
						num_r1 = num_r2;
					 
					}
					
			    
			    }
	
	    } 

		rn = r_start;	

		for(int p = 0;p<N_histo_points;p++){
		
			rn = r_start + delta_r*p;
			g_norm = 4.0*M_PI*rn*rn*double(Box->N)*(double(Box->N)/Box->V);
			
			g_radial_distribution[p].x = rn;
			g_radial_distribution[p].y = double(g_radial_histo[p])/g_norm;
			
		}	
	
		g_min = 10.0;
	
		for(int p=0;p<N_histo_points;p++){
			
			if(g_radial_distribution[p].y<g_min){
				
				p_min = p;
				g_min = g_radial_distribution[p].y;			
		    }
		    
					
		}
	
	
        g_Cut_Off = g_radial_distribution[p_min].x;
	
	
	
	
    }
