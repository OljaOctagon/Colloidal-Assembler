


 void move::Vol_deform (particles& Particles, box* Box, fileio& Fileio, int mc_time){
				
				
	
		//SAVE OLD COORDINATES
		
		
		//Particles.Set_Cell_List(Box);	
		 
	 
		Box->V_old = Box->V;
		 
		Box->Lx_old = Box->Lx;
		Box->Ly_old = Box->Ly;
		Box->Lz_old = Box->Lz;
		
		Box->alpha_old = Box->alpha;
		Box->gamma_old = Box->gamma; 
		 
		
		//CHANGE BOX: ADD A SMALL PERTURBATION TO ONE OF THE COORDINATES OF ONE OF THE VECTORS
		
		
		rand_V=gsl_ran_flat (r01,0,9);
		rand_v=gsl_rng_uniform(r01);
		
		
		if (rand_V<3){
			if(rand_V<1){
				box_ax_1.x = box_ax1.x + dmax_v*(rand_v-0.5);
			}
			if(randV >=1 and rand_V<2){
				box_ax_1.y = box_ax1.y + dmax_v*(rand_v-0.5);
			}
			if(rand_V>=2){
				box_ax_1.z = box_ax1.z + dmax_v*(rand_v-0.5);
			}			
		 }
		 
		if (rand_V>=3 and rand_V <6)  {
			if(rand_V<4){
				box_ax_2.x = box_ax2.x + dmax_v*(rand_v-0.5);
			}
			if(randV >=4 and rand_V<5){
				box_ax_2.y = box_ax2.y + dmax_v*(rand_v-0.5);
			}
			if(rand_V>=5){
				box_ax_2.z = box_ax2.z + dmax_v*(rand_v-0.5);
			}			
			
		}
		
		if (rand_V>=6){	
			if(rand_V<7){
				box_ax_3.x = box_ax3.x + dmax_v*(rand_v-0.5);
			}
			if(randV >=7 and rand_V<8){
				box_ax_3.y = box_ax3.y + dmax_v*(rand_v-0.5);
			}
			if(rand_V>=8){
				box_ax_3.z = box_ax3.z + dmax_v*(rand_v-0.5);
			}			
		
		}
		
	    Box->Set_lengths();

		
		
		
		/*
		Box->V = Box->V_old + dmax_V*(rand_V-0.5);
		
		Box->V_rel = double(Box->V)/double(Box->V_old);
		Box->VL_rel = pow(Box->V_rel,(1./3.)); 
		 
		Box->Lx = Box->Lx*Box->VL_rel;
		Box->Ly = Box->Ly*Box->VL_rel;
		Box->Lz = Box->Lz*Box->VL_rel;
		  
		Box->Lx_scale = Box->Lx/Box->Lx_old;
		Box->Ly_scale = Box->Ly/Box->Ly_old;
		Box->Lz_scale = Box->Lz/Box->Lz_old;
		*/ 
		 
		  
        b_factor_pre = exp(-1.0*Box->P_sigma*(Box->V - Box->V_old) + Box->N*log(Box->V_rel));
		b_factor = minimum(1,b_factor_pre);
		  
		XI = gsl_rng_uniform(r01);
		  
	    
	    //b_factor<XI: Trial move not accepted
		if(b_factor<XI){
				
		     Box->V = Box->V_old;
					
		     Box->Lx = Box->Lx_old;
		     Box->Ly = Box->Ly_old;
			 Box->Lz = Box->Lz_old;
			 
			 Box->alpha = Box->alpha_old;
			 Box->beta = Box->beta_old;
			 
		
		}
			
			

//change box volume and check wheter excepted due to acceptance criterion

//change to fractional coordinates


//do trickts (periodic images, lattice reduction....)


// calculate possible overlaps


//accept or reject move acoording to exit of overlap detection


// to back to real coordinates


//The end
