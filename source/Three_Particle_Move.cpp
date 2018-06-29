
#include "move.h"
#include <gsl/gsl_math.h>


void move::Three_Particle_Move(particles& Particles, box* Box, fileio& Fileio, int mc_time){


	double phi_t;
	double dU_old, dU_new;	
	double delta_U;
	dU_old=0;
	dU_new=0;
	double n1,n2;
	int list_id1, list_id2;
	

	id = gsl_rng_uniform_int(r, Box->N);
	
 
 	//Check with whom the particle is interacting


	Particles.Collision_List[id].Calculate_Bonds(Box, id, Particles.N_Particle, Particles.MAX_coll_p);
	
	if (Particles.Collision_List[id].Nm >=2){


		n1 = Particles.Collision_List[id].Nm;

		list_id1 = gsl_rng_uniform_int(r, Particles.Collision_List[id].Nm);

		do{
			list_id2 = gsl_rng_uniform_int(r, Particles.Collision_List[id].Nm);

		}while(list_id2==list_id1);	

		
		list_id1 = Particles.Collision_List[id].Elements[list_id1].nl_id;
		list_id2 = Particles.Collision_List[id].Elements[list_id2].nl_id;


		Particles.Collision_List[id].Calculate_OP(Box, id, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);
		dU_old = Calculate_Pair_Potential(id, Particles, Box, Particles.Collision_List);
	
		Particles.Collision_List[list_id1].Calculate_OP(Box, list_id1, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);		
		dU_old = dU_old + Calculate_Pair_Potential(list_id1, Particles, Box, Particles.Collision_List);	

		Particles.Collision_List[list_id2].Calculate_OP(Box, list_id2, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);		
		dU_old = dU_old + Calculate_Pair_Potential(list_id2, Particles, Box, Particles.Collision_List);	




		// rotate id particle 

		/*
		double sign1, sign2;
		sign1=1;
		sign2=1;
		double scp1, scp2;

		
		Particles.N_Particle[id]->Calculate_Long_Axis();
		Particles.N_Particle[list_id1]->Calculate_Long_Axis();
		Particles.N_Particle[list_id2]->Calculate_Long_Axis();

		scp1 = Particles.N_Particle[id]->long_axis.x*Particles.N_Particle[list_id1]->long_axis.x +
		 	   Particles.N_Particle[id]->long_axis.y*Particles.N_Particle[list_id1]->long_axis.y +
		       Particles.N_Particle[id]->long_axis.z*Particles.N_Particle[list_id1]->long_axis.z;

		scp2 = Particles.N_Particle[id]->long_axis.x*Particles.N_Particle[list_id2]->long_axis.x +
		 	   Particles.N_Particle[id]->long_axis.y*Particles.N_Particle[list_id2]->long_axis.y +
		       Particles.N_Particle[id]->long_axis.z*Particles.N_Particle[list_id2]->long_axis.z;       


		
		if (scp1<0){
			sign1=-1.0;
		}

		if (scp2<0){
			sign2=-1.0;
		}
		*/


		Particles.N_Particle[id]->phi = Particles.N_Particle[id]->phi + ((60.0*M_PI)/180.0);	
		

		if (Particles.N_Particle[id]->phi > 2*M_PI){
			Particles.N_Particle[id]->phi = Particles.N_Particle[id]->phi - 2.0*M_PI;
		}
			
		if (Particles.N_Particle[id]->phi < 0){	
			Particles.N_Particle[id]->phi = Particles.N_Particle[id]->phi + 2.0*M_PI;
		}
		

		phi_t = Particles.N_Particle[id]->phi;

		Rot_mat[0] = cos(phi_t);
		Rot_mat[1] = -sin(phi_t);
		Rot_mat[2] = 0;
			
		Rot_mat[3] = sin(phi_t);
		Rot_mat[4] = cos(phi_t);
		Rot_mat[5] = 0;
			
		Rot_mat[6] = 0;
		Rot_mat[7] = 0;
		Rot_mat[8] = 1;
		

		Rot_Update_Positions(Particles, id, Rot_mat);
		Particles.N_Particle[id]->Calculate_Axis();
		Particles.N_Particle[id]->Calculate_Patch_Position();

		// rotate list_id1 particle
	
		Particles.N_Particle[list_id1]->phi = Particles.N_Particle[list_id1]->phi + (60.0*M_PI)/180.0;


		if (Particles.N_Particle[list_id1]->phi > 2*M_PI){
			Particles.N_Particle[list_id1]->phi = Particles.N_Particle[list_id1]->phi - 2.0*M_PI;
		}
			
		if (Particles.N_Particle[list_id1]->phi < 0){	
			Particles.N_Particle[list_id1]->phi = Particles.N_Particle[list_id1]->phi + 2.0*M_PI;
		}


		Particles.N_Particle[list_id2]->phi = Particles.N_Particle[list_id2]->phi + (60.0*M_PI)/180.0;	

		
		if (Particles.N_Particle[list_id2]->phi > 2*M_PI){
			Particles.N_Particle[list_id2]->phi = Particles.N_Particle[list_id2]->phi - 2.0*M_PI;
		}
			
		if (Particles.N_Particle[list_id2]->phi < 0){	
			Particles.N_Particle[list_id2]->phi = Particles.N_Particle[list_id2]->phi + 2.0*M_PI;
		}

		
		phi_t = (60.0*M_PI)/180.0;

		Rot_mat[0] = cos(phi_t);
		Rot_mat[1] = -sin(phi_t);
		Rot_mat[2] = 0;
			
		Rot_mat[3] = sin(phi_t);
		Rot_mat[4] = cos(phi_t);
		Rot_mat[5] = 0;
			
		Rot_mat[6] = 0;
		Rot_mat[7] = 0;
		Rot_mat[8] = 1;


		Rot_Move_Map(Particles, id, list_id1, Box, Rot_mat);
		Particles.Check_Periodic_CM(list_id1, Box);		
		Update_Periodic_Positions(Particles, Box, list_id1);	

		Particles.N_Particle[list_id1]->Calculate_Axis();
		Particles.N_Particle[list_id1]->Calculate_Patch_Position();

		// rotate particle list_id2

		
		Rot_Move_Map(Particles, id, list_id2, Box, Rot_mat);

		Particles.Check_Periodic_CM(list_id2, Box);		
		Update_Periodic_Positions(Particles, Box, list_id2);	

		Particles.N_Particle[list_id2]->Calculate_Axis();
		Particles.N_Particle[list_id2]->Calculate_Patch_Position();


	
		// Calculate Collisions

	    exit_status = 0;
		col_count = 0;
		

	    Particles.Collision_List[id].Calculate_OP(Box, id, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);
		Collision_Test( Particles, Box, id, Particles.Collision_List);
		
		 
		Particles.Collision_List[list_id1].Calculate_OP(Box, list_id1, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);
		Collision_Test( Particles, Box, list_id1, Particles.Collision_List);
	
		
		Particles.Collision_List[list_id2].Calculate_OP(Box, list_id2, Particles.N_Particle, Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);
		Collision_Test( Particles, Box, list_id2, Particles.Collision_List);		
	

		
		if(exit_status > 0){

	
				
			Reset_Positions(Particles, id);	
			Particles.N_Particle[id]->Calculate_Axis();
			Particles.N_Particle[id]->Calculate_Patch_Position();

			Reset_Positions(Particles, list_id1);	
			Particles.N_Particle[list_id1]->Calculate_Axis();
			Particles.N_Particle[list_id1]->Calculate_Patch_Position();

			Reset_Positions(Particles, list_id2);	
			Particles.N_Particle[list_id2]->Calculate_Axis();
			Particles.N_Particle[list_id2]->Calculate_Patch_Position();

	
			
		}	
					
						
		if (exit_status == 0){
		
			
			

			dU_new =  Calculate_Pair_Potential(id, Particles, Box, Particles.Collision_List);
		    dU_new = dU_new +  Calculate_Pair_Potential(list_id1, Particles, Box, Particles.Collision_List);
		    dU_new = dU_new +  Calculate_Pair_Potential(list_id2, Particles, Box, Particles.Collision_List);
	 		
		    

			//Calculate_Pair_Potential( Particles, Box);		


			n2 = Particles.Collision_List[id].Nm;


			delta_U = dU_new - dU_old;


			b_factor_pre = ((n1*(n1-1))/(n2*(n2-1)))*exp(-1.0*beta*delta_U);
			
			b_factor = minimum(1,b_factor_pre);
			  
			XI = gsl_rng_uniform(r01);
			  	
			
			// Reject of XI > b_factor

			if (XI > b_factor){

			 	Reset_Positions(Particles, id);	
				Particles.N_Particle[id]->Calculate_Axis();
				Particles.N_Particle[id]->Calculate_Patch_Position();

				Reset_Positions(Particles, list_id1);	
				Particles.N_Particle[list_id1]->Calculate_Axis();
				Particles.N_Particle[list_id1]->Calculate_Patch_Position();

				Reset_Positions(Particles, list_id2);	
				Particles.N_Particle[list_id2]->Calculate_Axis();
				Particles.N_Particle[list_id2]->Calculate_Patch_Position();

					
			}

			 	  // Accept move if XI <= b_factor 

			if (XI < b_factor){
				
			
				Set_Positions(Particles, id);	
				Set_Positions(Particles, list_id1);	
				Set_Positions(Particles, list_id2);	

				Particles.Total_Energy = Particles.Total_Energy + delta_U;
					
			}


		}
			 
	 
	} 
			
	


}		