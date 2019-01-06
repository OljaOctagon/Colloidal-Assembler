




	void pmove::Rot_Update_Quarternions_VON_MISES(particles& Particles, int id){
			
			
		m_quarternion mu;
		m_quarternion q_s;
		double b, x0, c,u,v,z, t,s;
		double kappa;
		
		kappa=100.0;
			
			
		mu.w= Particles.N_Particle[id].q.w;
		mu.x= Particles.N_Particle[id].q.x;
		mu.y= Particles.N_Particle[id].q.y;
		mu.z= Particles.N_Particle[id].q.z;
		
		kappa= 	
			
		b = -kappa + sqrt(kappa*kappa + 1);
		x0 = (1.0 - b)/(1.0 + b);
		c = kappa*x0 + 2*log(1 - x0*x0);	
		
		
		do{
			
			// Generate a random value with the beta(a, a) distribution with a = 1.5.
			//(1.5 because the parameters for beta are (m - 1)/2 and m = 4.

			do{
			  u = gsl_ran_flat(r01,-1.0,1.0);
			  v = gsl_rng_uniform(r01);

			  s = u*u + v*v;

			  }while(s>1);
			
			z = 0.5 + u*v*sqrt(1 - s)/s;
			
							
			u = gsl_rng_uniform (r01)
			w = (1 - (1 + b)*z)/(1 - (1 - b)*z);
			t = kappa*w + 2*log(1 - x0 * w) - c;
			
		}while(t<log(u));
				
			
		do {
			u0 = gsl_ran_flat(r01,-1.0,1.0);
			u1 = gsl_ran_flat(r01,-1.0,1.0);
						
			sq_u01 = u0*u0 + u1*u1;
				
		} while(sq_u01>1.0);
					
		sp_x = 2.*u0*sqrt(1.0-sq_u01);
		sp_y = 2.*u1*sqrt(1.0-sq_u01);
		sp_z = 1.0 - 2.*sq_u01;
		
		
		q_t.w = w;
		q_t.x = sp_x*sqrt(1 - w*w) 
		q_t.y= sp_y*sqrt(1 - w*w) 
		q_t.z = sp_zsqrt(1 - w*w) 
			
		
		q_s.w = q_t.w*mu.w - q_t.x*mu.x - q_t.y*mu.y - q_t.z*mu.z;
		q_s.x = q_t.w*mu.x + q_t.x*mu.w + q_t.y*mu.z - q_t.z*mu.x 
		q_s.y = q_t.w*mu.y - q_t.x*mu.z + q_t.y*mu.w + q_t.z*mu.x
		q_s.z = q_t.w*mu.z + q_t.x*mu.y - q_t.y*mu.x + q_t.x*mu.w;
		
		Particles.N_Particle[id].q.w = q_s.w;
		Particles.N_Particle[id].q.x = q_s.x;
		Particles.N_Particle[id].q.y = q_s.y;
		Particles.N_Particle[id].q.z = q_s.z;
		
		ofstream orient_out("quarternion_von_mises.dat",  ios::out | ios::app);
		orient_out<<Particles.N_Particle[id].q.w<<"  "<<Particles.N_Particle[id].q.x<<"  "<<Particles.N_Particle[id].q.y<<"  "<<Particles.N_Particle[id].q.z<<endl;
		orient_out.close();
		
		
	
		
	}	

/*
  References:
 
  [The paper by Fisher is the definitive paper on the distribution on a
  sphere. Mardia and Jupp discuss it and the generalization to
  n-dimensions, but refer to Wood for details of how to generate a set of
  samples from the distribution. Wood's work is based on the paper by
  Ulrich. The report by Dhillon and Sra was used to code the implementation
  below, but they gave an algorithm for the general m-dimensional case. We
  have simplified the code for the case of three dimensions and also made a
  correction -- their Figure 4 does not make use of c after computing it.
  Wood's paper reveals that it should be used in the test on Z and U.]
 
  R. A. Fisher, "Dispersion on a sphere", Proceedings of the Royal Society
  of London Series A., 217, pp295-305, (1953).
 
  K. V. Mardia and P. E. Jupp, "Directional Statistics" (2nd edition), John
  Wiley (2000). ISBN 0-471-95333-4. [§9.3.]
 
  Gary Ulrich, "Computer Generation of Distributions on the m-Sphere",
  Applied Statistics, 33(2), pp158-163, (1984).
 
  A. T. A. Wood, "Simulation of the von-Mises Distribution", Communications
  in statistics : simulation & computation, 23, pp157-164, 1994.
 
  Inderjit S. Dhillon and Suvrit Sra, "Modeling Data using Directional
  Distributions", Technical Report TR-03-06, Department of Computer
  Sciences, The University of Texas at Austin, Austin, TX 78712, USA.
  25 January, 2003.
 Accessed at: http://www.cs.utexas.edu/research/publications/ in May 2008.


*/



