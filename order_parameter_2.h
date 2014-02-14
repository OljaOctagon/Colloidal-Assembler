

# ifndef _order_param_h_included_
# define _order_param_h_included_


#include <complex>
#include <cmath>
//#include <math.h>
#include <limits>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_histogram.h>
//#include <gsl/gsl_complex.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>


#include "cube.h"
#include "particles.h"
#include "box.h"
#include "alg_struct.h"


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using namespace std;

	class spher_harm_op{

		public:

		
		complex<double> complex_exp_value;
		complex<double>*** Yij_m;
		complex<double>** q_id_m;
		complex<double>** q_id_m_av;
		
		double* ql;

		double Q_l;
		double dist_cos_theta;
		double dist_phi;
		double legendre_spher_harm;
		double pre_factor;
		int m_double;
		int nb_id;
		int b_count;
		int number_of_bonds;
		int list_j;
		double Q_global;
		int MAX_coll_p;
		double qcorr_T;
		
		
		double* q_id_norm;
		complex<double>** qq_id_m;
		complex<double> qscalar_id_j;
		double** ql_ij_scalar;
		int* N_qbonds;
					
		int l_op;
		int m_total;		
		void Set(int l, int N, int MAX_coll_p_in, double qcorr_T_in);
		void Calculate_local_op(particles& Particles, box* Box);
		void Calculate_Corr_Funct(particles& Particles, box* Box);
		

			
	};


	class order_parameter{


		public:

		order_parameter();
		order_parameter(box* Box, particles& Particles, int MAX_coll_p_in);
		~order_parameter();

		int MAX_coll_p;


		double dist_x2;
		double dist_y2;
		double dist_z2;
		double dist;
		
		int list_j;
		double particle_dist, particle_dist_x, particle_dist_y, particle_dist_z;
		double particle_dist_squared;
		
		//local order parameters general 
		double cutoff_coll_list;
		double cutoff_op;
	
	
		// g_r Calculation
		double r_start;
		double r_stop;
		double max_r;
		double r_min;
		double r_dist1, r_dist2;
		double delta_r;
		int N_histo_points;
		double r_dist;
		double g_norm;
		double r_n;
		int* g_radial_histo;
		int* gx_histo;
		int* gy_histo;
		int* gz_histo;
		data_point_2D* g_radial_distribution;
		data_point_2D* gx_distribution;
		data_point_2D* gy_distribution;
		data_point_2D* gz_distribution;
		double Cut_Off;
		double g_Cut_Off;
		double g_min;
		int p_min;
		double num_r1, num_r2, rn;
		double interval_r1;
		double interval_r2;
		double theta_id, theta_j;
		m_vector rv_dist, newcoord_distance, nbd; 
		
		//spherical harmoincs
		spher_harm_op q4;
		spher_harm_op q6;
		double qcorr_T_in;

		//Cluster recognition
		int* is_nPhase;
		//char Cutoff_Type[128];
		string Cutoff_Type;
		
		double phase_T;
	
		//op3
		double a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z;
		double b1x, b1y, b1z, b2x, b2y, b2z, b3x, b3y, b3z;
		double cubic_op3;
		double scp[9];
		double cubic_op3_T;
		double cubic_op4_T;
		
		double* cubic_op3_av;
		double scpmax1, scpmax2;
		double scpmax;
		double q_scp;
		double cutoff_cubic_op3;
		double cubic_op3_order;
		int* N_obonds;
		
		
		double* cubic_op4_av;
		int* N_obonds_cubic4;
		double cubic_op4;
		
		
		//histogram
		size_t N_bins;
		gsl_histogram * histogram; 
		double xmin, xmax;
		int set_state, close_state, save_state;
		int update_state;

		//Functions
		void Calculate_Local_Order_Parameters(particles& Particles, cube* N_Particle, box* Box);
		void Caculate_is_nPhase(box* Box);
		double Calculate_g_r(particles& Particles, box* Box);
		void Calculate_g_xyz(particles& Particles, box* Box);
		void Calculate_cubic_op3(particles&, box* Box);
		void Calculate_cubic_op4(int Np);
		void unset_histogram();
		void update_histogram(double hv);
		void write_histogram_info();

	


		/*
		void g_r();
		void g_o();
		void g_xyz();
		*/ 


	};


# endif
