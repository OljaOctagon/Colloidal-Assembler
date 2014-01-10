
# ifndef _polyhedra_h_included_
# define _polyhedra_h_included_

#define _USE_MATH_DEFINES	

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include "alg_struct.h"



class polyhedra{

public:
  
  
                   
double x_center; // center of mass coordinates
double y_center; 
double z_center; 

double* x;
double* y;
double* z;

m_vector trans_init;

double* dist_x;
double* dist_y;
double* dist_z;

double* new_dist_x;
double* new_dist_y;
double* new_dist_z;

double* trans_old;
double* Rot_old;  

double* trans_periodic;


double cut_off;
double cut_off_squared;
double cell_factor;
double V;

m_vector ax_1;
m_vector ax_2;
m_vector ax_3;

m_vector ax_1_old;
m_vector ax_2_old;
m_vector ax_3_old;

double ax_norm_x, ax_norm_y, ax_norm_z;
double p4;
double alpha, beta, gamma; //Euler angles
double alpha_old, beta_old, gamma_old;

/*
int MAX_cell_members;
int MAX_verlet_members;
int N_Neighbours;
int MAX_coll_partners;
int MAX_first_shell_partners;
*/


bool copy_count;

int left_count, right_count, top_count, bottom_count, back_count, front_count;
int cm_left_count, cm_right_count, cm_top_count, cm_bottom_count, cm_back_count, cm_front_count;
int cm_out;
int sum_edge_out;
int* edge_out; 

 
 
m_quaternion q;



double rot_angle;

m_vector a;
m_vector a_old;

double theta;
double theta_old;


int cell_out;
int edge_N;

void Set_Initial_Quaternion();
void Set_Initial_Center_Position();
void Set_Initial_Axis();



};



# endif
