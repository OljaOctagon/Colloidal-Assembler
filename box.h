

# ifndef _box_h_included_
# define _box_h_included_


#define _USE_MATH_DEFINES	

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstring>
#include "cuboid.h"


using namespace std;

class box : public cuboid {

public:

box();
~box();

void edges_from_center();
void distance_from_center();
void Startconfig(int N_in, double P_sigma_in, double packing_fraction_in, double V_p);
void Startconfig_former(int N_in, double P_sigma_in);
void Set_Lengths(double Vp);

double Lx;
double Ly;
double Lz;
double Lx_scale, Ly_scale, Lz_scale;
double VL_rel;
double V_rel;
double Lx_old, Ly_old, Lz_old, V_old;
double g_factor;
double g_trans_vec;
double P;
double P_sigma;
double T;
int N;
double packing_fraction;
double V;          
//double V_p;
int N_c;
double N_c_float;

};

#endif
