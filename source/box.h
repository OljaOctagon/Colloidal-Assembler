

#ifndef _box_h_included_
#define _box_h_included_

#define _USE_MATH_DEFINES

#include "cuboid.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>

using namespace std;

class box : public polyhedra {

  public:
    box();
    ~box();

    void edges_from_center();
    void distance_from_center();
    void Startconfig(int N_in, double P_sigma_in, double mu_in_1, double mu_in_2,
                     double packing_fraction_in, double V_p);
    void Startconfig_former(int N_in, double P_sigma_in, double mu_in_1, double mu_in_2);
    void Set_Lengths(double Vp);
    void Calculate_Box_Axis();

    m_vector box_ax_1;
    m_vector box_ax_2;
    m_vector box_ax_3;
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
    double mu_1;
    double mu_2;
    double packing_fraction;
    double V;
    double V_p;
    int N_c;
    double N_c_float;
};

#endif
