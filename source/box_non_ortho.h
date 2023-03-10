
#ifndef _box_non_ortho_h_included_
#define _box_non_ortho_h_included_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
//#include "cuboid.h"
#include "polyhedra.h"

using namespace std;

class box : public polyhedra {

  public:
    box();
    ~box();

    void edges_from_center();
    void distance_from_center();
    void Startconfig(int N_in, double P_sigma_in, double packing_fraction_in,
                     double V_p);
    void Startconfig_former(int N_in, double P_sigma_in);
    void Set_Lengths(double Vp);
    void Set_Lengths(m_vector v1, m_vector v2, m_vector v3);

    m_vector v1;
    m_vector v2;
    m_vector v3;

    m_vector v1_old;
    m_vector v2_old;
    m_vector v3_old;

    m_vector vc;
    m_vector cross_p;
    double alpha, beta, gamma;
    double cut_off;
    double cut_off_squared;

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
    // double V_p;
    int N_c;
    double N_c_float;
};

#endif
