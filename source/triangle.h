#ifndef _triangle_h_included_
#define _triangle_h_included_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "polyhedra.h"

using namespace std;

class triangle : public polyhedra {

  public:
    double alpha;
    double sinus;
    double L;
    double H;
    double L_2;
    double H_2;
    double a_x;
    double h;
    double h_2;
    double diag2_short, diag2_long;
    double diag3_short, diag3_long;
    double p, q;
    double px_l, py_l, px_s, py_s;
    double qx_l, qy_l, qx_s, qy_s;
    double patch_delta;
    m_vector ra;
    m_vector rb;
    m_vector rc;
    double d0, d1, d2, d3, d4, d5;
    double patch_x;
    double patch_size;

    m_vector vc; 
    m_vector cross_p;
    string triangle_type;
    triangle();
    ~triangle();

    void edges_from_center();
    void distance_from_center();
    void Calculate_Axis();
    void Set_Axis();
    void Set_Lengths();
    void Set_Lengths(int e0, int e1, int e2, int e3, int e4, int e5);
    void Set_Start_Lattice_Position(int id, double box_Lx, int N_box);
    void Set_Start_Lattice_Position(int id, double box_Lx, double box_Ly,
                                    double box_Lz, int N_box);
    void Calculate_Face_Normals();
    double Calculate_Projection_to_Separating_Axis(m_vector laxis);
    void Calculate_Patch_Position();
    void Calculate_Long_Axis();
};

#endif

