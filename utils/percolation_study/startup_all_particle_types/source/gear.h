#ifndef _gear_h_included_
#define _gear_h_included_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "cube.h"
#include "polyhedra.h"

using namespace std;


class gear : public polyhedra {

  public:
    gear();
    ~gear();

    int Ng;
    cube *gear_wheels;
    double R;
    double Lx_max;
    double u;
    double ds;
    double h_t;
    double h;
    double p_Lx;
    double d;
    double As, At, Ap;

    m_vector *base_v;
    double phi_n;
    double d_phi;

    double Rot_mat[9];
    m_vector *rel_pos;

    void get_initial_relative_wheel_positions(int Ng);
    void get_circle_partition(int Ng);
    void get_relative_orientation();
    void Update_Rotation();

    void write_positions(ofstream &fout);
    void edges_from_center();
    void distance_from_center();
    void Calculate_Axis();
    void Set_Axis();
    void Set_Lengths();
    void Set_Start_Lattice_Position(int id, double box_Lx, int N_box);
    void Set_Start_Lattice_Position(int id, double box_Lx, double box_Ly,
                                    double box_Lz, int N_box);

    void Calculate_Face_Normals();
    double Calculate_Projection_to_Separating_Axis(m_vector laxis);
};

#endif
