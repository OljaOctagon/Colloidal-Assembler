#ifndef _parallelepiped_h_included_
#define _parallelepiped_h_included_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "polyhedra.h"

using namespace std;

class parallelepiped : public polyhedra {

  public:
    m_vector v1;
    m_vector v2;
    m_vector v3;

    m_vector vc;
    m_vector cross_p;

    parallelepiped();
    ~parallelepiped();

    void write_positions(ofstream &fout);
    void edges_from_center();
    void distance_from_center();
    void Calculate_Axis();
    void Set_Axis();
    void Set_Lengths();
    void Set_Lengths(double Lx, double Ly, double Lz, double alpha,
                     double gamma);
    void Set_Start_Lattice_Position(int id, double box_Lx, int N_box);
    void Calculate_Face_Normals();
    double Calculate_Projection_to_Separating_Axis(m_vector laxis);
};

#endif
