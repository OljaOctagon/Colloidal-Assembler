#ifndef _truncatedcube_h_included_
#define _truncatedcube_h_included_

#define _USE_MATH_DEFINES

#include "polyhedra.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace std;

class truncated_cube : public polyhedra {

  public:
    truncated_cube();
    ~truncated_cube();

    double edge_R;
    double lx, ly, lz;

    void write_positions(ofstream &fout);
    void edges_from_center();
    void distance_from_center();
    void Calculate_Axis();
    void Set_Axis();
    void Set_Start_Lattice_Position(int id, double box_Lx, int N_box);
    void Set_Lengths();
    double Calculate_Projection_to_Separating_Axis(m_vector laxis);
    void Calculate_Face_Normals();
};

#endif
