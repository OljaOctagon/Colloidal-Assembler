#ifndef _octahedron_h_included_
#define _octahedron_h_included_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "polyhedra.h"

using namespace std;

class octahedron : public polyhedra {

  public:
    octahedron();
    ~octahedron();

    double height;

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
