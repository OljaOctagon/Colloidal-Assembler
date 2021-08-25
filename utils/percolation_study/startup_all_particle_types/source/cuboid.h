#ifndef _cuboid_h_included_
#define _cuboid_h_included_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "polyhedra.h"

using namespace std;

class cuboid : public polyhedra {

  public:
    cuboid();
    ~cuboid();
    double Lx, Ly, Lz;
    double V;

    void edges_from_center();
    void distance_from_center();
};

#endif
