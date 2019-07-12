#ifndef _sphere_h_included_
#define _sphere_h_included_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "algstruct.h"

using namespace std;

class sphere {

  public:
    sphere();

    ~sphere();

    double R;

    double x_center; // center of mass coordinates
    double y_center;
    double z_center;

    m_vector trans_init;

    double *trans_periodic;

    double cut_off;
    double V;

    int MAX_cell_members;
    int MAX_verlet_members;
    int N_Neighbours;
    int MAX_coll_partners;

    bool copy_count;

    int left_count, right_count, top_count, bottom_count, back_count,
        front_count;
    int cm_left_count, cm_right_count, cm_top_count, cm_bottom_count,
        cm_back_count, cm_front_count;
    int cm_out;

    int cell_out;
};

#endif
