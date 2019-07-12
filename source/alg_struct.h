

#ifndef _alg_struct_h_included_
#define _alg_struct_h_included_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>

class m_vector {

  public:
    double x, y, z;

    m_vector();
    ~m_vector();

    double norm();
    double norm_m;
};

class m_quaternion {

  public:
    double x, y, z, w;

    m_quaternion();
    ~m_quaternion();

    double norm_q;

    double norm();
};

class data_point_2D {

  public:
    double x, y;
};

#endif
