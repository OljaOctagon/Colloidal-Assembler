#include <stdio.h>

#include "alg_struct.h"

m_vector::m_vector() {

    x = 0;
    y = 0;
    z = 0;
}

m_vector::~m_vector() {}

double m_vector::norm() {

    norm_m = sqrt(x * x + y * y + z * z);

    if (norm_m <= 1.0e-14) {

        norm_m = 1.0e-14;
    }

    return norm_m;
}

m_quaternion::m_quaternion() {

    x = 0;
    y = 0;
    z = 0;
    w = 0;
}

m_quaternion::~m_quaternion() {}

double m_quaternion::norm() {

    norm_q = sqrt(x * x + y * y + z * z + w * w);

    return norm_q;
}
