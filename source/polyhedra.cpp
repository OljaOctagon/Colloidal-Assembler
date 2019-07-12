#include "polyhedra.h"

void polyhedra::Set_Initial_Quaternion() {

    q.x = 0.0;
    q.y = 0.0;
    q.z = 0.0;
    q.w = 1.0;

    phi = 0.0;

    // N_patches=2;
}

void polyhedra::Set_Initial_Axis() {

    ax_1.x = 0.0;
    ax_1.y = 0.0;
    ax_1.z = 0.0;

    ax_2.x = 0.0;
    ax_2.y = 0.0;
    ax_2.z = 0.0;

    ax_3.x = 0.0;
    ax_3.y = 0.0;
    ax_3.z = 0.0;

    long_axis.x = 0.0;
    long_axis.y = 0.0;
    long_axis.z = 0.0;
}

void Calculate_Quaternions_2D() {

    // NEEDS TO BE FILLED :)
}

double polyhedra::Calculate_Projection_to_Separating_Axis(m_vector laxis) {
    double Rp;
    Rp = 0;

    std::cout << "Calculate Projection to Separating Axis POLYHEDRA"
              << std::endl;

    return Rp;
}
