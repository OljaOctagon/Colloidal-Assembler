
#include "box_non_ortho.h"

box::box() {

    packing_fraction = 0.0;
    P_sigma = 0.0;
    N = 0;
    T = 0.0;
    P = 0.0;
    V = 0.0;

    Lx = 0.0;
    Ly = 0.0;
    Lz = 0.0;

    x_center = 0.0;
    y_center = 0.0;
    z_center = 0.0;

    edge_N = 8;

    x = NULL;
    y = NULL;
    z = NULL;

    dist_x = NULL;
    dist_y = NULL;
    dist_z = NULL;

    new_dist_x = NULL;
    new_dist_y = NULL;
    new_dist_z = NULL;

    x = new double[edge_N];
    y = new double[edge_N];
    z = new double[edge_N];

    dist_x = new double[edge_N];
    dist_y = new double[edge_N];
    dist_z = new double[edge_N];

    new_dist_x = new double[edge_N];
    new_dist_y = new double[edge_N];
    new_dist_z = new double[edge_N];
}

void box::Startconfig(int N_in, double P_sigma_in, double packing_fraction_in,
                      double V_p) {

    packing_fraction = packing_fraction_in;
    P_sigma = P_sigma_in;
    N = N_in;
    // V_p = 1.0;

    T = 1.0;
    // P = P_sigma/(3.*sqrt(3.)/8.);
    P = P_sigma;
    // V = double(N)*V_p/packing_fraction;

    V = double(N) * V_p / packing_fraction;

    cout << "V: " << V << " "
         << "V_p:" << V_p << endl;

    // for cubes
    Lx = pow(V, (1. / 3.));
    Ly = Lx;
    Lz = Lx;

    x_center = Lx / 2.0;
    y_center = Ly / 2.0;
    z_center = Lz / 2.0;

    v1.x = Lx;
    v1.y = 0.0;
    v1.z = 0.0;

    v2.x = 0.0;
    v2.y = Ly;
    v2.z = 0.0;

    v3.x = 0.0;
    v3.y = 0.0;
    v3.z = Lz;

    alpha = M_PI / 2.0;
    beta = M_PI / 2.0;
    gamma = M_PI / 2.0;

    vc.x = v1.x + v2.x + v2.x;
    vc.y = v1.y + v2.y + v2.y;
    vc.z = v1.z + v2.z + v2.z;

    cut_off = vc.norm();
    cut_off_squared = cut_off * cut_off;
}

void box::Startconfig_former(int N_in, double P_sigma_in) {

    P_sigma = P_sigma_in;
    N = N_in;
    // V_p = 1.0;
    T = 1.0;
    // P = P_sigma/(3.*sqrt(3.)/8.);
    P = P_sigma;
}

box::~box() {

    packing_fraction = 0.0;
    P_sigma = 0.0;
    N = 0;
    T = 0.0;
    P = 0.0;
    V = 0.0;

    Lx = 0.0;
    Ly = 0.0;
    Lz = 0.0;
}

void box::edges_from_center() { // only valid if cubes are aligned with
                                // coordinate axes!

    x[0] = x_center + (-v1.x - v2.x - v3.x) / 2.0;
    y[0] = y_center + (-v1.y - v2.y - v3.y) / 2.0;
    z[0] = z_center + (-v1.z - v2.z - v3.z) / 2.0;

    x[1] = x_center + (v1.x - v2.x - v3.x) / 2.0;
    y[1] = y_center + (v1.y - v2.y - v3.y) / 2.0;
    z[1] = z_center + (v1.z - v2.z - v3.z) / 2.0;

    x[2] = x_center + (v1.x + v2.x - v3.x) / 2.0;
    y[2] = y_center + (v1.y + v2.y - v3.y) / 2.0;
    z[2] = z_center + (v1.z + v2.z - v3.z) / 2.0;

    x[3] = x_center + (-v1.x + v2.x - v3.x) / 2.0;
    y[3] = y_center + (-v1.y + v2.y - v3.y) / 2.0;
    z[3] = z_center + (-v1.z + v2.z - v3.z) / 2.0;

    x[4] = x_center + (-v1.x - v2.x + v3.x) / 2.0;
    y[4] = y_center + (-v1.y - v2.y + v3.y) / 2.0;
    z[4] = z_center + (-v1.z - v2.z + v3.z) / 2.0;

    x[5] = x_center + (v1.x - v2.x + v3.x) / 2.0;
    y[5] = y_center + (v1.y - v2.y + v3.y) / 2.0;
    z[5] = z_center + (v1.z - v2.z + v3.z) / 2.0;

    x[6] = x_center + (v1.x + v2.x + v3.x) / 2.0;
    y[6] = y_center + (v1.y + v2.y + v3.y) / 2.0;
    z[6] = z_center + (v1.z + v2.z + v3.z) / 2.0;

    x[7] = x_center + (-v1.x + v2.x + v3.x) / 2.0;
    y[7] = y_center + (-v1.y + v2.y + v3.y) / 2.0;
    z[7] = z_center + (-v1.z + v2.z + v3.z) / 2.0;
}

void box::distance_from_center() {

    for (int j = 0; j < edge_N; j++) {

        dist_x[j] = x[j] - x_center;
        dist_y[j] = y[j] - y_center;
        dist_z[j] = z[j] - z_center;
    }

    // cout<<"Distance from center"<<endl;
}

void box::Set_Lengths(m_vector v1, m_vector v2, m_vector v3) {

    Lx = v1.norm();
    Ly = v2.norm();
    Lz = v3.norm();

    alpha = acos((v1.x * v3.x + v1.y * v3.y + v1.z * v3.z) / (Lx * Lz));
    beta = acos((v2.x * v3.x + v2.y * v3.y + v2.z * v3.z) / (Ly * Lz));
    gamma = acos((v1.x * v2.x + v1.y * v2.y + v1.z * v2.z) / (Lx * Ly));

    vc.x = v1.x + v2.x + v2.x;
    vc.y = v1.y + v2.y + v2.y;
    vc.z = v1.z + v2.z + v2.z;

    cut_off = vc.norm();
    cut_off_squared = cut_off * cut_off;

    cross_p.x = v2.y * v3.z - v2.z * v3.y;
    cross_p.y = v2.z * v3.x - v2.x * v3.z;
    cross_p.z = v2.x * v3.y - v2.y * v3.x;

    V = v1.x * cross_p.x + v1.y * cross_p.y + v1.z * cross_p.z;
}

void box::Set_Lengths(double Vp) {

    // for isotropic box

    V = double(N * Vp) / packing_fraction;
    Lx = pow(V, (1. / 3.));
    Ly = Lx;
    Lz = Lx;

    v1.x = Lx;
    v1.y = 0.0;
    v1.z = 0.0;

    v2.x = 0.0;
    v2.y = Ly;
    v2.z = 0.0;

    v3.x = 0.0;
    v3.y = 0.0;
    v3.z = Lz;

    alpha = M_PI / 2.0;
    beta = M_PI / 2.0;
    gamma = M_PI / 2.0;

    vc.x = v1.x + v2.x + v2.x;
    vc.y = v1.y + v2.y + v2.y;
    vc.z = v1.z + v2.z + v2.z;

    cut_off = vc.norm();
    cut_off_squared = cut_off * cut_off;
}
