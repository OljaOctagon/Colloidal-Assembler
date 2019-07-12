#include "cuboid.h"

cuboid::cuboid() {

    edge_N = 8;

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

cuboid::~cuboid() {

    delete[] x;
    delete[] y;
    delete[] z;

    delete[] dist_x;
    delete[] dist_y;
    delete[] dist_z;

    delete[] new_dist_x;
    delete[] new_dist_y;
    delete[] new_dist_z;
}

void cuboid::edges_from_center() { // only valid if cubes are aligned with
                                   // coordinate axes!

    x[0] = x_center - 0.5 * Lx;
    y[0] = y_center - 0.5 * Ly;
    z[0] = z_center - 0.5 * Lz;

    x[1] = x_center + 0.5 * Lx;
    y[1] = y_center - 0.5 * Ly;
    z[1] = z_center - 0.5 * Lz;

    x[2] = x_center + 0.5 * Lx;
    y[2] = y_center + 0.5 * Ly;
    z[2] = z_center - 0.5 * Lz;

    x[3] = x_center - 0.5 * Lx;
    y[3] = y_center + 0.5 * Ly;
    z[3] = z_center - 0.5 * Lz;

    x[4] = x_center - 0.5 * Lx;
    y[4] = y_center - 0.5 * Ly;
    z[4] = z_center + 0.5 * Lz;

    x[5] = x_center + 0.5 * Lx;
    y[5] = y_center - 0.5 * Ly;
    z[5] = z_center + 0.5 * Lz;

    x[6] = x_center + 0.5 * Lx;
    y[6] = y_center + 0.5 * Ly;
    z[6] = z_center + 0.5 * Lz;

    x[7] = x_center - 0.5 * Lx;
    y[7] = y_center + 0.5 * Ly;
    z[7] = z_center + 0.5 * Lz;
}

void cuboid::distance_from_center() {

    for (int j = 0; j < 8; j++) {

        dist_x[j] = x[j] - x_center;
        dist_y[j] = y[j] - y_center;
        dist_z[j] = z[j] - z_center;
    }
}
