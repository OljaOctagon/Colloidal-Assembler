
#ifndef _polyhedra_h_included_
#define _polyhedra_h_included_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "alg_struct.h"

using namespace std;

class polyhedra {

  public:
    double x_center; // center of mass coordinates
    double y_center;
    double z_center;

    double *x;
    double *y;
    double *z;

    double *x_patch;
    double *y_patch;
    double *z_patch;

    double halo_energy;
    double halo_cutoff;

    int N_patches;
    int N_patch_types;
    int N_active;

    double *r_patch;
    double *patch_cutoff;
    double *patch_cutoff_squared;
    double **patch_energy;
    int *patch_type;

    double Lx;
    double Ly;
    double Lz;

    // Inner Radius
    double r;
    double cut_off_inner;
    double V;
    double A;

    m_vector trans_init;

    m_vector *rel_dist;

    double *dist_x;
    double *dist_y;
    double *dist_z;

    double *new_dist_x;
    double *new_dist_y;
    double *new_dist_z;

    double *trans_old;
    double *Rot_old;

    double *trans_periodic;

    double cut_off;
    double cut_off_small;
    double cut_off_squared;
    double cut_off_small_squared;
    double cell_factor;

    m_vector ax_1;
    m_vector ax_2;
    m_vector ax_3;

    m_vector long_axis;

    m_vector ax_1_old;
    m_vector ax_2_old;
    m_vector ax_3_old;

    double ax_norm_x, ax_norm_y, ax_norm_z;
    double p4;
    double alpha, beta, gamma; // Euler angles
    double alpha_old, beta_old, gamma_old;

    bool copy_count;

    int left_count, right_count, top_count, bottom_count, back_count,
        front_count;
    int cm_left_count, cm_right_count, cm_top_count, cm_bottom_count,
        cm_back_count, cm_front_count;
    int cm_out;
    int sum_edge_out;
    int *edge_out;

    m_quaternion q;

    double phi;
    double phi_old;

    double rot_angle;

    m_vector a;
    m_vector a_old;

    double theta;
    double theta_old;
    int N_all;

    int cell_out;
    int edge_N;

    double Rp;
    double rmax;
    double rmin;

    int N_independent_faces;
    int N_cross_edges;
    int edge_N_vec;

    m_vector *facenormal;
    m_vector *edges;
    m_vector *edge_cross;

    void Set_Initial_Quaternion();
    void Set_Initial_Center_Position();
    void Set_Initial_Axis();
    void Calculate_Quaternions_2D();

    virtual void write_positions(ofstream &fout) {}
    virtual void edges_from_center() {}
    virtual void distance_from_center() {}
    virtual void Calculate_Axis() {}
    virtual void Set_Axis() {}
    virtual void Set_Lengths() {}
    virtual void Set_Lengths(int a, int b, int c, int d, int e, int f) {}
    virtual void Set_Start_Lattice_Position(int id, double box_Lx, int N_box) {
    }
    virtual void Set_Start_Lattice_Position(int id, double box_Lx,
                                            double box_Ly, double box_Lz,
                                            int N_box) {}
    virtual void Calculate_Face_Normals() {}
    virtual double Calculate_Projection_to_Separating_Axis(m_vector laxis);
    virtual void Calculate_Patch_Position(){};
    virtual void Calculate_Long_Axis(){};
};

#endif
