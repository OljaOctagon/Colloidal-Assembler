#ifndef _move_h_included_
#define _move_h_included_

#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <math.h>

#include "box.h"
#include "particles.h"
//#include "box_non_ortho.h"
#include "alg_struct.h"
#include "fileio.h"
#include "polyhedra.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

//# include "cuboid.h"

using namespace std;

class pmove : public particles, public box {

  public:
    pmove();
    pmove(int size, int edge_length, double delta_tmax, double delta_rmax,
          double delta_Vmax, double Temperature, string is_translation_ON,
          string is_rotation_ON, string is_volumemove_ON,
          string is_grand_canonical_ON, string is_cluster_move_ON, int is_2D);
    ~pmove();

    gsl_rng *r;
    gsl_rng *r01;

    int random_init_r;
    int random_init_r01;

    double Rot_mat[9];
    double R_Rot_mat[9];

    double phi_reverse;

    m_vector trans_vec; // vector for translation
    m_vector zero_vec;

    double rand_x, rand_y, rand_z;
    double rand_V, rand_Va;
    double rand_Lx, rand_Ly, rand_Lz;
    double rand_alpha, rand_beta, rand_gamma;
    double delta_alpha, delta_gamma, delta_beta;

    int r_id;

    double XI;

    double rand_q_0, rand_q_1, rand_q_2, rand_q_3, rand_q_4;
    double dmax_alpha, dmax_beta, dmax_gamma;
    double dmax_t, dmax_q, dmax_V, dmax_L;
    double dmax_angle;
    double norm_q;
    double delta_q[4];
    double Nplus1;

    double rand_theta, rand_a1, rand_a2, rand_a3;
    double dmax_a, dmax_theta;
    double delta_theta;
    m_vector delta_a;
    double sigma_trans;

    double q_0, q_1, q_2, q_3; // quaternions helper vars

    double *P_mt;
    int mt_sum;
    int mt_w[5];
    size_t sizep;
    size_t value;
    gsl_ran_discrete_t *mt;
    string is_translation_ON;
    string is_rotation_ON;
    string is_volumemove_ON;
    string is_cluster_move_ON;
    string is_grand_canonical_ON;
    double beta_f;

    m_vector center_mass;
    int move_type;
    int movetype;
    int index1, index2, index3;
    int out_count;
    int coll_partner;
    int cp_id;

    // 2D

    int is_2D;

    // Calculate quaterinion distribution

    m_quaternion mu;
    m_quaternion q_s;
    double b, x0, c, u, v, z, w, t, s;
    double kappa;
    double sqrt_x;

    // Volmove

    double b_factor;
    double b_factor_pre;
    double kB;
    double T;
    double beta;

    double cos_a, sin_a;
    double sp_x, sp_y, sp_z;
    double sp_angle;
    double sp_ran;

    double Unity[9];

    m_vector trans;
    double Rot[9];
    int col_count;
    int exit_status;

    // Iterate
    int N_move_p;
    int move_p;
    int v_number;

    double m_scalar_product;
    int collision_partner;
    double dL;
    double R_id;
    double R_collision_partner;
    double f_exit;
    double f_exit_MAX;
    int collision_count;
    double L_axis_norm2;

    int accept_translate, accept_rotate, accept_iso_vol, accept_aniso_vol,
        accept_complete;
    int N_trans_moves, N_rotate_moves, N_iso_vol_moves, N_aniso_vol_moves;
    double accept_translate_procent, accept_rotate_procent,
        accept_iso_vol_procent, accept_aniso_vol_procent,
        accept_complete_procent;
    double accept_translate_procent_MAX, accept_rotate_procent_MAX,
        accept_iso_vol_procent_MAX, accept_aniso_vol_procent_MAX;
    double accept_translate_procent_MIN, accept_rotate_procent_MIN,
        accept_iso_vol_procent_MIN, accept_aniso_vol_procent_MIN;

    int cluster_counter;
    int cluster_size;
    double cluster_radius;
    double cluster_radius_square;
    double diffc_x, diffc_y, diffc_z;
    double diffc_square;
    int *Cluster_List;

    double square_term;
    double sq_u01;
    double sq_u23;
    double u0, u1, u2, u3;

    m_quaternion q_t;
    m_quaternion q_diff;

    double q_diff_norm;
    double q_MAX;

    // Minimum
    double return_value;

    m_vector *L_axis; // only for cubes
    m_vector particle_distance;

    double u_x, u_y, u_z;

    // Patches

    double Total_Energy;
    double Total_Energy_old;
    double Halo_Energy;
    double Energy_Neighbours_id;

    double **Patch_Energies;
    double **Patch_Energies_old;
    int potential_exit_status;

    double patch_dist_x_ti;
    double patch_dist_y_ti;
    double patch_dphi_ti;

    double *patch_dist_x;
    double *patch_dist_y;
    double *patch_dphi;


    // Cluster move

    int N_List;
    int N_Bonds;
    int N_failed_links;
    int **pseudo_cluster_info;
    int *List;
    bool *is_element;
    double *p_f;
    double *p_r;
    double q_factor;

    double *q_f;
    double *q_r;

    m_vector *map_dist;
    m_vector map_dist_center;
    m_vector map_rot_center;

    m_vector r_2trans_vec;
    m_vector f_2trans_vec;
    m_vector r_trans_vec;

    double e_old, e_new, e_pair, e_single;

    double minimum(double a, double b);
    void Iterate(particles &Particles, box *Box, fileio &Fileio, int mc_time);
    void Calculate_Acceptances(int mc_time);
    void Calculate_Separating_Axis(particles &Particles, int id, int j,
                                   m_vector *L_axis);
    double Scalar_Product(m_vector a, m_vector b);

    void Trans_Update_Positions(particles &Particles, int id,
                                m_vector &trans_vec);
    void Rot_Update_Positions(particles &Particles, int id,
                              double (&Rot_mat)[9]);
    void Rot_Update_Quarternions_RANDOM(particles &Particles, int id);
    void Rot_Update_Axis_Angle_RANDOM(particles &Particles, int id);
    void Rot_Update_Quarternions_VON_MISES(particles &Particles, int id);

    void Update_Periodic_Positions(int id, double *periodic_vector, double Lx,
                                   double Ly, double Lz);
    void Translate(particles &Particles, box *Box, fileio &Fileio, int id,
                   int mc_time);
    void Rotate(particles &Particles, box *Box, fileio &Fileio, int id,
                int mc_time);

    void Iso_Vol_Change(particles &Particles, box *Box, fileio &Fileio,
                        int mc_time);
    void Aniso_Vol_Change(particles &Particles, box *Box, fileio &Fileio,
                          int mc_time);

    int Random(int N);
    void Set_Random_State(gsl_rng *r_start, gsl_rng *r01_start);

    void Collide_Periodic(particles &Particles, box *Box, fileio &Fileio,
                          int id, int time, collision_list *Collision_List);
    void Col_Test(particles &Particles, box *Box, int id,
                  collision_list *Collision_List);
    void Collision_Test(particles &Particles, box *Box, int id,
                        collision_list *Collision_List);
    void Collision_Test_Halo(particles &Particles, box *Box, int id,
                             collision_list *Collision_List);

    double Calculate_Cross_Product_x(m_vector a, m_vector b);
    double Calculate_Cross_Product_y(m_vector a, m_vector b);
    double Calculate_Cross_Product_z(m_vector a, m_vector b);
    int Calculate_Separating_Axis_GENERAL(particles &Particles, int id, int j,
                                          m_vector *L_axis);

    int Calculate_Separating_Axis_RHOMBI(particles &Particles, int id,
                                         int j, m_vector *L_axis);

    void Update_Periodic_Positions(particles &Particles, box *Box, int id);
    void Reset_Positions(particles &Particles, int id);
    void Set_Positions(particles &Particles, int id);

    void Calculate_Cluster_List(particles &Particles, box *Box);

    void Trans_Update_Positions_SPHERE(particles &Particles, int id,
                                       m_vector &trans_vec);
    void Update_Periodic_Positions_SPHERE(particles &Particles, box *Box,
                                          int id);
    void Reset_Positions_SPHERE(particles &Particles, int id);
    void Set_Positions_SPHERE(particles &Particles, int id);
    void Collision_Test_SPHERE(particles &Particles, box *Box, int id,
                               collision_list *Collision_List);

    // for 2D system

    void Rot_Update_Random_2D(particles &Particles, int id);
    void Rotate2D(particles &Particles, box *Box, fileio &Fileio, int id,
                  int mc_time);

    void Reset_Pair_Potential(particles &Particles, box *Box);
    void Set_Pair_Potential(particles &Particles, box *Box);
    void Reset_Pair_Potential(int id, particles &Particles, box *Box);
    void Set_Pair_Potential(int id, particles &Particles, box *Box);

    double Calculate_Pair_Potential(int id, particles &Particles, box *Box,
                                    collision_list *Collision_List);
    void Calculate_Pair_Potential(particles &Particles, box *Box,
                                  collision_list *Collision_List);
    void Calculate_Pair_Potential(particles &Particles, box *Box);
    double Calculate_Patch_Distance(int id1, int id2, int pid1, int pid2,
                                    particles &Particles, box *Box);
    void Calculate_Pair_Potential_Neighbours(particles &Particles, box *Box,
                                             int id1);

    void Particle_Insertion(particles &Particles, box *Box, fileio &Fileio,
                            int mc_time);
    void Particle_Deletion(particles &Particles, box *Box, fileio &Fileio,
                           int mc_time);

    double Calculate_Potential(int a, int b, particles &Particles, box *Box);
    void Reset_Pseudo_Cluster(box *Box);
    void Trans_Pseudocluster_Recursion(int id_j, int cn, int fl,
                                       particles &Particles, box *Box);
    void Trans_Cluster_Move(particles &Particles, box *Box, fileio &Fileio,
                            int mc_time);

    void Pseudocluster_Recursion(int id_j, int cn, int fl,
                                 particles &Particles, box *Box);
    void Rot_Cluster_Move(particles &Particles, box *Box, fileio &Fileio,
                          int mc_time);
    void Check_Periodic_Center_of_Mass(m_vector &center_mass, box *Box);

    void Three_Particle_Move(particles &Particles, box *Box, fileio &Fileio,
                             int mc_time);
    void Rot_Move_Map(particles &Particles, int id1, int id2, box *Box,
                      double (&Rot_mat)[9]);
    void Rot_Move_Map(particles &Particles, int id1, box *Box,
                      m_vector &center_mass, double (&Rot_mat)[9]);
};

#endif
