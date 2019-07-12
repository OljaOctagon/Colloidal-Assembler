

#ifndef _fileio_h_included_
#define _fileio_h_included_

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

#include "box.h"
#include "cube.h"
#include "particles.h"
//#include "box_non_ortho.h"
#include "cluster.h"
#include "gear.h"
#include "order_parameter.h"
#include "polyhedra.h"

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

using namespace std;

class fileio : public particles {

  public:
    string string_1;
    string string_g;
    string string_2;
    string string_3;
    string string_m;
    string string_4;
    string string_e;
    string string_5;

    string string_v_list;
    string string_V_List;
    string string_txt;
    string string_cell;
    string string_Cell;
    string string_C_List;
    string string_c_list;
    string string_cut_off;
    string string_Cut_Off;
    string string_Box;
    string string_para_out;
    string string_para_in;
    string string_orientations;
    string string_positions;
    string string_D2;
    string string_D4;
    string string_i2;
    string string_Vol_Packing;
    string string_acceptances;
    string string_NPT;
    string string_cluster_list;
    string string_local_OP;
    string string_g_r;
    string string_dat;
    string string_bin;
    string string_random1;
    string string_random2;

    int M;

    double P_sigma_in;
    int N_in;
    double mu_in;
    double packing_fraction_in;
    int number_of_cells_in;
    int MAX_coll_p_in;
    int MAX_fshell_p_in;

    string binary_on_in;
    double phi_binary_in;

    int read_state;
    int save_state;
    int close_state;

    int checkpoint_frequency_in;

    int time_in;
    int calc_frequency_in;
    int frame_frequency_in;
    double delta_tmax;
    double delta_rmax;
    double delta_Vmax;
    int seed1_in;
    int seed2_in;

    string is_translation_ON;
    string is_rotation_ON;
    string is_volumemove_ON;
    string is_grand_canonical_ON;
    string is_cluster_move_ON;

    int is_2D;

    string particle_type_in;

    double Temperature;

    void Set_Vars();

    void Write_Positions(int mc_time, box *Box, particles &Particles);
    void Write_Positions_Gears(int mc_time, box *Box, particles &Particles);
    void Write_Positions_SPHERE(int mc_time, box *Box, particles &Particles);
    void Write_Cell_Lists(int id, int mc_time, particles &Particles);
    void Write_Cell_Positions(int number_of_cells, cell *Cell, int mc_time,
                              box *Box);
    void Write_Order_Parameters(int mc_time, double D2_order, double D4_order);
    void Write_Vol_Packing(int mc_time, box *Box);
    void Write_Local_Order_Parameters(int mc_time, particles &Particles,
                                      box *Box,
                                      order_parameter &Order_Parameter,
                                      cluster &Cluster);
    void Write_Cluster_List(int *Cluster_List, int cluster_size,
                            particles &Particles,
                            order_parameter &Order_Parameter);
    void Write_Box_Positions(int mc_time, box *Box);
    void Write_g_r(int mc_time, order_parameter &Order_Parameter);

    void Save_Config(box *Box, particles &Particles, gsl_rng *r, gsl_rng *r01,
                     double dmax_t, double dmax_alpha, double dmax_beta,
                     double dmax_gamma, double dmax_V, double dmax_L,
                     int mc_time);
    void Save_Box(box *Box, int mc_time);
    void Save_Positions(box *Box, particles &Particles, int mc_time);
    void Save_Orientations(box *Box, particles &Particles, int mc_time);
    void Save_Orientations_EULER(box *Box, particles &Particles);
    void Save_Orientations_Axis_Angle(box *Box, particles &Particles,
                                      int mc_time);
    void Save_Random_State(gsl_rng *r, gsl_rng *r01, int mc_time);
    void Save_Patches(box *box, particles &Particles, int mc_time);

    void Read_Config(box *Box, particles &Particles, gsl_rng *r, gsl_rng *r01,
                     int mc_time);
    void Read_Parameters();
    void Read_Positions(box *Box, particles &Particles, int mc_time);
    void Read_Orientations(box *Box, particles &Particles, int mc_time);
    void Read_Orientations_EULER(box *Box, particles &Particles);
    void Read_Orientations_Axis_Angle(box *Box, particles &Particles,
                                      int mc_time);
    void Read_Box(box *Box, int mc_time);
    void Read_Random_State(gsl_rng *r, gsl_rng *r01, int mc_time);
    void Read_Patches(box *Box, particles &Particles, int mc_time);

    void Write_Acceptances(int mc_time, double accept_translate_procent,
                           double accept_rotate_procent,
                           double accept_iso_vol_procent,
                           double accept_complete_procent);
    void Write_NPT(int mc_time, box *Box);
    void Write_Energy(box *Box, particles &Particles, int mc_time);
};

#endif
