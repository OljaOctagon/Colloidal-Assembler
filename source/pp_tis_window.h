
#ifndef _tis_window_h_included_
#define _tis_window_h_included_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

class tis_window {

  public:
    // tis_window();
    //~tis_window();

    int w_number;
    int wv_current;
    int wv_0;
    int wv_ip1;
    int wv_i;

    int w_term_time;

    int is_proper_path;
    int is_lambda_0_path;
    int is_lambda_p1_path;
    int exceded_max_path_length;
    int comming_from_basinA;
    int has_crossed_wv_im1;
    int time_tau_b;
    int time_tau_f;
    int time_N;
    int has_crossed_wv_0;
    int no_crossing_with_lambda_i;
    int no_crossing_with_lambda_0;
    // int is_lambda_i_trajectory;

    int N_positive_crossing;
    int N_check_points;
    int check_point_frequency;
    string tis_string;
    string string_cp;

    bool exit_code;

    // gsl_histogram * histogram;
    size_t N_bins;
    gsl_histogram *histogram;

    // random control

    int N_random_list;
    int random_list[100];
    int read_state;
    int random_counter;
    string string_random1;
    string string_random2;
    string string_m;
    string string_bin;
    string string_3;
    int tis_time;

    int update_state;
    int set_state;
    int save_state;
    int close_state;
    int has_crossed_wv_i;
    double r_size;
    double xmin, xmax;
    double delta_lambda;
    double Flux;
    int traj_info;
    int wv_max;
    int has_crossed_NOW;
    int crossing_time;

    int comming_from_inside;
    int comming_from_outside;

    char traj_type[128];

    tis_window();

    tis_window(int wv_0_in, int wv_i_in, int wv_ip1_in, int term_time,
               char traj_type_in[128]);
    void Get_path_info(int wv, int mc_time);
    void Get_path_info_BACKWARDS(int wv, int mc_time);

    void Write_path_info(int N_check_points, int check_point_frequency,
                         int mc_time);
    bool tis_controller(int wv, int mc_time, int N_check_points,
                        int check_point_frequency);

    void Get_FLUX_path_info(int wv);
    void Get_FLUX_path_info_BACKWARDS(int wv);

    void Write_FLUX_path_info(int mc_time, int delta_t);
    bool FLUX_controller(int wv, int mc_time, int delta_t);

    void set_histogram(int N_bins, int xmin, int xmax);
    void unset_histogram();
    void update_histogram(int wv);
    void write_histogram_info();

    void Get_Ptraj_info(int wv, int mc_time);
    void Get_Ptraj_info_BACKWARDS(int wv, int mc_time);

    void Read_Ptraj_info(int checkpoint_time);
    void Write_Ptraj_info();
    void Ptraj_controller(int wv, int mc_time);

    void random_init();
    void random_control(gsl_rng *r, gsl_rng *r01, int time);
};

#endif
