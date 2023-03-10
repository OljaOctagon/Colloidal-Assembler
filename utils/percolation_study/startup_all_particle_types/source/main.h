
#define _USE_MATH_DEFINES

#include <cmath>
#include <complex>
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
#include "particles.h"
//#include "box_non_ortho.h"
#include "fileio.h"
#include "order_parameter.h"
#include "pmove.h"
#include "polyhedra.h"
#include "tis_window.h"

using namespace std;

double D2_order;
double D4_order;

char runtype[128];
char new_runtype[128];
char former_runtype[128];
char tis_runtype[128];
char flux_runtype[128];

int term_time;
int calc_frequency;
int frame_frequency;
int checkpoint_frequency;
int N_check_points;
// int check_point_frequency;

char traj_type_in[128];

gsl_rng *r;
gsl_rng *r01;

int checkpoint_time;
int start_time;

bool exit_code;
bool term_code;

int wv_0_in, wv_ip1_in, wv_current_in, wv_i_in, wv_im1_in;
int current_window_value;
int wn_in;
int seed_r, seed_r01;
