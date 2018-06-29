

#define _USE_MATH_DEFINES	


#include <ctime>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <sstream>
#include <complex>
#include <cmath>
#include <fstream>
#include <cstring>
#include <gsl/gsl_rng.h>

# include "particles.h"
# include "box.h"
//# include "box_non_ortho.h"
# include "move.h"
# include "fileio.h"
# include "polyhedra.h"
# include "order_parameter_2.h"



using namespace std;

double D2_order;
double D4_order;
  



char runtype[128];
char new_runtype[128];
char former_runtype[128];

int term_time;
int calc_frequency;
int frame_frequency;
int checkpoint_frequency;
int exit_code;

gsl_rng * r;
gsl_rng * r01;

int checkpoint_time;
int start_time;

	
