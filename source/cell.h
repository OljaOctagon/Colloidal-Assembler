
# ifndef _cell_h_included_
# define _cell_h_included_

#define _USE_MATH_DEFINES	

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "cuboid.h"
#include "box.h"
//#include "box_non_ortho.h"
#include "cube.h"


using namespace std;

class cell : public cuboid {

public:

int left_count, right_count, top_count, bottom_count, front_count, back_count;
int Cell_x,Cell_y,Cell_z;
int*** Neighbour;


int outside;
cell();
~cell();


};

#endif

