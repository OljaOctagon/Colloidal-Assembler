#define _USE_MATH_DEFINES	

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>


#include "polyhedra.h"

using namespace std;



class parallelepiped : public polyhedra
{ 

public:  

double alpha;
double a_x;
double h;
double h_2;
double diag2_short, diag2_long;
double diag3_short, diag3_long;
double p,q;
double px_l, py_l, px_s, py_s;
double qx_l, qy_l, qx_s, qy_s;
double alpha_2, gamma_2;


parallelepiped();
~ parallelepiped();


void edges_from_center();
void distance_from_center();
void Calculate_Axis(); 
void Set_Axis();
void Set_Lengths();
void Set_Lengths(double Lx, double Ly, double Lz, double alpha, double gamma);
void Set_Start_Lattice_Position(int id, double box_Lx, int N_box);
void Calculate_Face_Normals();
double Calculate_Projection_to_Separating_Axis(m_vector laxis);
		
		
};
		
# endif
