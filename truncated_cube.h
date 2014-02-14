# ifndef _truncated_cube_h_included_
# define _truncated_cube_h_included_



#define _USE_MATH_DEFINES	

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>


#include "polyhedra.h"

using namespace std;



class truncated_cube : public polyhedra
{ 


public:  

truncated_cube();
~ truncated_cube();


double Lx;
double edge_R;
double V;
int  edge_N_vec;

double lx,ly,lz;

int N_independent_faces;
int N_cross_edges;
m_vector* facenormal;
m_vector* edges;


void edges_from_center();
void distance_from_center();
void Calculate_Axis(); 
void Set_Axis();

void Set_Start_Lattice_Position(int id, double box_Lx, int N_box);
void Set_Lengths();	
double Calculate_Projection_to_Separating_Axis(m_vector laxis);

void Calculate_Face_Normals();  

};



# endif
