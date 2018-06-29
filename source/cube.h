# ifndef _cube_h_included_
# define _cube_h_included_



#define _USE_MATH_DEFINES	

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>


#include "polyhedra.h"

using namespace std;



class cube : public polyhedra
{ 


public:  

cube();
~ cube();


m_vector* rel_pos;

void write_positions(ofstream& fout);
void edges_from_center();
void distance_from_center();
void Calculate_Axis(); 
void Set_Axis();
void Set_Lengths();
void Set_Start_Lattice_Position(int id, double box_Lx, int N_box);
void Set_Start_Lattice_Position(int id, double box_Lx, double box_Ly, double box_Lz, int N_box);


void Set_Lengths(double a);

void relative_distances();		
void Calculate_Face_Normals();
double Calculate_Projection_to_Separating_Axis(m_vector laxis);
		
		
};
		
# endif
