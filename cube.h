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


double Lx;
double V;
int N_independent_faces;

m_vector* facenormal;

void edges_from_center();
void distance_from_center();
void Calculate_Axis(); 
void Set_Lengths();
void Set_Start_Lattice_Position(int id, double box_Lx, int N_box);
void Calculate_Face_Normals();
		
		
};
		
# endif
