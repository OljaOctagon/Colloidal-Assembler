# ifndef _octahedron_h_included_
# define _octahedron_h_included_



#define _USE_MATH_DEFINES	

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>


#include "polyhedra.h"

using namespace std;



class octahedron : public polyhedra
{ 


public:  

octahedron();
~ octahedron();


double Lx;

double V;
int  edge_N_vec;
double height;

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
