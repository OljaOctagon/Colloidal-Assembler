# ifndef _hexbipyramid_h_included_
# define _hexbipyramid_h_included_



#define _USE_MATH_DEFINES	

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>


#include "polyhedra.h"

using namespace std;



class hexbipyramid : public polyhedra
{ 


public:  

hexbipyramid();
~ hexbipyramid();


double height;
double height_t;
double height_o;
double Lx_t;
double Vt;
double h_Lx;
double h_Lx_t;
double cf;
	    
void write_positions(ofstream& fout);	
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

