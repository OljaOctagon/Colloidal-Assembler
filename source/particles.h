#ifndef _particles_h_included_
#define _particles_h_included_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iomanip>
#include <iostream>
// include solid 3.5.6 library
#include <sstream>
#include <string>

using namespace std;

#include "cube.h"
#include "gear.h"
#include "hexbipyramid.h"
#include "octahedron.h"
#include "parallelepiped.h"
#include "polyhedra.h"
#include "rhombohedron.h"
#include "truncated_cube.h"
#include "rectangle.h"

//# include "cuboid.h"

#include "box.h"
//# include "box_non_ortho.h"
#include "alg_struct.h"
#include "cell.h"

class list_elements {

  public:
    list_elements();
    ~list_elements();

    m_vector distance;
    double distance_norm;
    bool is_neighbour;

    int nl_id;
};

class collision_list {

  public:
    collision_list();
    ~collision_list();
    void Set(int MAX_coll_p);

    list_elements *Elements;
    int Nm;  // Number of members of Collison List
    int NNm; // Number of Neighbours according to chosen Cut_Off;
    int cell_j;
    int coll_member_counter;
    int p_id, n_id;
    double particle_dist_x, particle_dist_y, particle_dist_z;
    double particle_dist_squared;
    int MAX_coll_p;

    // cubic cutoff
    cube Cutoff_box;
    double Scaleing;
    m_vector nbd;
    m_vector NC;
    double LXC;

    // void Calculate(box* Box, int id, int* Id_Cell_List, int** Cell_List,
    // cell* Cell, cube* N_Particle);

    void Calculate(box *Box, int id, int *Id_Cell_List, int **Cell_List,
                   cell *Cell, polyhedra **N_Particle, double Cut_Off,
                   int MAX_coll_p);

    void Calculate_Neighbours(double Cut_Off);
    void Calculate_Neighbours_Cubic(double Cut_Off, polyhedra **N_particle,
                                    int id);
    void Calculate_OP(box *Box, int id, polyhedra **N_Particle, double Cut_Off,
                      int MAX_coll_p);
    void Calculate_Bonds(box *Box, int id, polyhedra **N_Particle,
                         int MAX_coll_p);
};

class particles : public polyhedra {

  public:
    particles();
    particles(particles const &) {
        std::cout << "particles copy constructor called!!! This is bad!\n"
                  << std::flush;
    }
    particles(int number_of_cells_in, int size, int MAX_coll_p_in,
              int MAX_fshell_in, string particle_type, string binary_on_in, string ternary_on_in);
    ~particles();

    cube *Cube;
    cube *Cube_old;
    octahedron *Octahedron;
    octahedron *Octahedron_old;
    truncated_cube *Truncated_Cube;
    truncated_cube *Truncated_Cube_old;
    rhombohedron *Rhombohedron;
    rhombohedron *Rhombohedron_old;
    hexbipyramid *Hexbipyramid;
    hexbipyramid *Hexbipyramid_old;
    rectangle *Rectangle;
    rectangle *Rectangle_old;

    gear *Gear;
    gear *Gear_old;

    parallelepiped *Parallelepiped;
    parallelepiped *Parallelepiped_old;

    polyhedra **N_Particle;
    polyhedra **N_Particle_old;

    collision_list *Collision_List;
    collision_list *Collision_List_old;

    // neighbour_list* Neighbour_List;
    gsl_rng *randpatch;

    double Total_Energy;

    double phi_binary;
    string binary_on;
    string ternary_on;

    // double V_p;
    double N_p_float;
    int N_p;
    double l_diff;
    int id;
    int gid;
    double cos_a, cos_b, cos_c, sin_a, sin_b, sin_c;
    int Res_Size;

    int **Neighbour_List;
    int **Cell_List;
    int **Cell_List_old;
    int *Id_Cell_List;
    int *Id_Cell_List_old;

    // int** Collision_List;
    int *Bond_List;
    double inital_cell_volume;
    double V_fraction_Cell;
    int number_of_cells;
    int max_neighbours;
    int Id_Cell_List_test;
    double d_trans_vec;
    double particle_dist, particle_dist_x, particle_dist_y, particle_dist_z;
    double particle_dist_squared;
    int nx, ny, nz;
    int n_id;
    int Cell_x, Cell_y, Cell_z;
    int cell_id;
    cell *Cell;
    cell *Cell_old;

    int c_id, d_id;

    double N_c_float;
    int N_c;
    int Id_Cell_x, Id_Cell_y, Id_Cell_z;
    int scout;
    int id_num;
    double dist_trans_vec;
    int memb_n_id;
    int verlet_counter, cell_counter;
    int Cell_nx, Cell_ny, Cell_nz;
    double id_x_center, id_y_center, id_z_center;
    int MAX_bond_number;
    int b_count;
    int max_id;
    double Rot_mat_INIT[9];
    double Rot_mat_INIT_EULER[9];
    int MAX_coll_p;
    int MAX_fshell_p;
    int MAX_cell_members;

    // static const char* object2_filename;

    void Startconfig(box *Box);
    void Startconfig_former(box *Box);
    void Set_former_Config(box *Box);

    void Check_Periodic_BC(int id, box *Box);
    void Check_Periodic_CM(int id, box *Box);

    void Make_Cell_List(box *Box);
    void Update_Cell_List(int id, box *Box);
    void Update_Cell_List(int *List, int N_List, box *Box);


    void Check_Cell(int id, int cell_id, box *Box);
    void Make_Cell_Neighbour_List();


    // void Calculate_Collision_Partners(int id, box* Box);
    // void Calculate_Axis(int id);

    void Set_Cell_List(box *Box);
    void Set_Cell_List(int *List, int N_List, box *Box);
    void Set_Cell_List(box *Box, int id);


    void Reset_Cell_List(box *Box);
    void Reset_Cell_List(box *Box, int id);
    void Reset_Cell_List(int *List, int N_List, box *Box);

    void Startconfig_SPHERE(box *Box);

    bool ALLOCATED;
};

#endif
