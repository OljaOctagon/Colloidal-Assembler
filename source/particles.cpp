#include "particles.h"
#include <iostream>

particles::~particles() {

    if (ALLOCATED) {
        delete[] N_Particle;
        delete[] Collision_List;
        delete[] Collision_List_old;
        ALLOCATED = false;
        // std::cout << "DEALLOCATED" << "\n" << std::flush;
    }

    // std::cout << "particles destructor called!!!\n";
}

particles::particles() {

    // N_Particle = 0;
    // N_Particle_old = 0;
    Cell = 0;
    Cell_old = 0;

    ALLOCATED = false;

    // std::cout << "particles empty constructor called!!!\n";
}

particles::particles(int number_of_cells_in, int size, int MAX_coll_p_in,
                     int MAX_fshell_in, string particle_type,
                     string binary_on_in, double phi_binary_in) {
    ALLOCATED = true;

    // std::cout << "particles custom constructor called!!!\n";

    randpatch = gsl_rng_alloc(gsl_rng_ranlxd2);

    number_of_cells = number_of_cells_in;
    binary_on = binary_on_in;
    phi_binary = phi_binary_in;

    Res_Size = 2500;
    size = Res_Size;

    max_id = size - 1;
    MAX_coll_p = MAX_coll_p_in;
    MAX_fshell_p = MAX_fshell_p;

    N_Particle = new polyhedra *[size];
    N_Particle_old = new polyhedra *[size];

    if (particle_type.compare("cube") == 0) {
        Cube = new cube[size];
        Cube_old = new cube[size];

        for (int i = 0; i < size; i++) {

            N_Particle[i] = &Cube[i];
            N_Particle_old[i] = &Cube_old[i];
        }
    }

    if (particle_type.compare("octahedron") == 0) {

        Octahedron = new octahedron[size];
        Octahedron_old = new octahedron[size];

        for (int i = 0; i < size; i++) {

            N_Particle[i] = &Octahedron[i];
            N_Particle_old[i] = &Octahedron_old[i];
        }
    }

    if (particle_type.compare("truncated_cube") == 0) {

        Truncated_Cube = new truncated_cube[size];
        Truncated_Cube_old = new truncated_cube[size];

        for (int i = 0; i < size; i++) {

            N_Particle[i] = &Truncated_Cube[i];
            N_Particle_old[i] = &Truncated_Cube_old[i];
        }
    }


    if (particle_type.compare("rhombohedron") == 0) {

      Rhombohedron = new rhombohedron[size];
      Rhombohedron_old = new rhombohedron[size];

      for (int i = 0; i < size; i++) {

        N_Particle[i] = &Rhombohedron[i];
        N_Particle_old[i] = &Rhombohedron_old[i];
      }
    }

    if (particle_type.compare("rectangle") == 0) {

        Rectangle =  new rectangle[size];
        Rectangle_old = new rectangle[size];

        for (int i = 0; i < size; i++) {

            N_Particle[i] = &Rectangle[i];
            N_Particle_old[i] = &Rectangle_old[i];
        }
    }

    if (particle_type.compare("hexbipyramid") == 0) {

        Hexbipyramid = new hexbipyramid[size];
        Hexbipyramid_old = new hexbipyramid[size];

        for (int i = 0; i < size; i++) {

            N_Particle[i] = &Hexbipyramid[i];
            N_Particle_old[i] = &Hexbipyramid_old[i];
        }
    }

    if (particle_type.compare("parallelepiped") == 0) {

        Parallelepiped = new parallelepiped[size];
        Parallelepiped_old = new parallelepiped[size];

        for (int i = 0; i < size; i++) {

            N_Particle[i] = &Parallelepiped[i];
            N_Particle_old[i] = &Parallelepiped_old[i];
        }
    }

    if (particle_type.compare("gear") == 0) {

        Gear = new gear[size];
        Gear_old = new gear[size];

        for (int i = 0; i < size; i++) {

            N_Particle[i] = &Gear[i];
            N_Particle_old[i] = &Gear_old[i];
        }
    }

    Cell = new cell[number_of_cells];
    Cell_old = new cell[number_of_cells];

    if (particle_type.compare("gear") == 0) {

        Gear = new gear[size];
        Gear_old = new gear[size];

        for (int i = 0; i < size; i++) {

            N_Particle[i] = &Gear[i];
            N_Particle_old[i] = &Gear_old[i];
        }
    }

    // Initialization of N_Particle
    double rand_s;
    int a, b, c, d;

    for (int id = 0; id < size; id++) {
        N_Particle[id]->Set_Initial_Quaternion();
        N_Particle_old[id]->Set_Initial_Quaternion();

        N_Particle[id]->Set_Initial_Axis();
        N_Particle_old[id]->Set_Initial_Axis();

        N_Particle[id]->Set_Lengths();
        N_Particle_old[id]->Set_Lengths();

        // choose left or right-handedness.

        if (binary_on.compare("on") == 0) {
            rand_s = gsl_rng_uniform(randpatch);
            // binary L's
            if (rand_s < phi_binary) {
                a = 0;
                b = 2;
                c = 1;
                d = 1;
            }

            if (rand_s >= phi_binary) {
                a = 2;
                b = 0;
                c = 1;
                d = 1;
            }

            N_Particle[id]->Set_Lengths(a, b, c, d);
            N_Particle_old[id]->Set_Lengths(a, b, c, d);
        }
    }
}

void particles::Startconfig(box *Box) {

    // 3D///
    // N_p = rint(pow(double(Box->N),1./3.));
    // N_p_float = double(N_p);
    // 2D///
    N_p = rint(pow(double(Box->N), 1. / 2.));
    N_p_float = double(N_p);

    l_diff = (Box->Lx - N_p_float) / N_p_float;
    Box->N = int(Box->N);

    Cell[0].V = Box->V / double(number_of_cells);

    Cell[0].Lx = pow(Cell[0].V, 1. / 3.);
    Cell[0].Ly = Cell[0].Lx;
    Cell[0].Lz = Cell[0].Lx;

    Cell_old[0].V = Box->V / double(number_of_cells);
    Cell_old[0].Lx = pow(Cell[0].V, 1. / 3.);
    Cell_old[0].Ly = Cell[0].Lx;
    Cell_old[0].Lz = Cell[0].Lx;

    for (int c_id = 1; c_id < number_of_cells; c_id++) {

        Cell[c_id].V = Cell[0].V;
        Cell[c_id].Lx = Cell[0].Lx;
        Cell[c_id].Ly = Cell[0].Ly;
        Cell[c_id].Lz = Cell[0].Lz;

        Cell_old[c_id].V = Cell[0].V;
        Cell_old[c_id].Lx = Cell[0].Lx;
        Cell_old[c_id].Ly = Cell[0].Ly;
        Cell_old[c_id].Lz = Cell[0].Lz;

        Cell[c_id].left_count = 0;
        Cell[c_id].right_count = 0;

        Cell[c_id].front_count = 0;
        Cell[c_id].back_count = 0;

        Cell[c_id].top_count = 0;
        Cell[c_id].bottom_count = 0;
    }

    MAX_cell_members = 5 * int(ceil(Cell[0].V) / double(N_Particle[0]->V));

    int size;

    size = Res_Size;

    Collision_List = new collision_list[size];
    for (int id = 0; id < size; id++) {
        Collision_List[id].Set(MAX_coll_p);
    }

    Collision_List_old = new collision_list[size];
    for (int id = 0; id < size; id++) {
        Collision_List_old[id].Set(MAX_coll_p);
    }

    Cell_List = new int *[number_of_cells];
    for (int j = 0; j < number_of_cells; j++) {
        Cell_List[j] = new int[MAX_cell_members];
    }

    Cell_List_old = new int *[number_of_cells];
    for (int j = 0; j < number_of_cells; j++) {
        Cell_List_old[j] = new int[MAX_cell_members];
    }

    Id_Cell_List = new int[size];

    Id_Cell_List_old = new int[size];

    for (int id = 0; id < Box->N; id++) {

        N_Particle[id]->Set_Start_Lattice_Position(id, Box->Lx, Box->Ly,
                                                   Box->Lz, Box->N);

        /*
        N_Particle[id]->x_center = 0.5 +  l_diff/2.0 + double(id %N_p +
        l_diff*(id % N_p)); N_Particle[id]->y_center = 0.5 +  l_diff/2.0 +
        double((id/N_p)%N_p + l_diff*((id/N_p)%N_p)); N_Particle[id]->z_center
        = 0.5 +  l_diff/2.0 + double(id/int(pow((double)N_p,2)) +
        l_diff*(id/int((double)pow(N_p,2))));
    */

        N_Particle_old[id]->x_center = N_Particle[id]->x_center;
        N_Particle_old[id]->y_center = N_Particle[id]->y_center;
        N_Particle_old[id]->z_center = N_Particle[id]->z_center;

        N_Particle[id]->copy_count = 0;

        N_Particle[id]->edges_from_center();
        N_Particle_old[id]->edges_from_center();

        N_Particle[id]->trans_periodic[0] = 0.0;
        N_Particle[id]->trans_periodic[1] = 0.0;
        N_Particle[id]->trans_periodic[2] = 0.0;

        // calculate ax vectors

        N_Particle[id]->Calculate_Axis();
        N_Particle[id]->Calculate_Patch_Position();

        N_Particle[id]->cm_left_count = 0;
        N_Particle[id]->cm_right_count = 0;

        N_Particle[id]->cm_top_count = 0;
        N_Particle[id]->cm_bottom_count = 0;

        N_Particle[id]->cm_back_count = 0;
        N_Particle[id]->cm_front_count = 0;

        N_Particle[id]->cm_out = 0;

        N_Particle[id]->q.x = 0.0;
        N_Particle[id]->q.y = 0.0;
        N_Particle[id]->q.z = 0.0;
        N_Particle[id]->q.w = 1.0;
        N_Particle[id]->phi = 0.0;

        N_Particle_old[id]->q.x = N_Particle[id]->q.x;
        N_Particle_old[id]->q.y = N_Particle[id]->q.y;
        N_Particle_old[id]->q.z = N_Particle[id]->q.z;
        N_Particle_old[id]->q.w = N_Particle[id]->q.w;
        N_Particle[id]->phi_old = N_Particle[id]->phi;
    }

    // Make_Cell_List(Box, Cell, N_Particle);
    // Set_Cell_List(Box, Cell, Cell_old, N_Particle);
    // Make_Cell_Neighbour_List(Cell);

    // for Resorvoir

    for (int id = Box->N; id < Res_Size; id++) {

        N_Particle[id]->x_center = -10.0;
        N_Particle[id]->y_center = -10.0;
        N_Particle[id]->z_center = -10.0;

        N_Particle[id]->phi = 0.0;
        N_Particle[id]->edges_from_center();
        N_Particle[id]->Calculate_Axis();
        N_Particle[id]->Calculate_Patch_Position();

        N_Particle_old[id]->x_center = N_Particle[id]->x_center;
        N_Particle_old[id]->y_center = N_Particle[id]->y_center;
        N_Particle_old[id]->z_center = N_Particle[id]->z_center;

        for (int j = 0; j < N_Particle[id]->edge_N; j++) {

            N_Particle_old[id]->x[j] = N_Particle[id]->x[j];
            N_Particle_old[id]->y[j] = N_Particle[id]->y[j];
            N_Particle_old[id]->z[j] = N_Particle[id]->z[j];
        }

        N_Particle_old[id]->q.x = N_Particle[id]->q.x;
        N_Particle_old[id]->q.y = N_Particle[id]->q.y;
        N_Particle_old[id]->q.z = N_Particle[id]->q.z;
        N_Particle_old[id]->q.w = N_Particle[id]->q.w;
        N_Particle[id]->phi_old = N_Particle[id]->phi;
    }
}

void particles::Set_former_Config(box *Box) {

    double phi_t;

    for (int id = 0; id < Box->N; id++) {

        phi_t = N_Particle[id]->phi;

        Rot_mat_INIT[0] = cos(phi_t);
        Rot_mat_INIT[1] = -sin(phi_t);
        Rot_mat_INIT[2] = 0;

        Rot_mat_INIT[3] = sin(phi_t);
        Rot_mat_INIT[4] = cos(phi_t);
        Rot_mat_INIT[5] = 0;

        Rot_mat_INIT[6] = 0;
        Rot_mat_INIT[7] = 0;
        Rot_mat_INIT[8] = 1;

        /*
         Rot_mat_INIT[0] = 1.0 - 2.0*N_Particle[id]->q.y*N_Particle[id]->q.y
         - 2.0*N_Particle[id]->q.z*N_Particle[id]->q.z; Rot_mat_INIT[1]
         = 2.0*N_Particle[id]->q.x*N_Particle[id]->q.y
         - 2.0*N_Particle[id]->q.w*N_Particle[id]->q.z; Rot_mat_INIT[2]
         = 2.0*N_Particle[id]->q.x*N_Particle[id]->q.z
         + 2.0*N_Particle[id]->q.w*N_Particle[id]->q.y;

         Rot_mat_INIT[3] = 2.0*N_Particle[id]->q.x*N_Particle[id]->q.y
         + 2.0*N_Particle[id]->q.w*N_Particle[id]->q.z; Rot_mat_INIT[4] = 1.0
         - 2.0*N_Particle[id]->q.x*N_Particle[id]->q.x
         - 2.0*N_Particle[id]->q.z*N_Particle[id]->q.z; Rot_mat_INIT[5]
         = 2.0*N_Particle[id]->q.y*N_Particle[id]->q.z
         - 2.0*N_Particle[id]->q.w*N_Particle[id]->q.x;

         Rot_mat_INIT[6] = 2.0*N_Particle[id]->q.x*N_Particle[id]->q.z
         - 2.0*N_Particle[id]->q.w*N_Particle[id]->q.y; Rot_mat_INIT[7]
         = 2.0*N_Particle[id]->q.y*N_Particle[id]->q.z
         + 2.0*N_Particle[id]->q.w*N_Particle[id]->q.x; Rot_mat_INIT[8] = 1.0
         - 2.0*N_Particle[id]->q.x*N_Particle[id]->q.x
         - 2.0*N_Particle[id]->q.y*N_Particle[id]->q.y;
        */

        N_Particle[id]->edges_from_center();
        N_Particle[id]->distance_from_center();

        for (int j = 0; j < N_Particle[id]->edge_N; j++) {

            N_Particle[id]->new_dist_x[j] =
                Rot_mat_INIT[0] * N_Particle[id]->dist_x[j] +
                Rot_mat_INIT[1] * N_Particle[id]->dist_y[j] +
                Rot_mat_INIT[2] * N_Particle[id]->dist_z[j];
            N_Particle[id]->new_dist_y[j] =
                Rot_mat_INIT[3] * N_Particle[id]->dist_x[j] +
                Rot_mat_INIT[4] * N_Particle[id]->dist_y[j] +
                Rot_mat_INIT[5] * N_Particle[id]->dist_z[j];
            N_Particle[id]->new_dist_z[j] =
                Rot_mat_INIT[6] * N_Particle[id]->dist_x[j] +
                Rot_mat_INIT[7] * N_Particle[id]->dist_y[j] +
                Rot_mat_INIT[8] * N_Particle[id]->dist_z[j];

            N_Particle[id]->x[j] =
                N_Particle[id]->x_center + N_Particle[id]->new_dist_x[j];
            N_Particle[id]->y[j] =
                N_Particle[id]->y_center + N_Particle[id]->new_dist_y[j];
            N_Particle[id]->z[j] =
                N_Particle[id]->z_center + N_Particle[id]->new_dist_z[j];
        }

        N_Particle_old[id]->x_center = N_Particle[id]->x_center;
        N_Particle_old[id]->y_center = N_Particle[id]->y_center;
        N_Particle_old[id]->z_center = N_Particle[id]->z_center;

        for (int j = 0; j < N_Particle[id]->edge_N; j++) {

            N_Particle_old[id]->x[j] = N_Particle[id]->x[j];
            N_Particle_old[id]->y[j] = N_Particle[id]->y[j];
            N_Particle_old[id]->z[j] = N_Particle[id]->z[j];
        }

        N_Particle_old[id]->q.x = N_Particle[id]->q.x;
        N_Particle_old[id]->q.y = N_Particle[id]->q.y;
        N_Particle_old[id]->q.z = N_Particle[id]->q.z;
        N_Particle_old[id]->q.w = N_Particle[id]->q.w;
        N_Particle[id]->phi_old = N_Particle[id]->phi;


        N_Particle_old[id]->patch_type[0] = N_Particle[id]->patch_type[0];
        N_Particle_old[id]->patch_type[1] = N_Particle[id]->patch_type[1];
        N_Particle_old[id]->patch_type[2] = N_Particle[id]->patch_type[2];
        N_Particle_old[id]->patch_type[3] = N_Particle[id]->patch_type[3];

        N_Particle[id]->Calculate_Axis();
        N_Particle[id]->Calculate_Patch_Position();
    }
}

void particles::Check_Periodic_CM(int id, box *Box) {

    N_Particle[id]->cm_left_count = 0;
    N_Particle[id]->cm_right_count = 0;
    N_Particle[id]->cm_front_count = 0;
    N_Particle[id]->cm_back_count = 0;
    N_Particle[id]->cm_top_count = 0;
    N_Particle[id]->cm_bottom_count = 0;

    N_Particle[id]->cm_out = 0;

    if (N_Particle[id]->x_center > Box->x[1]) {
        N_Particle[id]->cm_right_count = 1;
    }

    if (N_Particle[id]->x_center < Box->x[0]) {
        N_Particle[id]->cm_left_count = 1;
    }

    if (N_Particle[id]->y_center > Box->y[3]) {
        N_Particle[id]->cm_back_count = 1;
    }

    if (N_Particle[id]->y_center < Box->y[0]) {
        N_Particle[id]->cm_front_count = 1;
    }

    if (N_Particle[id]->z_center > Box->z[4]) {
        N_Particle[id]->cm_top_count = 1;
    }

    if (N_Particle[id]->z_center < Box->z[0]) {
        N_Particle[id]->cm_bottom_count = 1;
    }

    N_Particle[id]->cm_out =
        N_Particle[id]->cm_right_count + N_Particle[id]->cm_left_count +
        N_Particle[id]->cm_back_count + N_Particle[id]->cm_front_count +
        N_Particle[id]->cm_top_count + N_Particle[id]->cm_bottom_count;
}

void particles::Check_Periodic_BC(int id, box *Box) {

    N_Particle[id]->left_count = 0;
    N_Particle[id]->right_count = 0;
    N_Particle[id]->front_count = 0;
    N_Particle[id]->back_count = 0;
    N_Particle[id]->top_count = 0;
    N_Particle[id]->bottom_count = 0;

    N_Particle[id]->sum_edge_out = 0;

    for (int k = 0; k < N_Particle[id]->edge_N; k++) {
        N_Particle[id]->edge_out[k] = 0;
    }

    for (int k = 0; k < N_Particle[id]->edge_N; k++) {

        if (N_Particle[id]->x[k] > Box->x[1]) {
            N_Particle[id]->right_count = 1;
            N_Particle[id]->edge_out[k] = 1;
        }

        if (N_Particle[id]->x[k] < Box->x[0]) {
            N_Particle[id]->left_count = 1;
            N_Particle[id]->edge_out[k] = 1;
        }

        if (N_Particle[id]->y[k] > Box->y[2]) {
            N_Particle[id]->back_count = 1;
            N_Particle[id]->edge_out[k] = 1;
        }

        if (N_Particle[id]->y[k] < Box->y[0]) {
            N_Particle[id]->front_count = 1;
            N_Particle[id]->edge_out[k] = 1;
        }

        if (N_Particle[id]->z[k] > Box->z[4]) {
            N_Particle[id]->top_count = 1;
            N_Particle[id]->edge_out[k] = 1;
        }

        if (N_Particle[id]->z[k] < Box->z[0]) {
            N_Particle[id]->bottom_count = 1;
            N_Particle[id]->edge_out[k] = 1;
        }

        N_Particle[id]->sum_edge_out =
            N_Particle[id]->sum_edge_out + N_Particle[id]->edge_out[k];
    }
}
