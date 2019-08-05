#include "particles.h"

void particles::Check_Cell(int id, int cell_id, box *Box) {

    Cell[cell_id].left_count = 0;
    Cell[cell_id].right_count = 0;

    Cell[cell_id].front_count = 0;
    Cell[cell_id].back_count = 0;

    Cell[cell_id].top_count = 0;
    Cell[cell_id].bottom_count = 0;

    N_Particle[cell_id]->cell_out = 0;

    Cell[cell_id].edges_from_center();

    id_x_center = N_Particle[id]->x_center;
    id_y_center = N_Particle[id]->y_center;
    id_z_center = N_Particle[id]->z_center;

    if (id_x_center > Cell[cell_id].x[1]) {
        Cell[c_id].right_count = 1;
    }

    if (id_x_center < Cell[cell_id].x[0]) {
        Cell[c_id].left_count = 1;
    }

    if (id_y_center > Cell[cell_id].y[2]) {
        Cell[c_id].back_count = 1;
    }

    if (id_y_center < Cell[cell_id].y[0]) {
        Cell[c_id].front_count = 1;
    }

    if (id_z_center > Cell[cell_id].z[4]) {
        Cell[c_id].top_count = 1;
    }

    if (id_z_center < Cell[cell_id].z[0]) {
        Cell[c_id].bottom_count = 1;
    }

    N_Particle[id]->cell_out =
        Cell[cell_id].right_count + Cell[cell_id].left_count +
        Cell[cell_id].back_count + Cell[cell_id].front_count +
        Cell[cell_id].bottom_count + Cell[cell_id].top_count;
}

void particles::Make_Cell_List(box *Box) {

    ////3D//////
    // N_c = rint(pow(number_of_cells,1./3.));
    // N_c_float = double(N_c);

    /// 2D//////
    N_c = rint(pow(number_of_cells, 1. / 2.));
    N_c_float = double(N_c);

    for (int id = 0; id < number_of_cells; id++) {

        Cell[id].x_center =
            Box->x[0] + Cell[0].Lx / 2.0 + Cell[0].Lx * double(id % N_c);
        Cell[id].y_center = Box->y[0] + Cell[0].Ly / 2.0 +
                            Cell[0].Ly * double((id / N_c) % N_c);
        Cell[id].z_center = Box->z[0] + Cell[0].Lz / 2.0 +
                            Cell[0].Lz * double(id / (N_c * N_c));
        Cell[id].edges_from_center();
    }

    // Link id of particle to id of cell
    for (int id = 0; id < Box->N; id++) {

        Check_Periodic_CM(id, Box);

        Id_Cell_x =
            floor((N_Particle[id]->x_center - Box->x[0]) / double(Cell[0].Lx));
        Id_Cell_y =
            floor((N_Particle[id]->y_center - Box->y[0]) / double(Cell[0].Ly));
        Id_Cell_z =
            floor((N_Particle[id]->z_center - Box->z[0]) / double(Cell[0].Lz));

        Id_Cell_List[id] = Id_Cell_x + Id_Cell_y * N_c + Id_Cell_z * N_c * N_c;
    }

    // Create Cell Lists: assign to each cell all particles in it,
    // Cell_List[c_id][0] gives the number of particles in cell

    for (int c_id = 0; c_id < number_of_cells; c_id++) {

        // cout<<"MAX_cell_members"<<MAX_cell_members<<endl;

        for (int j = 1; j < MAX_cell_members; j++) {
            Cell_List[c_id][j] = -100;
        }

        Cell_List[c_id][0] = 0;
        cell_counter = 0;

        for (int id = 0; id < Box->N; id++) {

            if (Id_Cell_List[id] == c_id) {
                Cell_List[c_id][cell_counter + 1] = id;
                cell_counter = cell_counter + 1;
            }
        }

        Cell_List[c_id][0] = cell_counter;
    }
}

void particles::Make_Cell_Neighbour_List() {

    N_c = rint(pow(number_of_cells, 1. / 3.));
    N_c_float = double(N_c);

    // assign Neighbours (with periodic Boundary Conditions)
    for (int c_id = 0; c_id < number_of_cells; c_id++) {

        Cell[c_id].Cell_x = c_id % N_c;
        Cell[c_id].Cell_y = (c_id / N_c) % N_c;
        Cell[c_id].Cell_z = c_id / (N_c * N_c);

        for (int n_x = 0; n_x < 3; n_x++) {
            nx = n_x - 1;
            for (int n_y = 0; n_y < 3; n_y++) {
                ny = n_y - 1;
                for (int n_z = 0; n_z < 3; n_z++) {
                    nz = n_z - 1;
                    Cell_nx = (Cell[c_id].Cell_x + nx + N_c) % N_c;
                    Cell_ny = (Cell[c_id].Cell_y + ny + N_c) % N_c;
                    Cell_nz = (Cell[c_id].Cell_z + nz + N_c) % N_c;

                    Cell[c_id].Neighbour[n_x][n_y][n_z] =
                        Cell_nx + Cell_ny * N_c + Cell_nz * N_c * N_c;
                }
            }
        }
    }
}

void particles::Update_Cell_List(int id, box *Box) {

    c_id = Id_Cell_List[id];

    Id_Cell_List_old[id] = c_id;

    // Test if particle is outside its cell
    Check_Periodic_CM(id, Box);
    Check_Cell(id, c_id, Box);

    if (N_Particle[id]->cell_out >= 1) {

        // find cell id of new cell

        Cell_nx =
            floor((N_Particle[id]->x_center - Box->x[0]) / double(Cell[0].Lx));
        Cell_ny =
            floor((N_Particle[id]->y_center - Box->y[0]) / double(Cell[0].Ly));
        Cell_nz =
            floor((N_Particle[id]->z_center - Box->z[0]) / double(Cell[0].Lz));

        n_id = Cell_nx + Cell_ny * N_c + Cell_nz * N_c * N_c;

        // store old cell-lists
        for (int j = 0; j < MAX_cell_members; j++) {
            Cell_List_old[c_id][j] = Cell_List[c_id][j];
        }

        for (int j = 0; j < MAX_cell_members; j++) {
            Cell_List_old[n_id][j] = Cell_List[n_id][j];
        }

        // where is id in Cell List
        scout = 0;
        id_num = 1;

        do {

            if (Cell_List[c_id][id_num] == id) {
                scout = 1;
            }

            id_num = id_num + 1;

        } while (scout == 0);

        id_num = id_num - 1;

        // change cell counter to +/-
        Cell_List[c_id][0] = Cell_List[c_id][0] - 1;
        Cell_List[n_id][0] = Cell_List[n_id][0] + 1;

        cell_counter = Cell_List[c_id][0];
        Cell_List[c_id][cell_counter + 1] = -100;

        // kick id out of Cell_List[c_id]
        for (int i = id_num; i <= Cell_List[c_id][0]; i++) {
            Cell_List[c_id][i] = Cell_List_old[c_id][i + 1];
        }

        // add id at end of new Cell_List
        cell_counter = Cell_List[n_id][0];

        Cell_List[n_id][cell_counter] = id;

        Id_Cell_List[id] = n_id;
    }
}

void particles::Set_Cell_List(box *Box) {

    // Cell coordinates
    for (int id = 0; id < number_of_cells; id++) {

        Cell_old[id].x_center = Cell[id].x_center;
        Cell_old[id].y_center = Cell[id].y_center;
        Cell_old[id].z_center = Cell[id].z_center;

        for (int j = 0; j < Cell[id].edge_N; j++) {

            Cell_old[id].x[j] = Cell[id].x[j];
            Cell_old[id].y[j] = Cell[id].y[j];
            Cell_old[id].z[j] = Cell[id].z[j];
        }
    }

    // Link id of particle to id of cell
    for (int id = 0; id < Box->N; id++) {
        Id_Cell_List_old[id] = Id_Cell_List[id];
        // cout<<"Id_Cell_List["<<id<<"]:"<<Id_Cell_List[id]<<endl;
    }

    // Create Cell Lists: assign to each cell all particles in it,
    // Cell_List[c_id][0] gives the number of particles in cell
    for (int c_id = 0; c_id < number_of_cells; c_id++) {
        for (int j = 0; j < MAX_cell_members; j++) {
            Cell_List_old[c_id][j] = Cell_List[c_id][j];
        }
    }
}

void particles::Reset_Cell_List(box *Box) {

    // Link id of particle to id of cell
    for (int id = 0; id < Box->N; id++) {
        Id_Cell_List[id] = Id_Cell_List_old[id];
    }

    // Create Cell Lists: assign to each cell all particles in it,
    // Cell_List[c_id][0] gives the number of particles in cell
    for (int c_id = 0; c_id < number_of_cells; c_id++) {
        for (int j = 0; j < MAX_cell_members; j++) {
            Cell_List[c_id][j] = Cell_List_old[c_id][j];
        }
    }
}

void particles::Reset_Cell_List(box *Box, int id, int &c_id, int &n_id,
                                int &id_num) {

    if (N_Particle[id]->cell_out >= 1) {

        n_id = Id_Cell_List[id];
        c_id = Id_Cell_List_old[id];

        cell_counter = Cell_List[n_id][0];
        Cell_List[n_id][cell_counter] = -100;

        Cell_List[c_id][0] = Cell_List_old[c_id][0];
        Cell_List[n_id][0] = Cell_List_old[n_id][0];

        Id_Cell_List[id] = Id_Cell_List_old[id];

        for (int i = id_num; i < MAX_cell_members; i++) {
            Cell_List[c_id][i] = Cell_List_old[c_id][i];
        }
    }
}

list_elements::list_elements() {}

list_elements::~list_elements() {}

collision_list::collision_list() {}

void collision_list::Set(int MAX_coll_p) {

    Elements = new list_elements[MAX_coll_p];
    Nm = 0;

    for (int i = 0; i < MAX_coll_p; i++) {

        Elements[i].nl_id = -100;

        Elements[i].distance.x = 0;
        Elements[i].distance.y = 0;
        Elements[i].distance.z = 0;

        Elements[i].is_neighbour = 0;
    }
}

collision_list::~collision_list() { delete[] Elements; }

void collision_list::Calculate_Neighbours(double Cut_Off) {

    NNm = 0;

    for (int i = 0; i < Nm; i++) {

        Elements[i].is_neighbour = 0;

        if (Elements[i].distance_norm <= Cut_Off) {
            Elements[i].is_neighbour = 1;

            NNm = NNm + 1;
        }
    }
}

void collision_list::Calculate_Neighbours_Cubic(double Cut_Off,
                                                polyhedra **N_Particle,
                                                int id) {

    NNm = 0;

    Cutoff_box.Lx = Cut_Off;

    // Calculate eges for Cutoff_box

    Cutoff_box.x_center = N_Particle[id]->x_center;
    Cutoff_box.y_center = N_Particle[id]->y_center;
    Cutoff_box.z_center = N_Particle[id]->z_center;

    Scaleing = Cutoff_box.Lx / N_Particle[id]->Lx;

    for (int i = 0; i <= N_Particle[id]->edge_N; i++) {
        Cutoff_box.x[i] =
            (N_Particle[id]->x[i] - N_Particle[id]->x_center) * Scaleing;
    }

    // Calculate Axis of particle id
    N_Particle[id]->Calculate_Axis();

    for (int i = 0; i < Nm; i++) {

        Elements[i].is_neighbour = 0;

        nbd.x = Elements[i].distance.x;
        nbd.y = Elements[i].distance.y;
        nbd.z = Elements[i].distance.z;

        // calculate distance.x, distance.y, distance.z in body centered
        // coordinate system of particle id
        NC.x = N_Particle[id]->ax_1.x * nbd.x +
               N_Particle[id]->ax_2.x * nbd.y + N_Particle[id]->ax_3.x * nbd.z;
        NC.y = N_Particle[id]->ax_1.y * nbd.x +
               N_Particle[id]->ax_2.y * nbd.y + N_Particle[id]->ax_3.y * nbd.z;
        NC.z = N_Particle[id]->ax_1.z * nbd.x +
               N_Particle[id]->ax_2.z * nbd.y + N_Particle[id]->ax_3.z * nbd.z;

        LXC = Cutoff_box.Lx;

        if ((fabs(NC.x) <= LXC) && (fabs(NC.y) <= LXC) &&
            (fabs(NC.z) <= LXC)) {
            Elements[i].is_neighbour = 1;
            NNm = NNm + 1;
        }
    }
}

void collision_list::Calculate(box *Box, int id, int *Id_Cell_List,
                               int **Cell_List, cell *Cell,
                               polyhedra **N_Particle, double Cut_Off,
                               int MAX_coll_p) {

    n_id = Id_Cell_List[id];
    // cout<<"n_id "<<n_id<<endl;
    // reset Collision_List
    double Cut_Off_Squared;

    Cut_Off_Squared = Cut_Off * Cut_Off;
    // cout<<"Cut_Off_Squared: "<<Cut_Off_Squared<<endl;

    for (int i = 0; i < MAX_coll_p; i++) {
        Elements[i].nl_id = -100;
    }

    Nm = 0;
    coll_member_counter = 0;

    for (int n_x = 0; n_x < 3; n_x++) {
        for (int n_y = 0; n_y < 3; n_y++) {
            for (int n_z = 0; n_z < 3; n_z++) {

                // cout<<"Cell[n_id].Neighbour[n_x][n_y][n_z]:
                // "<<Cell[n_id].Neighbour[n_x][n_y][n_z]<<endl;
                p_id = Cell[n_id].Neighbour[n_x][n_y][n_z];

                for (int j = 1; j <= Cell_List[p_id][0]; j++) {

                    // cout<<"Cell_List["<<p_id<<"][0]"<<Cell_List[p_id][0]<<endl;

                    cell_j = Cell_List[p_id][j];

                    // for(int cell_j=0; cell_j<Box->N;cell_j++){

                    particle_dist_x = N_Particle[id]->x_center -
                                      N_Particle[cell_j]->x_center;
                    particle_dist_y = N_Particle[id]->y_center -
                                      N_Particle[cell_j]->y_center;
                    particle_dist_z = N_Particle[id]->z_center -
                                      N_Particle[cell_j]->z_center;

                    particle_dist_x =
                        particle_dist_x -
                        Box->Lx * rint(particle_dist_x / Box->Lx);
                    particle_dist_y =
                        particle_dist_y -
                        Box->Ly * rint(particle_dist_y / Box->Ly);
                    particle_dist_z =
                        particle_dist_z -
                        Box->Lz * rint(particle_dist_z / Box->Lz);

                    particle_dist_squared = particle_dist_x * particle_dist_x +
                                            particle_dist_y * particle_dist_y +
                                            particle_dist_z * particle_dist_z;

                    if (particle_dist_squared <= (Cut_Off_Squared) &&
                        cell_j != id) {

                        Elements[coll_member_counter].nl_id = cell_j;
                        Elements[coll_member_counter].distance.x =
                            particle_dist_x;
                        Elements[coll_member_counter].distance.y =
                            particle_dist_y;
                        Elements[coll_member_counter].distance.z =
                            particle_dist_z;
                        Elements[coll_member_counter].distance_norm =
                            Elements[coll_member_counter].distance.norm();

                        coll_member_counter = coll_member_counter + 1;
                    }

                    //}
                }
            }
        }
    }

    Nm = coll_member_counter;
}

void collision_list::Calculate_OP(box *Box, int id, polyhedra **N_Particle,
                                  double Cut_Off, int MAX_coll_p) {
    double Cut_Off_Squared;
    Cut_Off_Squared = Cut_Off * Cut_Off;
    Nm = 0;
    coll_member_counter = 0;

    for (int i = 0; i < MAX_coll_p; i++) {
        Elements[i].nl_id = -100;
    }

    for (cell_j = 0; cell_j < Box->N; cell_j++) {
        particle_dist_x =
            N_Particle[id]->x_center - N_Particle[cell_j]->x_center;
        particle_dist_y =
            N_Particle[id]->y_center - N_Particle[cell_j]->y_center;
        particle_dist_z =
            N_Particle[id]->z_center - N_Particle[cell_j]->z_center;

        particle_dist_x =
            particle_dist_x - Box->Lx * rint(particle_dist_x / Box->Lx);
        particle_dist_y =
            particle_dist_y - Box->Ly * rint(particle_dist_y / Box->Ly);
        particle_dist_z =
            particle_dist_z - Box->Lz * rint(particle_dist_z / Box->Lz);

        particle_dist_squared = particle_dist_x * particle_dist_x +
                                particle_dist_y * particle_dist_y +
                                particle_dist_z * particle_dist_z;

        if (particle_dist_squared <= (Cut_Off_Squared) && cell_j != id) {

            Elements[coll_member_counter].nl_id = cell_j;
            Elements[coll_member_counter].distance.x = particle_dist_x;
            Elements[coll_member_counter].distance.y = particle_dist_y;
            Elements[coll_member_counter].distance.z = particle_dist_z;
            Elements[coll_member_counter].distance_norm =
                Elements[coll_member_counter].distance.norm();

            coll_member_counter = coll_member_counter + 1;
        }
    }

    Nm = coll_member_counter;
}

void collision_list::Calculate_Bonds(box *Box, int id1, polyhedra **N_Particle,
                                     int MAX_coll_p) {

    Nm = 0;
    double patch_distance_squared;
    coll_member_counter = 0;

    for (int i = 0; i < MAX_coll_p; i++) {
        Elements[i].nl_id = -100;
    }

    for (int id2 = 0; id2 < Box->N; id2++) {

        for (int pid1 = 0; pid1 < N_Particle[id1]->N_patches; pid1++) {
            for (int pid2 = 0; pid2 < N_Particle[id1]->N_patches; pid2++) {

                particle_dist_x = N_Particle[id1]->x_patch[pid1] -
                                  N_Particle[id2]->x_patch[pid2];
                particle_dist_y = N_Particle[id1]->y_patch[pid1] -
                                  N_Particle[id2]->y_patch[pid2];
                particle_dist_z = N_Particle[id1]->z_patch[pid1] -
                                  N_Particle[id2]->z_patch[pid2];

                particle_dist_x = particle_dist_x -
                                  Box->Lx * rint(particle_dist_x / Box->Lx);
                particle_dist_y = particle_dist_y -
                                  Box->Ly * rint(particle_dist_y / Box->Ly);
                particle_dist_z = particle_dist_z -
                                  Box->Lz * rint(particle_dist_z / Box->Lz);

                patch_distance_squared = particle_dist_x * particle_dist_x +
                                         particle_dist_y * particle_dist_y +
                                         particle_dist_z * particle_dist_z;

                int type_1 = N_Particle[id1]->patch_type[pid1];
                int type_2 = N_Particle[id2]->patch_type[pid2];

                if (patch_distance_squared <
                        N_Particle[id1]->patch_cutoff_squared[pid1] &&
                    id1 != id2 &&
                    N_Particle[id1]->patch_energy[type_1][type_2] != 0) {

                    Elements[coll_member_counter].nl_id = id2;
                    Elements[coll_member_counter].distance.x = particle_dist_x;
                    Elements[coll_member_counter].distance.y = particle_dist_y;
                    Elements[coll_member_counter].distance.z = particle_dist_z;
                    Elements[coll_member_counter].distance_norm =
                        Elements[coll_member_counter].distance.norm();

                    coll_member_counter = coll_member_counter + 1;

                    }
                }
            }
        }
    }

    Nm = coll_member_counter;
}
