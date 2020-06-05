
#include "pmove.h"

double pmove::Calculate_Cross_Product_x(m_vector a, m_vector b) {
    double cx;
    cx = a.y * b.z - a.z * b.y;
    return cx;
}

double pmove::Calculate_Cross_Product_y(m_vector a, m_vector b) {
    double cy;
    cy = a.z * b.x - a.x * b.z;
    return cy;
}

double pmove::Calculate_Cross_Product_z(m_vector a, m_vector b) {
    double cz;
    cz = a.x * b.y - a.y * b.x;
    return cz;
}

int pmove::Calculate_Separating_Axis_GENERAL(particles &Particles, int id,
                                             int j, m_vector *L_axis) {

    int N_indep, N_cross;
    double norm_L;
    int N_separating_axis;

    N_indep = Particles.N_Particle[id]->N_independent_faces;
    N_cross = Particles.N_Particle[id]->N_cross_edges;

    N_separating_axis = 2 * N_indep + N_cross * N_cross;
    // cout<<"N_separating_axis: "<<N_separating_axis<<endl;

    for (int k = 0; k < N_indep; k++) {

        L_axis[k].x = Particles.N_Particle[id]->facenormal[k].x;
        L_axis[k].y = Particles.N_Particle[id]->facenormal[k].y;
        L_axis[k].z = Particles.N_Particle[id]->facenormal[k].z;
    }

    for (int k = N_indep; k < (2 * N_indep); k++) {

        L_axis[k].x = Particles.N_Particle[j]->facenormal[k - N_indep].x;
        L_axis[k].y = Particles.N_Particle[j]->facenormal[k - N_indep].y;
        L_axis[k].z = Particles.N_Particle[j]->facenormal[k - N_indep].z;
    }

    int k;

    k = 2 * N_indep;

    for (int s = 0; s < N_cross; s++) {
        for (int m = 0; m < N_cross; m++) {

            L_axis[k].x =
                Calculate_Cross_Product_x(Particles.N_Particle[id]->edges[s],
                                          Particles.N_Particle[j]->edges[m]);
            L_axis[k].y =
                Calculate_Cross_Product_y(Particles.N_Particle[id]->edges[s],
                                          Particles.N_Particle[j]->edges[m]);
            L_axis[k].z =
                Calculate_Cross_Product_z(Particles.N_Particle[id]->edges[s],
                                          Particles.N_Particle[j]->edges[m]);

            norm_L = L_axis[k].norm();

            L_axis[k].x = L_axis[k].x / norm_L;
            L_axis[k].y = L_axis[k].y / norm_L;
            L_axis[k].z = L_axis[k].z / norm_L;

            k = k + 1;
        }
    }

    return N_separating_axis;
}

int pmove::Calculate_Separating_Axis_RHOMBI(particles &Particles, int id,
                                             int j, m_vector *L_axis) {

    int N_indep, N_cross;
    double norm_L;
    int N_separating_axis;

    N_indep = 2;
    N_cross = 0;

    N_separating_axis = 2 * N_indep;
   

    for (int k = 0; k < N_indep; k++) {

        L_axis[k].x = Particles.N_Particle[id]->facenormal[k+1].x;
        L_axis[k].y = Particles.N_Particle[id]->facenormal[k+1].y;
        L_axis[k].z = Particles.N_Particle[id]->facenormal[k+1].z;
    }

    for (int k = N_indep; k < (2 * N_indep); k++) {

        L_axis[k].x = Particles.N_Particle[j]->facenormal[k+1 - N_indep].x;
        L_axis[k].y = Particles.N_Particle[j]->facenormal[k+1 - N_indep].y;
        L_axis[k].z = Particles.N_Particle[j]->facenormal[k+1 - N_indep].z;
    }

   
    return N_separating_axis;
}



double pmove::Scalar_Product(m_vector a, m_vector b) {

    m_scalar_product = a.x * b.x + a.y * b.y + a.z * b.z;
    return m_scalar_product;
}

void pmove::Collision_Test(particles &Particles, box *Box, int id,
                           collision_list *Collision_List) {

    int cp_id = 0;
    int N_it;
    N_it = 0;
    Halo_Energy = 0;
    double f_halo;
    double f_exit_MIN = 3;
    double f_halo_MAX;

    Particles.N_Particle[id]->Calculate_Axis();
    Particles.N_Particle[id]->Calculate_Face_Normals();

    do {

        if (Collision_List[id].Nm > 0) {

            collision_partner = Collision_List[id].Elements[cp_id].nl_id;

            if (collision_partner != id) {

                Particles.N_Particle[collision_partner]
                    ->Calculate_Face_Normals();

                //N_it = Calculate_Separating_Axis_GENERAL(
                //           Particles, id, collision_partner, L_axis) -
                //       1;

                N_it = Calculate_Separating_Axis_RHOMBI(Particles, id, collision_partner, L_axis) - 1;


                int j;
                j = -1;

                f_exit_MAX = -3.0;
                f_halo_MAX = -3.0;
                f_exit = 0.0;
                f_halo = 0.0;

                do {

                    j = j + 1;

                    L_axis_norm2 = Scalar_Product(L_axis[j], L_axis[j]);

                    if (L_axis_norm2 > 0.5) {

                        dL = Scalar_Product(
                            Collision_List[id].Elements[cp_id].distance,
                            L_axis[j]);
                        dL = fabs(dL);

                        R_id = Particles.N_Particle[id]
                                   ->Calculate_Projection_to_Separating_Axis(
                                       L_axis[j]);
                        R_collision_partner =
                            Particles.N_Particle[collision_partner]
                                ->Calculate_Projection_to_Separating_Axis(
                                    L_axis[j]);

                        // if(Particles.N_Particle[id]->rmin <
                        // Particles.N_Particle[collision_partner]->rmin){
                        //	 dL=fabs(Particles.N_Particle[collision_partner]->rmax
                        //- Particles.N_Particle[id]->rmin);

                        //}

                        // if(Particles.N_Particle[collision_partner]->rmin <
                        // Particles.N_Particle[id]->rmin){
                        // dL=fabs(Particles.N_Particle[id]->rmax -
                        // Particles.N_Particle[collision_partner]->rmin);

                        //}

                        double a;
                        a = (Particles.N_Particle[id]->halo_cutoff / 2.0) +
                            Particles.N_Particle[id]->Lx;

                        f_exit = dL - R_id - R_collision_partner;
                        f_halo = dL - a * R_id - a * R_collision_partner;

                        if (f_exit > f_exit_MAX) {
                            f_exit_MAX = f_exit;
                        }

                        if (f_halo > f_halo_MAX) {
                            f_halo_MAX = f_halo;
                        }

                        // smallest positive value: minimum positive distance
                        if (f_exit > 0 && f_exit < f_exit_MIN) {
                            f_exit_MIN = f_exit;
                        }
                    }

                    // cout<<"f_exit_MAX "<<f_exit_MAX<<endl;

                    //}while((j!=N_it)&&(f_exit_MAX<0.0));
                } while (j != N_it);

                if (f_exit_MAX < 0.0) {

                    col_count = 1;
                }

                if (f_exit_MAX >= 0.0) {

                    col_count = 0;

                    if (f_halo_MAX < 0) {
                        Halo_Energy += Particles.N_Particle[id]->halo_energy;
                    }
                }

                exit_status = exit_status + col_count;
            }

            cp_id = cp_id + 1;
        }

    } while ((exit_status == 0) && (cp_id < Collision_List[id].Nm));
}

void pmove::Reset_Pair_Potential(particles &Particles, box *Box) {

    Total_Energy = Total_Energy_old;

    /*
    for (int id1 = 0 ;id1<Box->N;id1++){
            for (int id2 = 0 ;id2<Box->N; id2++){

                    Patch_Energies[id1][id2] = Patch_Energies_old[id1][id2];


            }


    }

    */
}

void pmove::Set_Pair_Potential(particles &Particles, box *Box) {

    Total_Energy_old = Total_Energy;

    /*
    for (int id1=0;id1<Box->N;id1++){
            for(int id2=0;id2<Box->N;id2++){

                    Patch_Energies_old[id1][id2] = Patch_Energies[id1][id2];

            }

    }
    */
}

double pmove::Calculate_Patch_Distance(int id1, int id2, int patch_id1,
                                       int patch_id2, particles &Particles,
                                       box *Box) {

    double particle_dist_x, particle_dist_y, particle_dist_z;
    double particle_dist_sq;

    particle_dist_x = Particles.N_Particle[id1]->x_patch[patch_id1] -
                      Particles.N_Particle[id2]->x_patch[patch_id2];
    particle_dist_y = Particles.N_Particle[id1]->y_patch[patch_id1] -
                      Particles.N_Particle[id2]->y_patch[patch_id2];
    particle_dist_z = Particles.N_Particle[id1]->z_patch[patch_id1] -
                      Particles.N_Particle[id2]->z_patch[patch_id2];

    particle_dist_x =
        particle_dist_x - Box->Lx * rint(particle_dist_x / Box->Lx);
    particle_dist_y =
        particle_dist_y - Box->Ly * rint(particle_dist_y / Box->Ly);
    particle_dist_z =
        particle_dist_z - Box->Lz * rint(particle_dist_z / Box->Lz);

    particle_dist_sq = particle_dist_x * particle_dist_x +
                       particle_dist_y * particle_dist_y +
                       particle_dist_z * particle_dist_z;

    return particle_dist_sq;
}

void pmove::Calculate_Pair_Potential(particles &Particles, box *Box) {

    Total_Energy = 0;
    int l, m;
    double patch_energy_ij;
    patch_energy_ij = 0;

    double patch_distance_squared;

    for (int id1 = 0; id1 < Box->N; id1++) {
        for (int id2 = 0; id2 < Box->N; id2++) {

            for (int pid1 = 0; pid1 < Particles.N_Particle[id1]->N_patches;
                 pid1++) {

                for (int pid2 = 0; pid2 < Particles.N_Particle[id2]->N_patches;
                     pid2++) {

                    patch_distance_squared = Calculate_Patch_Distance(
                        id1, id2, pid1, pid2, Particles, Box);

                    // ONLY if all the radii are the same width

                    if (patch_distance_squared <
                            Particles.N_Particle[id1]
                                ->patch_cutoff_squared[pid1] &&
                        id1 != id2) {

                      l = Particles.N_Particle[id1]->patch_type[pid1];
                      m = Particles.N_Particle[id2]->patch_type[pid2];

                        patch_energy_ij =
                            Particles.N_Particle[id1]->patch_energy[l][m];
                        Total_Energy += patch_energy_ij;
                    }
                }
            }
        }
    }

    Total_Energy = double(Total_Energy) / 2.0;
    Particles.Total_Energy = Total_Energy;
}

double pmove::Calculate_Pair_Potential(int id1, particles &Particles, box *Box,
                                       collision_list *Collision_List) {

    double delta_U;
    double patch_distance_squared;
    int id2;
    int l, m;
    double patch_energy_ij;
    patch_energy_ij = 0;

    delta_U = 0;

    for (int k = 0; k < Collision_List[id1].Nm; k++) {

        id2 = Collision_List[id1].Elements[k].nl_id;
        // Collision_List[id1].Elements[k].is_neighbour=0;

        for (int pid1 = 0; pid1 < Particles.N_Particle[id1]->N_patches;
             pid1++) {

            for (int pid2 = 0; pid2 < Particles.N_Particle[id2]->N_patches;
                 pid2++) {

                patch_distance_squared = Calculate_Patch_Distance(
                    id1, id2, pid1, pid2, Particles, Box);

                // ONLY if all the radii are the same width
                // if ( patch_distance_squared <
                // Particles.N_Particle[id1]->patch_cutoff_squared[pid1] &&
                // id1!= id2){

                if (patch_distance_squared <
                        Particles.N_Particle[id1]
                            ->patch_cutoff_squared[pid1] &&
                    id1 != id2) {

                      /////////////////////////////////////
                      ofstream fout("patch_state.dat", ios::out | ios::app);
                      double particle_dist_x;
                      double particle_dist_y; 
                      double dphi;

                      particle_dist_x = Particles.N_Particle[id1]->x_patch[pid1] -
                        Particles.N_Particle[id2]->x_patch[pid2];
                      particle_dist_y = Particles.N_Particle[id1]->y_patch[pid1] -
                        Particles.N_Particle[id2]->y_patch[pid2];

                      particle_dist_x =
                        particle_dist_x - Box->Lx * rint(particle_dist_x / Box->Lx);
                      particle_dist_y =
                        particle_dist_y - Box->Ly * rint(particle_dist_y / Box->Ly);

                      dphi = Particles.N_Particle[id1]->phi - Particles.N_Particle[id2]->phi; 

                      fout <<particle_dist_x<<" "<<particle_dist_y<<" "<<dphi<<endl;
                      fout.close();

                      //////////////////////////////////////////////


                    l = Particles.N_Particle[id1]->patch_type[pid1];
                    m = Particles.N_Particle[id2]->patch_type[pid2];

                    patch_energy_ij =
                        Particles.N_Particle[id1]->patch_energy[l][m];
                    delta_U += patch_energy_ij;
                }
            }
        }
    }

    return delta_U;
}
