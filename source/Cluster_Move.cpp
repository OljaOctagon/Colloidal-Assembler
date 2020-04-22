#include "pmove.h"
#include <gsl/gsl_math.h>

double pmove::Calculate_Potential(int a, int b, particles &Particles,
                                  box *Box) {

    double patch_distance_squared;
    double epsilon_ij;
    int ela;
    int elb;
    double patch_energy_ab;

    epsilon_ij = 0;
    patch_energy_ab = 0;

    // TODO: transform to do-while looop as particles can be bonded only once.
    // save info abut is bonded state 
    for (int pid1 = 0; pid1 < Particles.N_Particle[a]->N_patches; pid1++) {

        for (int pid2 = 0; pid2 < Particles.N_Particle[b]->N_patches; pid2++) {

            patch_distance_squared =
                Calculate_Patch_Distance(a, b, pid1, pid2, Particles, Box);

            // ONLY if all the radii are the same width

            if (patch_distance_squared <
                    Particles.N_Particle[a]->patch_cutoff_squared[pid1] &&
                a != b) {

                ela = Particles.N_Particle[a]->patch_type[pid1];
                elb = Particles.N_Particle[b]->patch_type[pid2];

                patch_energy_ab =
                    Particles.N_Particle[a]->patch_energy[ela][elb];

                epsilon_ij = epsilon_ij + patch_energy_ab;
            }
        }
    }

    return epsilon_ij;
}

void pmove::Reset_Pseudo_Cluster(box *Box) {

    for (int i = 0; i < Box->N; i++) {
        is_element[i] = false;
        for (int j = 0; j < Box->N; j++) {
            pseudo_cluster_info[i][j] = -1;
        }
        List[i] = -2;
    }
}

void pmove::Pseudocluster_Recursion(int id_j, int cn, int fl,
                                    particles &Particles, box *Box) {

    int j;

    if (is_element[id_j] == false) {
        List[N_List] = id_j;
        N_List = N_List + 1;
        is_element[id_j] = true;
    }


    Particles.Collision_List[id_j].Calculate(Box, id_j, Particles.Id_Cell_List,
                                           Particles.Cell_List, Particles.Cell,
                                           Particles.N_Particle,Particles.N_Particle[0]->cut_off,
                                           Particles.MAX_coll_p);

    //Particles.Collision_List[id_j].Calculate_OP(
    //    Box, id_j, Particles.N_Particle, Particles.N_Particle[0]->cut_off,
    //    Particles.MAX_coll_p);

    for (int cp = 0; cp < Particles.Collision_List[id_j].Nm; cp++) {

        j = Particles.Collision_List[id_j].Elements[cp].nl_id;

        if (pseudo_cluster_info[id_j][j] == -1) {

            e_pair = Calculate_Potential(id_j, j, Particles, Box);

            pseudo_cluster_info[id_j][j] = 1;
            pseudo_cluster_info[j][id_j] = 1;

            b_factor = 1 - exp(beta_f * (e_pair));
            b_factor = GSL_MAX(0, b_factor);

            XI = gsl_rng_uniform(r01);

            if (XI < b_factor) {

                cn = N_Bonds;
                N_Bonds = N_Bonds + 1;

                Pseudocluster_Recursion(j, cn, fl, Particles, Box);
            }
        }
    }
}

void pmove::Rot_Cluster_Move(particles &Particles, box *Box, fileio &Fileio,
                             int mc_time) {

    //double Total_Energy_old;
    double delta_U;
   
    N_List = 0;
    N_Bonds = 0;
    r_id = -1;
    double phi_t;
    double phi_r;

    // choose particle id

    id = gsl_rng_uniform_int(r, Box->N);

    Reset_Pseudo_Cluster(Box);

    // calculate translation

    // Rot_Update_Quarternions_VON_MISES(Particles, id);
    double rand_phi;

    // rand_phi = gsl_ran_gaussian(r01, kappa);
    rand_phi = gsl_rng_uniform(r01);

    phi_t = kappa * (rand_phi - 0.5);

    Rot_mat[0] = cos(phi_t);
    Rot_mat[1] = -sin(phi_t);
    Rot_mat[2] = 0;

    Rot_mat[3] = sin(phi_t);
    Rot_mat[4] = cos(phi_t);
    Rot_mat[5] = 0;

    Rot_mat[6] = 0;
    Rot_mat[7] = 0;
    Rot_mat[8] = 1;

    int links, failed_links;
    links = 0;
    failed_links = 0;

    Pseudocluster_Recursion(id, links, failed_links, Particles, Box);

    // Collision Test

    double theta_x, theta_y, theta_av_x, theta_av_y;
    double xp1, xp2, yp1, yp2;

    xp1 = 0;
    xp2 = 0;
    yp1 = 0;
    yp2 = 0;

    for (int k = 0; k < N_List; k++) {
        int j;
        j = List[k];
        theta_x = (Particles.N_Particle[j]->x_center / Box->Lx) * 2 * M_PI;
        theta_y = (Particles.N_Particle[j]->y_center / Box->Ly) * 2 * M_PI;

        xp1 = xp1 + (Box->Lx / (2 * M_PI)) * cos(theta_x);
        xp2 = xp2 + (Box->Lx / (2 * M_PI)) * sin(theta_x);

        yp1 = yp1 + (Box->Ly / (2 * M_PI)) * cos(theta_y);
        yp2 = yp2 + (Box->Ly / (2 * M_PI)) * sin(theta_y);
    }

    xp1 = xp1 / double(N_List);
    xp2 = xp2 / double(N_List);
    yp1 = yp1 / double(N_List);
    yp2 = yp2 / double(N_List);
    N_List = int(N_List);

    theta_av_x = atan2(-xp2, -xp1) + M_PI;
    theta_av_y = atan2(-yp2, -yp1) + M_PI;

    center_mass.x = (Box->Lx / (2 * M_PI)) * (theta_av_x);
    center_mass.y = (Box->Ly / (2 * M_PI)) * (theta_av_y);
    center_mass.z = Box->Lz / 2.0;


    // Calculate old energy difference 
    double dU_old;
    dU_old = 0;

    for (int k = 0; k < N_List; k++){
      int j;
      j = List[k];
      dU_old += Calculate_Pair_Potential(j, Particles, Box,
                                        Particles.Collision_List);

    }


    for (int k = 0; k < N_List; k++) {
        int j;
        j = List[k];
        Rot_Move_Map(Particles, j, Box, center_mass, Rot_mat);
        Particles.Check_Periodic_CM(j, Box);
        Update_Periodic_Positions(Particles, Box, j);
        Particles.N_Particle[j]->Calculate_Axis();
        Particles.N_Particle[j]->Calculate_Patch_Position();
        Particles.Update_Cell_List(j, Box);
    }

    int k = 0;

    exit_status = 0;
    col_count = 0;

    int Max_Length;
    Max_Length = rint(1. / gsl_rng_uniform(r01));

    if (N_List > Max_Length) {
        exit_status = 1;
    }

    do {
        Particles.Collision_List[List[k]].Calculate(Box, List[k], Particles.Id_Cell_List,
                                               Particles.Cell_List, Particles.Cell,
                                               Particles.N_Particle,Particles.N_Particle[0]->cut_off,
                                               Particles.MAX_coll_p);


        // Particles.Collision_List[List[k]].Calculate_OP(
        //    Box, List[k], Particles.N_Particle,
        //    Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);

        // Collision Test
        Collision_Test(Particles, Box, List[k], Particles.Collision_List);
        k++;

    } while ((exit_status == 0) && (k < N_List));


    if (exit_status > 0) {

        // Reset all positions!
        for (int k = 0; k < N_List; k++) {
            int j;
            j = List[k];
            Reset_Positions(Particles, j);
            Particles.N_Particle[j]->Calculate_Axis();
            Particles.N_Particle[j]->Calculate_Patch_Position();
            // Set_Positions(Particles, j);
        }
        Particles.Reset_Cell_List(Box);
        Particles.Set_Cell_List(Box);

    }

    if (exit_status == 0) {


        //Total_Energy_old = Particles.Total_Energy;
        double dU_new;
        dU_new = 0;

        for (int k = 0; k < N_List; k++){
            int j;
            j = List[k];
            dU_new += Calculate_Pair_Potential(j, Particles, Box,
                                   Particles.Collision_List);

        }

        // TODO make less expensive! one by one calculation
        //Calculate_Pair_Potential(Particles, Box);
        //delta_U = Total_Energy - Total_Energy_old;

        //double delta_U_test;
        delta_U  = dU_new - dU_old; 

        //cout<<"rot    "<<delta_U<<"   "<<delta_U_test<<endl;  

        b_factor = GSL_MIN(1, exp((beta_f - beta) * delta_U));
        XI = gsl_rng_uniform(r01);

        if (b_factor < XI) {

            for (int k = 0; k < N_List; k++) {
                int j;
                j = List[k];
                Reset_Positions(Particles, j);
                Particles.N_Particle[j]->Calculate_Axis();
                Particles.N_Particle[j]->Calculate_Patch_Position();

            }
            Particles.Reset_Cell_List(Box);
            Particles.Set_Cell_List(Box);
        }


        if (b_factor >= XI) {


            for (int k = 0; k < N_List; k++) {
                int j;
                j = List[k];

                Particles.N_Particle[j]->phi =
                    Particles.N_Particle[j]->phi + phi_t;

                if (Particles.N_Particle[j]->phi > 2 * M_PI) {
                    Particles.N_Particle[j]->phi =
                        Particles.N_Particle[j]->phi - 2.0 * M_PI;
                }

                if (Particles.N_Particle[j]->phi < 0) {
                    Particles.N_Particle[j]->phi =
                        Particles.N_Particle[j]->phi + 2.0 * M_PI;
                }

                Set_Positions(Particles, j);
            }

            Particles.Total_Energy = Particles.Total_Energy - dU_old + dU_new;

        }
    }

    Reset_Pseudo_Cluster(Box);
}



void pmove::Trans_Cluster_Move(particles &Particles, box *Box, fileio &Fileio,
                               int mc_time) {
    N_List = 0;
    N_Bonds = 0;
    N_failed_links = 0;
    //double Total_Energy_old;
    double delta_U;

    // choose particle id

    id = gsl_rng_uniform_int(r, Box->N);

    Reset_Pseudo_Cluster(Box);

    // calculate translation

    rand_x = gsl_ran_gaussian(r01, sigma_trans);
    rand_y = gsl_ran_gaussian(r01, sigma_trans);
    rand_z = gsl_ran_gaussian(r01, sigma_trans);

    trans_vec.x = rand_x;
    trans_vec.y = rand_y;

    if (is_2D == 1) {
        trans_vec.z = 0.0;
    }

    else {
        trans_vec.z = rand_z;
    }

    int links, failed_links;
    links = 0;
    failed_links = 0;

    Pseudocluster_Recursion(id, links, failed_links, Particles, Box);

    // Calculate old energy difference 
    double dU_old;
    dU_old = 0;

    for (int k = 0; k < N_List; k++){
      int j;
      j = List[k];
      dU_old += Calculate_Pair_Potential(j, Particles, Box,
                                         Particles.Collision_List);

    }
    // Update positions  

    for (int k = 0; k < N_List; k++) {

        Trans_Update_Positions(Particles, List[k], trans_vec);
        Particles.Check_Periodic_CM(List[k], Box);
        Update_Periodic_Positions(Particles, Box, List[k]);
        Particles.N_Particle[List[k]]->Calculate_Axis();
        Particles.N_Particle[List[k]]->Calculate_Patch_Position();
        Particles.Update_Cell_List(List[k], Box);

    }

    int k = 0;

    exit_status = 0;
    col_count = 0;

    int Max_Length;
    Max_Length = rint(1. / gsl_rng_uniform(r01));

    if (N_List > Max_Length) {

        exit_status = 1;
    }

    do {
        Particles.Collision_List[List[k]].Calculate(Box, List[k],
                                                    Particles.Id_Cell_List,
                                                    Particles.Cell_List, Particles.Cell,
                                                    Particles.N_Particle,
                                                    Particles.N_Particle[0]->cut_off,
                                                    Particles.MAX_coll_p);

        //Particles.Collision_List[List[k]].Calculate_OP(
        //    Box, List[k], Particles.N_Particle,
        //    Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);

        // Collision Test
        Collision_Test(Particles, Box, List[k], Particles.Collision_List);
        k++;

    } while ((exit_status == 0) && (k < N_List));

    if (exit_status > 0) {

        for (int k = 0; k < N_List; k++) {
            Reset_Positions(Particles, List[k]);
            Particles.N_Particle[List[k]]->Calculate_Axis();
            Particles.N_Particle[List[k]]->Calculate_Patch_Position();

        }

        Particles.Reset_Cell_List(Box);
        Particles.Set_Cell_List(Box);

    }

    if (exit_status == 0) {

        double dU_new;
        dU_new = 0;

        for (int k = 0; k < N_List; k++){
          int j;
          j = List[k];
          dU_new += Calculate_Pair_Potential(j, Particles, Box,
                                            Particles.Collision_List);

        }

        //double delta_U_test;
        delta_U = dU_new - dU_old; 

        // Total_Energy_old = Particles.Total_Energy;
        // TODO: make faster!
        //Calculate_Pair_Potential(Particles, Box);
        //delta_U = Total_Energy - Total_Energy_old;

        //cout<<"trans    "<<delta_U<<"   "<<delta_U_test<<endl;  

        b_factor = GSL_MIN(1, exp((beta_f - beta) * delta_U));
        XI = gsl_rng_uniform(r01);

        if (b_factor < XI) {

            for (int k = 0; k < N_List; k++) {
                int j;
                j = List[k];
                Reset_Positions(Particles, j);
                Particles.N_Particle[j]->Calculate_Axis();
                Particles.N_Particle[j]->Calculate_Patch_Position();

            }
            Particles.Reset_Cell_List(Box);
            Particles.Set_Cell_List(Box);

            Particles.Total_Energy = Particles.Total_Energy - dU_old + dU_new;
            //Particles.Total_Energy = Total_Energy_old;
            //Total_Energy = Total_Energy_old;
        }

        if (b_factor >= XI) {
            for (int k = 0; k < N_List; k++) {
                int j;
                j = List[k];
                Set_Positions(Particles, j);
            }


            //Particles.Total_Energy = Total_Energy;
        }
    }

    Reset_Pseudo_Cluster(Box);
}

//// Dynamic linking scheme: complicated and broken


/*

void pmove::Trans_Pseudocluster_Recursion(int id_j, int cn, int fl,
                                          particles &Particles, box *Box) {

    int j;

    if (is_element[id_j] == false) {
        List[N_List] = id_j;
        N_List = N_List + 1;
        is_element[id_j] = true;
    }

    Particles.Collision_List[id_j].Calculate_OP(
        Box, id_j, Particles.N_Particle, Particles.N_Particle[0]->cut_off,
        Particles.MAX_coll_p);

    for (int cp = 0; cp < Particles.Collision_List[id_j].Nm; cp++) {
        Particles.Collision_List_old[id_j].Elements[cp] =
            Particles.Collision_List[id_j].Elements[cp];
    }

    Particles.Collision_List_old[id_j].Nm = Particles.Collision_List[id_j].Nm;

    Trans_Update_Positions(Particles, id_j, trans_vec);
    Particles.Check_Periodic_CM(id_j, Box);
    Update_Periodic_Positions(Particles, Box, id_j);
    Particles.N_Particle[id_j]->Calculate_Axis();
    Particles.N_Particle[id_j]->Calculate_Patch_Position();

    Particles.Collision_List[id_j].Calculate_OP(
        Box, id_j, Particles.N_Particle, Particles.N_Particle[0]->cut_off,
        Particles.MAX_coll_p);

    Reset_Positions(Particles, id_j);
    Particles.N_Particle[id_j]->Calculate_Axis();
    Particles.N_Particle[id_j]->Calculate_Patch_Position();

    for (int cp = 0; cp < Particles.Collision_List[id_j].Nm; cp++) {

        Trans_Update_Positions(Particles, id_j, trans_vec);
        Particles.Check_Periodic_CM(id_j, Box);
        Update_Periodic_Positions(Particles, Box, id_j);
        Particles.N_Particle[id_j]->Calculate_Axis();
        Particles.N_Particle[id_j]->Calculate_Patch_Position();

        j = Particles.Collision_List[id_j].Elements[cp].nl_id;

        if (pseudo_cluster_info[id_j][j] == -1) {

            e_single = Calculate_Potential(id_j, j, Particles, Box);

            pseudo_cluster_info[id_j][j] = 1;
            pseudo_cluster_info[j][id_j] = 1;
            // if (is_interacting=true){

            Particles.Collision_List[j].Calculate_OP(
                Box, j, Particles.N_Particle, Particles.N_Particle[0]->cut_off,
                Particles.MAX_coll_p);

            for (int cp = 0; cp < Particles.Collision_List[j].Nm; cp++) {
                Particles.Collision_List_old[j].Elements[cp] =
                    Particles.Collision_List[j].Elements[cp];
            }

            Particles.Collision_List_old[j].Nm =
                Particles.Collision_List[j].Nm;

            Trans_Update_Positions(Particles, j, trans_vec);
            Particles.Check_Periodic_CM(j, Box);
            Update_Periodic_Positions(Particles, Box, j);
            Particles.N_Particle[j]->Calculate_Axis();
            Particles.N_Particle[j]->Calculate_Patch_Position();

            e_pair = Calculate_Potential(id_j, j, Particles, Box);

            q_factor = exp(beta * (e_pair - e_single));
            b_factor = 1 - exp(beta * (e_pair - e_single));

            b_factor = GSL_MAX(0, b_factor);
            q_factor = GSL_MIN(1, q_factor);

            XI = gsl_rng_uniform(r01);

            Reset_Positions(Particles, id_j);
            Particles.N_Particle[id_j]->Calculate_Axis();
            Particles.N_Particle[id_j]->Calculate_Patch_Position();

            if (XI >= b_factor) {

                Reset_Positions(Particles, j);
                Particles.N_Particle[j]->Calculate_Axis();
                Particles.N_Particle[j]->Calculate_Patch_Position();

                if (is_element[j] == true) {

                    fl = N_failed_links;

                    q_f[fl] = q_factor;

                    Trans_Update_Positions(Particles, id_j, r_trans_vec);
                    Particles.Check_Periodic_CM(id_j, Box);
                    Update_Periodic_Positions(Particles, Box, id_j);
                    Particles.N_Particle[id_j]->Calculate_Axis();
                    Particles.N_Particle[id_j]->Calculate_Patch_Position();

                    e_single = Calculate_Potential(id_j, j, Particles, Box);

                    q_factor = exp(beta * (e_pair - e_single));
                    q_factor = GSL_MIN(1, q_factor);
                    q_r[fl] = q_factor;

                    Reset_Positions(Particles, id_j);
                    Particles.N_Particle[id_j]->Calculate_Axis();
                    Particles.N_Particle[id_j]->Calculate_Patch_Position();

                    N_failed_links = N_failed_links + 1;
                }
            }

            if (XI < b_factor) {

                // distinguish number of bonds and number of particles
                // that go into list

                cn = N_Bonds;
                p_f[cn] = b_factor;
                // VIRTUAL opposite move
                // go with linker particle in oppositie direction

                Trans_Update_Positions(Particles, id_j, r_trans_vec);
                Particles.Check_Periodic_CM(id_j, Box);
                Update_Periodic_Positions(Particles, Box, id_j);
                Particles.N_Particle[id_j]->Calculate_Axis();
                Particles.N_Particle[id_j]->Calculate_Patch_Position();

                // reset positions of linkee

                Reset_Positions(Particles, j);
                Particles.N_Particle[j]->Calculate_Axis();
                Particles.N_Particle[j]->Calculate_Patch_Position();

                // calculate single potential again and reverse linking
                // probablility
                e_single = Calculate_Potential(id_j, j, Particles, Box);


                b_factor = 1 - exp(beta * (e_pair - e_single));
                b_factor = GSL_MAX(0, b_factor);
                p_r[cn] = b_factor;

                // move linker back to new move

                Reset_Positions(Particles, id_j);
                Particles.N_Particle[id_j]->Calculate_Axis();
                Particles.N_Particle[id_j]->Calculate_Patch_Position();

                N_Bonds = N_Bonds + 1;

                // move linkee back to new move

                Trans_Pseudocluster_Recursion(j, cn, fl, Particles, Box);
            }
        }
    }
}
*/
