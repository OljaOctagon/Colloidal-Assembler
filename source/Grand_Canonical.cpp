#include "pmove.h"

void pmove::Particle_Insertion(particles &Particles, box *Box, fileio &Fileio,
                               int mc_time) {

    id = Box->N;
    double phi_t;
    double dU_old, dU_new;
    double delta_U;
    dU_old = 0;
    dU_new = 0;

    /*
    Particles.Collision_List[id].Calculate_OP(Box, id, Particles.N_Particle,
    Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p); dU_old =
    Calculate_Pair_Potential(id, Particles, Box, Particles.Collision_List);

    exit_status = 0;
    col_count = 0;

    Collision_Test( Particles, Box, id, Particles.Collision_List);
    dU_old = dU_old + Halo_Energy;
    */

    double x_trial, y_trial, z_trial;

    rand_x = gsl_rng_uniform(r01);
    rand_y = gsl_rng_uniform(r01);

    // pick position at random in Box

    x_trial = Box->x[0] + Box->Lx * rand_x;
    y_trial = Box->y[0] + Box->Ly * rand_y;
    z_trial = Particles.N_Particle[0]->z_center;

    // pick random rotation for the particle

    Particles.N_Particle[id]->phi_old = Particles.N_Particle[id]->phi;
    Particles.N_Particle[id]->phi = 2.0 * M_PI * (gsl_rng_uniform(r01));

    // set particle patch topology randomly but according to phi_binary
    double rand_s;
    int patch_1, patch_2, patch_3, patch_4, patch_5, patch_6;
    int p_type;

    double mu_current;
    mu_current = Box->mu_1;


    if (Particles.binary_on.compare("on") == 0) {

        patch_1 = 0;
        patch_2 = 0;
        patch_3 = 0;
        patch_4 = 0;
        patch_5 = 0;
        patch_6 = 0;
        rand_s = gsl_rng_uniform(r01);


        // particles types are sampled uniformly. The different chemical potentials
        // regulate the proportions.
        if (rand_s < 0.5) {
            patch_1 = 2;
            patch_2 = 1;
            patch_3 = 0;
            patch_4 = 0;
            patch_5 = 0;
            patch_6 = 0;

            mu_current = Box->mu_1;

        }

        if (rand_s >= 0.5) {
            patch_1 = 1;
            patch_2 = 3;
            patch_3 = 0;
            patch_4 = 0;
            patch_5 = 0;
            patch_6 = 0;

            mu_current = Box->mu_2;

        }


        Particles.N_Particle[id]->Set_Lengths(patch_1, patch_2, patch_3,
                                              patch_4, patch_5, patch_6);
        Particles.N_Particle_old[id]->Set_Lengths(patch_1, patch_2, patch_3,
                                                  patch_4, patch_5, patch_6);
    }

    if (Particles.ternary_on.compare("on") == 0) {

      patch_1 = 0;
      patch_2 = 0;
      patch_3 = 0;
      patch_4 = 0;
      patch_5 = 0;
      patch_6 = 0;
      rand_s = gsl_rng_uniform(r01);


      // particles types are sampled uniformly. The different chemical potentials
      // regulate the proportions.
      if (rand_s < 1./3.) {
        patch_1 = 2;
        patch_2 = 1;
        patch_3 = 0;
        patch_4 = 0;
        patch_5 = 0;
        patch_6 = 0;

        mu_current = Box->mu_1;

      }

      if ((rand_s >= 1./3.) && (rand_s<2./3.)) {
        patch_1 = 1;
        patch_2 = 3;
        patch_3 = 0;
        patch_4 = 0;
        patch_5 = 0;
        patch_6 = 0;

        mu_current = Box->mu_2;

      }

      if (rand_s >= 2./3.){ 
        patch_1 = 0;
        patch_2 = 1;
        patch_3 = 1;
        patch_4 = 0;
        patch_5 = 0;
        patch_6 = 0;

        mu_current = Box->mu_3;

      }

      Particles.N_Particle[id]->Set_Lengths(patch_1, patch_2, patch_3,
                                            patch_4, patch_5, patch_6);
      Particles.N_Particle_old[id]->Set_Lengths(patch_1, patch_2, patch_3,
                                                patch_4, patch_5, patch_6);
    }

    phi_t = Particles.N_Particle[id]->phi;

    Rot_mat[0] = cos(phi_t);
    Rot_mat[1] = -sin(phi_t);
    Rot_mat[2] = 0;

    Rot_mat[3] = sin(phi_t);
    Rot_mat[4] = cos(phi_t);
    Rot_mat[5] = 0;

    Rot_mat[6] = 0;
    Rot_mat[7] = 0;
    Rot_mat[8] = 1;

    Particles.N_Particle[id]->x_center = x_trial;
    Particles.N_Particle[id]->y_center = y_trial;
    Particles.N_Particle[id]->z_center = z_trial;

    Rot_Update_Positions(Particles, id, Rot_mat);
    Particles.N_Particle[id]->Calculate_Axis();
    Particles.N_Particle[id]->Calculate_Patch_Position();

    Particles.Collision_List[id].Calculate_OP(Box, id, Particles.N_Particle,
                                              Particles.N_Particle[0]->cut_off,
                                              Particles.MAX_coll_p);

    exit_status = 0;
    col_count = 0;

    // Collision Test
    Collision_Test(Particles, Box, id, Particles.Collision_List);

    if (exit_status > 0) {

        Reset_Positions(Particles, id);
        Particles.N_Particle[id]->Calculate_Axis();
        Particles.N_Particle[id]->Calculate_Patch_Position();
    }

    if (exit_status == 0) {

        // Calculate Pair Potential

        Box->N = Box->N + 1;

        dU_new = Calculate_Pair_Potential(id, Particles, Box,
                                          Particles.Collision_List);
        dU_new = dU_new + Halo_Energy;

        // Set_Pair_Potential(Particles, Box);

        delta_U = dU_new;

        // calculate Boltzmann-factor
        // cout<<"chemical potential"<<Box->mu<<endl;
        // cout<<"delta_U"<<delta_U<<endl;
        // cout<<"Total_Energy"<<Total_Energy<<endl;

        double Vf;
        Vf = Box->V / double(Box->N + 1);

        b_factor_pre = Vf * exp(-1.0 * beta * delta_U + beta * mu_current);
        b_factor = minimum(1, b_factor_pre);

        // cout<<"b_factor_pre "<<b_factor_pre<<endl;

        XI = gsl_rng_uniform(r01);

        // cout<<"bfactor  "<<b_factor<<endl;

        if (b_factor < XI) {

            Box->N = Box->N - 1;

            Reset_Positions(Particles, id);
            Particles.N_Particle[id]->Calculate_Axis();
            Particles.N_Particle[id]->Calculate_Patch_Position();
        }

        if (b_factor >= XI) {

            Set_Positions(Particles, id);
            Particles.Total_Energy = Particles.Total_Energy + delta_U;
            Box->packing_fraction =
                (Particles.N_Particle[0]->A * double(Box->N)) / Box->A;
        }
    }
}

void pmove::Particle_Deletion(particles &Particles, box *Box, fileio &Fileio,
                              int mc_time) {

    id = gsl_rng_uniform_int(r, Box->N);

    // Calculate old pair Potential

    double phi_t;
    double dU_old, dU_new;
    double delta_U;
    double mu_current;
    dU_old = 0;
    dU_new = 0;

    mu_current = Box->mu_1;

    Particles.Collision_List[id].Calculate_OP(Box, id, Particles.N_Particle,
                                              Particles.N_Particle[0]->cut_off,
                                              Particles.MAX_coll_p);
    dU_old =
        Calculate_Pair_Potential(id, Particles, Box, Particles.Collision_List);

    exit_status = 0;
    col_count = 0;

    Collision_Test(Particles, Box, id, Particles.Collision_List);
    dU_old = dU_old + Halo_Energy;

    dU_new = 0;

    delta_U = dU_new - dU_old;

    // calculate Boltzmann-factor

    // get the different chemical potential for binary mixtures

    double p0;
    if (Particles.binary_on.compare("on") == 0) {
        p0 = Particles.N_Particle[id]->patch_type[0];
        //NOTE: dependent on particle type!!!  
        if (p0 == 2){
            mu_current = Box->mu_1;
        }
        if (p0 == 1){
            mu_current = Box->mu_2;
        }

        if (p0 == 0){
           mu_current = Box->mu_3;
        }
    }

    b_factor_pre = (double(Box->N + 1) / double(Box->V)) *
                   exp(-beta * delta_U - beta * mu_current);
    b_factor = minimum(1, b_factor_pre);

    XI = gsl_rng_uniform(r01);

    if (b_factor < XI) {
    }

    if (b_factor >= XI) {

        Box->packing_fraction =
            (Particles.N_Particle[0]->A * double(Box->N - 1)) / Box->A;

        for (int k = id + 1; k < Box->N; k++) {
            Particles.N_Particle[k - 1]->x_center =
                Particles.N_Particle[k]->x_center;
            Particles.N_Particle[k - 1]->y_center =
                Particles.N_Particle[k]->y_center;
            Particles.N_Particle[k - 1]->z_center =
                Particles.N_Particle[k]->z_center;

            Particles.N_Particle[k - 1]->phi = Particles.N_Particle[k]->phi;

            Particles.N_Particle[k - 1]->q.x =
                Particles.N_Particle_old[k]->q.x;
            Particles.N_Particle[k - 1]->q.y =
                Particles.N_Particle_old[k]->q.y;
            Particles.N_Particle[k - 1]->q.z =
                Particles.N_Particle_old[k]->q.z;
            Particles.N_Particle[k - 1]->q.w =
                Particles.N_Particle_old[k]->q.w;

            Particles.N_Particle[k - 1]->patch_type[0] =
                Particles.N_Particle[k]->patch_type[0];
            Particles.N_Particle[k - 1]->patch_type[1] =
                Particles.N_Particle[k]->patch_type[1];
            Particles.N_Particle[k - 1]->patch_type[2] =
                Particles.N_Particle[k]->patch_type[2];
            Particles.N_Particle[k - 1]->patch_type[3] =
                Particles.N_Particle[k]->patch_type[3];

            for (int l = 0; l < Particles.N_Particle[k]->edge_N; l++) {
                Particles.N_Particle[k - 1]->x[l] =
                    Particles.N_Particle[k]->x[l];
                Particles.N_Particle[k - 1]->y[l] =
                    Particles.N_Particle[k]->y[l];
                Particles.N_Particle[k - 1]->z[l] =
                    Particles.N_Particle[k]->z[l];
            }

            Set_Positions(Particles, k - 1);
            Particles.N_Particle[k - 1]->Calculate_Axis();
            Particles.N_Particle[k - 1]->Calculate_Patch_Position();
        }

        Particles.N_Particle[Box->N - 1]->x_center = -10.0;
        Particles.N_Particle[Box->N - 1]->y_center = -10.0;
        Particles.N_Particle[Box->N - 1]->z_center = -10.0;

        Particles.N_Particle[Box->N - 1]->edges_from_center();
        Particles.N_Particle[Box->N - 1]->Calculate_Axis();
        Particles.N_Particle[Box->N - 1]->Calculate_Patch_Position();

        Particles.N_Particle[Box->N - 1]->q.x = 0.0;
        Particles.N_Particle[Box->N - 1]->q.y = 0.0;
        Particles.N_Particle[Box->N - 1]->q.z = 0.0;
        Particles.N_Particle[Box->N - 1]->q.w = 1.0;
        Particles.N_Particle[Box->N - 1]->phi = 0;

        Set_Positions(Particles, Box->N - 1);

        Box->N = Box->N - 1;

        Particles.Total_Energy = Particles.Total_Energy + delta_U;
    }
}
