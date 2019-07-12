

// Rot_Update_Random_2D
// Rot_mat
// N_Particle.phi

// reacalculate to quaternions

// Initialize Rotate2D function

// in para_in say wheter =2D

#include "pmove.h"

void pmove::Rot_Update_Random_2D(particles &Particles, int id) {

    double rand_phi;

    // rand_phi = gsl_ran_gaussian(r01, kappa);
    rand_phi = gsl_rng_uniform(r01);

    phi_reverse = Particles.N_Particle[id]->phi - kappa * (rand_phi - 0.5);

    Particles.N_Particle[id]->phi =
        Particles.N_Particle[id]->phi + kappa * (rand_phi - 0.5);

    if (Particles.N_Particle[id]->phi > 2 * M_PI) {
        Particles.N_Particle[id]->phi =
            Particles.N_Particle[id]->phi - 2.0 * M_PI;
    }

    if (Particles.N_Particle[id]->phi < 0) {
        Particles.N_Particle[id]->phi =
            Particles.N_Particle[id]->phi + 2.0 * M_PI;
    }

    if (phi_reverse > 2 * M_PI) {
        phi_reverse = phi_reverse - 2.0 * M_PI;
    }

    if (phi_reverse < 0) {
        phi_reverse = phi_reverse + 2.0 * M_PI;
    }
}

void pmove::Rotate2D(particles &Particles, box *Box, fileio &Fileio, int id,
                     int mc_time) {

    double phi_t;
    double dU_old, dU_new;
    double delta_U;
    dU_old = 0;
    dU_new = 0;

    Particles.Collision_List[id].Calculate_OP(Box, id, Particles.N_Particle,
                                              Particles.N_Particle[0]->cut_off,
                                              Particles.MAX_coll_p);

    exit_status = 0;
    col_count = 0;

    Collision_Test(Particles, Box, id, Particles.Collision_List);
    dU_old =
        Calculate_Pair_Potential(id, Particles, Box, Particles.Collision_List);
    dU_old = dU_old + Halo_Energy;
    // cout<<"Halo_Energy "<<Halo_Energy<<endl;
    // cout<<"dU_old "<<dU_old<<endl;

    Particles.N_Particle[id]->phi_old = Particles.N_Particle[id]->phi;

    // Rot_Update_Quarternions_VON_MISES(Particles, id);
    Rot_Update_Random_2D(Particles, id);

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

    // update positions

    // for gears uncomment edges from center und comment in rot_update,
    // checkperiodic and update periodic positions

    // Particles.N_Particle[id]->edges_from_center();

    Rot_Update_Positions(Particles, id, Rot_mat);
    // Which edges are outside the box

    // Particles.Check_Periodic_CM(id, Box);

    // is particle center outside the box?
    // Update_Periodic_Positions(Particles, Box, id);

    // update axis
    Particles.N_Particle[id]->Calculate_Axis();
    Particles.N_Particle[id]->Calculate_Patch_Position();

    // update cell list
    // Particles.Update_Cell_List(id, Box);

    // Calculate Collision Partners
    // Particles.Collision_List[id].Calculate(Box, id, Particles.Id_Cell_List,
    // Particles.Cell_List, Particles.Cell, Particles.N_Particle,
    // Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);
    // whitout cell-lists!!!
    Particles.Collision_List[id].Calculate_OP(Box, id, Particles.N_Particle,
                                              Particles.N_Particle[0]->cut_off,
                                              Particles.MAX_coll_p);

    exit_status = 0;
    col_count = 0;

    // Collision Test
    Collision_Test(Particles, Box, id, Particles.Collision_List);

    if (exit_status >= 1) {

        Reset_Positions(Particles, id);
        Particles.N_Particle[id]->Calculate_Axis();
        Particles.N_Particle[id]->Calculate_Patch_Position();
    }

    if (exit_status == 0) {

        dU_new = Calculate_Pair_Potential(id, Particles, Box,
                                          Particles.Collision_List);
        dU_new = dU_new + Halo_Energy;

        delta_U = dU_new - dU_old;

        // MC Step for Potential

        b_factor_pre = exp(-1.0 * beta * delta_U);
        b_factor = minimum(1, b_factor_pre);

        XI = gsl_rng_uniform(r01);

        // Reject of XI > b_factor

        if (XI > b_factor) {

            Reset_Positions(Particles, id);
            Particles.N_Particle[id]->Calculate_Axis();
            Particles.N_Particle[id]->Calculate_Patch_Position();
        }

        // Accept move if XI <= b_factor

        if (XI <= b_factor) {

            Particles.N_Particle[id]->phi_old = Particles.N_Particle[id]->phi;
            Set_Positions(Particles, id);

            Particles.Total_Energy = Particles.Total_Energy + delta_U;

            N_rotate_moves = N_rotate_moves + 1;
        }
    }
}
