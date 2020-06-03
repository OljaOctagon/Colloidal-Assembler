#include "pmove.h"

void pmove::Translate(particles &Particles, box *Box, fileio &Fileio, int id,
                      int mc_time) {

    double dU_old, dU_new;
    double delta_U;

    dU_old = 0;
    dU_new = 0;

    Particles.Collision_List[id].Calculate(Box, id, Particles.Id_Cell_List,
                                           Particles.Cell_List, Particles.Cell,
                                           Particles.N_Particle,Particles.N_Particle[0]->cut_off,
                                           Particles.MAX_coll_p);

    dU_old =
        Calculate_Pair_Potential(id, Particles, Box, Particles.Collision_List);


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


    Trans_Update_Positions(Particles, id, trans_vec);

    // is center of mass outside of the box?
    Particles.Check_Periodic_CM(id, Box);
    Update_Periodic_Positions(Particles, Box, id);
    Particles.N_Particle[id]->Calculate_Axis();
    Particles.N_Particle[id]->Calculate_Patch_Position();

    Particles.Update_Cell_List(id, Box); 

    Particles.Collision_List[id].Calculate(Box, id, Particles.Id_Cell_List,
                                            Particles.Cell_List, Particles.Cell,
                                            Particles.N_Particle,Particles.N_Particle[0]->cut_off,
                                            Particles.MAX_coll_p);
    exit_status = 0;
    col_count = 0;

    Collision_Test(Particles, Box, id, Particles.Collision_List);

    if (exit_status >= 1) {
        Reset_Positions(Particles, id);

        Particles.N_Particle[id]->Calculate_Axis();
        Particles.N_Particle[id]->Calculate_Patch_Position();

        Particles.Reset_Cell_List(Box, id);
    }


    if (exit_status == 0) {

        dU_new = Calculate_Pair_Potential(id, Particles, Box,
                                          Particles.Collision_List);

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

            Particles.Reset_Cell_List(Box, id);

        }

        // Accept move if XI <= b_factor

        if (XI <= b_factor) {

            Set_Positions(Particles, id);
            Particles.Total_Energy = Particles.Total_Energy - dU_old + dU_new;
            accept_translate = accept_translate + 1;
            N_trans_moves = N_trans_moves + 1;
            Particles.Set_Cell_List(Box,id);



        }
    }

}
