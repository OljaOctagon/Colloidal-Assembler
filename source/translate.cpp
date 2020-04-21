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
    //Particles.Collision_List[id].Calculate_OP(Box, id, Particles.N_Particle,
    //                                          Particles.N_Particle[0]->cut_off,
    //                                          Particles.MAX_coll_p);
    dU_old =
        Calculate_Pair_Potential(id, Particles, Box, Particles.Collision_List);

    exit_status = 0;
    col_count = 0;

    Collision_Test(Particles, Box, id, Particles.Collision_List);
    dU_old = dU_old + Halo_Energy;
    // cout<<" Old Halo_Energy "<<Halo_Energy<<endl;
    // cout<<"dU_old"<<dU_old<<endl;

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

    ////cout<<"x_center "<<Particles.N_Particle[id]->x_center<<"  "<<id<<endl;

    // Update position

    Trans_Update_Positions(Particles, id, trans_vec);

    // is center of mass outside of the box?

    ////cout<<"x_center "<<Particles.N_Particle[id]->x_center<<"  "<<id<<endl;

    Particles.Check_Periodic_CM(id, Box);

    ////cout<<"x_center "<<Particles.N_Particle[id]->x_center<<"  "<<id<<endl;

    Update_Periodic_Positions(Particles, Box, id);

    /*
    Particles.N_Particle[id]->x_center = Particles.N_Particle[id]->x_center -
    Box->Lx/2.; Particles.N_Particle[id]->y_center =
    Particles.N_Particle[id]->y_center - Box->Lx/2.;

    Particles.N_Particle[id]->x_center = Particles.N_Particle[id]->x_center -
    Box->Lx*rint(Particles.N_Particle[id]->x_center/Box->Lx);
    Particles.N_Particle[id]->y_center = Particles.N_Particle[id]->y_center -
    Box->Lx*rint(Particles.N_Particle[id]->y_center/Box->Lx);

    Particles.N_Particle[id]->x_center = Particles.N_Particle[id]->x_center +
    Box->Lx/2.; Particles.N_Particle[id]->y_center =
    Particles.N_Particle[id]->y_center + Box->Lx/2.;
    */

    // be careful only sets gears to right orientation! Need to fix the other
    // particles as well. To do other particles uncomment the part in
    // Trans_Update_Positions and
    // Particles.N_Particle[id]->edges_from_center()!
    // Particles.N_Particle[id]->edges_from_center();

    Particles.N_Particle[id]->Calculate_Axis();
    Particles.N_Particle[id]->Calculate_Patch_Position();

    Particles.Update_Cell_List(id, Box);

    Particles.Collision_List[id].Calculate(Box, id, Particles.Id_Cell_List,
                                            Particles.Cell_List, Particles.Cell,
                                            Particles.N_Particle,Particles.N_Particle[0]->cut_off,
                                            Particles.MAX_coll_p);
     //Particles.Collision_List[id].Calculate_OP(Box, id, Particles.N_Particle,
     // Particles.N_Particle[0]->cut_off,
     //                                       Particles.MAX_coll_p);


    exit_status = 0;
    col_count = 0;

    Collision_Test(Particles, Box, id, Particles.Collision_List);

    if (exit_status >= 1) {

        Reset_Positions(Particles, id);

        Particles.N_Particle[id]->Calculate_Axis();
        Particles.N_Particle[id]->Calculate_Patch_Position();

        Particles.Reset_Cell_List(Box, id, Particles.c_id, Particles.n_id,
        Particles.id_num);
    }

    ////cout<<"x_center "<<Particles.N_Particle[id]->x_center<<"  "<<id<<endl;

    if (exit_status == 0) {

        // Set_Pair_Potential(Particles, Box);

        dU_new = Calculate_Pair_Potential(id, Particles, Box,
                                          Particles.Collision_List);
        dU_new = dU_new + Halo_Energy;
        // cout<<" New Halo_Energy "<<Halo_Energy<<endl;
        // cout<<"dU_new "<<dU_new<<endl;

        delta_U = dU_new - dU_old;

        // Calculate_Pair_Potential(Particles, Box);

        // MC Step for Potential
        b_factor_pre = exp(-1.0 * beta * delta_U);
        // b_factor_pre = exp(-1.0*beta*(Total_Energy-Total_Energy_old));
        b_factor = minimum(1, b_factor_pre);
        XI = gsl_rng_uniform(r01);

        /*
        ////cout<<"TRANSLATE"<<endl;
        ////cout<<"Total E new...."<<Total_Energy<<endl;
        ////cout<<"Total E old...."<<Total_Energy_old<<endl;
        ////cout<<"Delta E........"<<Total_Energy-Total_Energy_old<<endl;
        ////cout<<"b_factor......."<<b_factor<<endl;
        ////cout<<"XI............."<<XI<<endl;
       */

        // Reject of XI > b_factor

        if (XI > b_factor) {

            //////cout<<"Rejected!"<<endl;

            Reset_Positions(Particles, id);

            Particles.N_Particle[id]->Calculate_Axis();
            Particles.N_Particle[id]->Calculate_Patch_Position();

            // Reset_Pair_Potential(Total_Energy, id, Particles, Box);
            // Reset_Pair_Potential(Particles, Box);
            //////cout<<"Total_Energy after reject"<<Total_Energy<<endl;
        }

        // Accept move if XI <= b_factor

        if (XI <= b_factor) {

            //////cout<<"Accepted!"<<endl;

            Set_Positions(Particles, id);
            // Set_Pair_Potential(Total_Energy, id, Particles, Box);
            // Set_Pair_Potential(Particles, Box);
            //////cout<<"Total_Energy after accept"<<Total_Energy<<endl;
            // Particles.Total_Energy = Total_Energy;
            Particles.Total_Energy = Particles.Total_Energy - dU_old + dU_new;
            // cout<<" dU_old "<<dU_old<<endl;
            // cout<<" dU_new "<<dU_new<<endl;
            // cout<<"Total_Energy after accept"<<Particles.Total_Energy<<endl;
            accept_translate = accept_translate + 1;
            N_trans_moves = N_trans_moves + 1;
        }
    }

    ////cout<<"x_center "<<Particles.N_Particle[id]->x_center<<"  "<<id<<endl;
}
