#include "main.h"

int main(int argc, char *argv[]) {

    r = gsl_rng_alloc(gsl_rng_ranlxd2);
    r01 = gsl_rng_alloc(gsl_rng_ranlxd2);

    strcpy(runtype, argv[1]);

    strcpy(new_runtype, "-n");
    strcpy(former_runtype, "-f");

    checkpoint_time = atoi(argv[2]);
    start_time = checkpoint_time;

    term_time = atoi(argv[3]);

    wv_0_in = 0;
    wv_i_in = 1;
    wv_ip1_in = 2;

    fileio Fileio;
    Fileio.Set_Vars();
    Fileio.Read_Parameters();

    calc_frequency = Fileio.calc_frequency_in;
    frame_frequency = Fileio.frame_frequency_in;
    checkpoint_frequency = Fileio.checkpoint_frequency_in;

    box *Box;
    Box = new box;

    strcpy(traj_type_in, "forward");

    particles Particles(Fileio.number_of_cells_in, Fileio.N_in,
                        Fileio.MAX_coll_p_in, Fileio.MAX_fshell_p_in,
                        Fileio.particle_type_in, Fileio.binary_on_in,
                        Fileio.phi_binary_in);

    // init new run
    if (strcmp(runtype, new_runtype) == 0) {

        // initialize box
        Box->Startconfig(Fileio.N_in, Fileio.P_sigma_in,
            Fileio.mu_in_1, Fileio.mu_in_2,
                         Fileio.packing_fraction_in,
                         Particles.N_Particle[0]->A);
        start_time = 0;
        Box->edges_from_center();
        // initialize lattice
        Particles.Startconfig(Box);
        // Particles.Make_Cell_List(Box);
        // Particles.Set_Cell_List(Box);
        // Particles.Make_Cell_Neighbour_List();

        gsl_rng_set(r01, 0xf143);
        gsl_rng_set(r, 0x01a23);
    }

    if (strcmp(runtype, former_runtype) == 0) {

        // initialize former box
        // Box->Startconfig(Fileio.N_in, Fileio.P_sigma_in, Fileio.mu_in,
        // Fileio.packing_fraction_in, Particles.N_Particle[0]->V);
        Fileio.Read_Box(Box, checkpoint_time);
        Box->Startconfig_former(Fileio.N_in, Fileio.P_sigma_in,
            Fileio.mu_in_1, Fileio.mu_in_2);
        Box->edges_from_center();

        // random state initialize
        // Fileio.Read_Random_State(r,r01, checkpoint_time);

        int R1, R2;

        R1 = atoi(argv[4]);
        R2 = atoi(argv[5]);

        gsl_rng_set(r01, R1);
        gsl_rng_set(r, R2);

        // initialize former positions

        Particles.Startconfig(Box);

        Fileio.Read_Positions(Box, Particles, checkpoint_time);
        Fileio.Read_Orientations(Box, Particles, checkpoint_time);
        Fileio.Read_Patches(Box, Particles, checkpoint_time);
        Particles.Set_former_Config(Box);

        // Particles.Make_Cell_List(Box);
        // Particles.Set_Cell_List(Box);
        // Particles.Make_Cell_Neighbour_List();
    }

    int time;

    time = start_time;

    order_parameter Order_Parameter(Box, Particles, Fileio.MAX_coll_p_in);
    cluster Cluster(Box);

    pmove Move(Particles.Res_Size, Particles.edge_N, Fileio.delta_tmax,
               Fileio.delta_rmax, Fileio.delta_Vmax, Fileio.Temperature,
               Fileio.is_translation_ON, Fileio.is_rotation_ON,
               Fileio.is_volumemove_ON, Fileio.is_grand_canonical_ON,
               Fileio.is_cluster_move_ON, Fileio.is_2D);

    Move.Set_Random_State(r, r01);

    int MC_cycle_time;
    time = start_time;

    Fileio.Write_Positions(time, Box, Particles);

    Fileio.Write_Box_Positions(time, Box);

    double Total_Halo_Energy;
    Total_Halo_Energy = 0;
    cout << "Number of particles........." << Box->N << endl;
    cout << "Number resovoir particles..." << Particles.Res_Size << endl;
    cout << "Temperature................." << Move.T << endl;
    cout << "Chemical potential mu_1......." << Box->mu_1 << endl;
    cout << "Chemical potential mu_2......." << Box->mu_2 << endl;

    for (int id = 0; id < Box->N; id++) {
        Particles.Collision_List[id].Calculate_OP(
            Box, id, Particles.N_Particle, Particles.N_Particle[0]->cut_off,
            Particles.MAX_coll_p);
        Move.Collision_Test(Particles, Box, id, Particles.Collision_List);
        Total_Halo_Energy += Move.Halo_Energy;
    }

    Move.Calculate_Pair_Potential(Particles, Box);
    Move.Set_Pair_Potential(Particles, Box);

    Particles.Total_Energy = Move.Total_Energy + Total_Halo_Energy / 2.0;

    Box->P_sigma = Fileio.P_sigma_in;

    cout << "Packing_fraction: " << Box->packing_fraction << endl;

    do {

        time = time + 1;

        MC_cycle_time = Move.mt_sum;

        for (int cycle_time = 0; cycle_time < MC_cycle_time; cycle_time++) {

            Move.Iterate(Particles, Box, Fileio, time);
        }

        if (time % frame_frequency == 0) {

            Fileio.Write_Positions(time, Box, Particles);
            Fileio.Write_Box_Positions(time, Box);
            // Fileio.Write_Positions_Gears(time, Box, Particles);
        }

        if (time % calc_frequency == 0) {

            cout << "calc..." << endl;
            Move.Calculate_Acceptances(time);

            // Order_Parameter.Calculate_Local_Order_Parameters(Particles,
            // Particles.N_Particle, Box);

            for (int id = 0; id < Box->N; id++) {
                Particles.Collision_List[id].Calculate_Bonds(
                    Box, id, Particles.N_Particle, Particles.MAX_coll_p);
            }

            for (int id = 0; id < Box->N; id++) {
                Particles.Collision_List[id].Calculate_Neighbours(
                    Particles.N_Particle[id]->cut_off);
            }

            Order_Parameter.Caculate_is_nPhase(Box);
            bool is_mass_calc;
            is_mass_calc = true;

            Cluster.Calculate(Order_Parameter, Particles, Box, time,
                              is_mass_calc);
            Order_Parameter.Calculate_2D_Psi(Cluster.Largest_Cluster_List,
                                             Cluster.Size_R, Particles, Box);
            Order_Parameter.Calculate_2D_Psi(Particles, Box);
            Order_Parameter.Calculate_Colors(Cluster.Largest_Cluster_List,
                                             Cluster.Size_R, Particles, Box);

            Fileio.Write_Acceptances(time, Move.accept_translate_procent,
                                     Move.accept_rotate_procent,
                                     Move.accept_iso_vol_procent,
                                     Move.accept_complete_procent);
            Fileio.Write_NPT(time, Box);
            Fileio.Write_Energy(Box, Particles, time);
            Fileio.Write_Local_Order_Parameters(time, Particles, Box,
                                                Order_Parameter, Cluster);
        }

        if ((time % checkpoint_frequency == 0) || (time == start_time + 1)) {
            cout << "Time...." << time << endl;
            Fileio.Save_Config(Box, Particles, Move.r, Move.r01, Move.dmax_t,
                               Move.dmax_alpha, Move.dmax_beta,
                               Move.dmax_gamma, Move.dmax_V, Move.dmax_L,
                               time);
            N_check_points = N_check_points + 1;
        }

        if ((exit_code == true) && (time % checkpoint_frequency == 0)) {
            term_code = true;
        }

    } while ((time < term_time) && (exit_code == false));

    Fileio.Save_Config(Box, Particles, Move.r, Move.r01, Move.dmax_t,
                       Move.dmax_alpha, Move.dmax_beta, Move.dmax_gamma,
                       Move.dmax_V, Move.dmax_L, term_time);

    cout << "End of Simulation. BYE BYE" << endl;
}
