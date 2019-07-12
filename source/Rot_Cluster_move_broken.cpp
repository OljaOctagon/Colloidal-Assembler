void pmove::Rot_Pseudocluster_Recursion(int id_j, int cn, int fl,
                                        particles &Particles, box *Box) {

    int j;

    if (is_element[id_j] == false) {
        List[N_List] = id_j;
        N_List = N_List + 1;
        is_element[id_j] = true;
    }

    if (N_List == 1) {
        r_id = id_j;
    }

    Particles.Collision_List[id_j].Calculate_OP(
        Box, id_j, Particles.N_Particle, Particles.N_Particle[0]->cut_off,
        Particles.MAX_coll_p);

    for (int cp = 0; cp < Particles.Collision_List[id_j].Nm; cp++) {
        Particles.Collision_List_old[id_j].Elements[cp] =
            Particles.Collision_List[id_j].Elements[cp];
    }

    Particles.Collision_List_old[id_j].Nm = Particles.Collision_List[id_j].Nm;

    Rot_Move_Map(Particles, r_id, r_id, Box, Rot_mat);
    Particles.Check_Periodic_CM(r_id, Box);
    Update_Periodic_Positions(Particles, Box, r_id);
    Particles.N_Particle[r_id]->Calculate_Axis();
    Particles.N_Particle[r_id]->Calculate_Patch_Position();

    Rot_Move_Map(Particles, r_id, id_j, Box, Rot_mat);
    Particles.Check_Periodic_CM(id_j, Box);
    Update_Periodic_Positions(Particles, Box, id_j);
    Particles.N_Particle[id_j]->Calculate_Axis();
    Particles.N_Particle[id_j]->Calculate_Patch_Position();

    Particles.Collision_List[id_j].Calculate_OP(
        Box, id_j, Particles.N_Particle, Particles.N_Particle[0]->cut_off,
        Particles.MAX_coll_p);

    Reset_Positions(Particles, r_id);
    Particles.N_Particle[r_id]->Calculate_Axis();
    Particles.N_Particle[r_id]->Calculate_Patch_Position();

    Reset_Positions(Particles, id_j);
    Particles.N_Particle[id_j]->Calculate_Axis();
    Particles.N_Particle[id_j]->Calculate_Patch_Position();

    for (int cp = 0; cp < Particles.Collision_List[id_j].Nm; cp++) {

        Rot_Move_Map(Particles, r_id, r_id, Box, Rot_mat);
        Particles.Check_Periodic_CM(r_id, Box);
        Update_Periodic_Positions(Particles, Box, r_id);
        Particles.N_Particle[r_id]->Calculate_Axis();
        Particles.N_Particle[r_id]->Calculate_Patch_Position();

        Rot_Move_Map(Particles, r_id, id_j, Box, Rot_mat);
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

            Rot_Move_Map(Particles, r_id, j, Box, Rot_mat);
            Particles.Check_Periodic_CM(j, Box);
            Update_Periodic_Positions(Particles, Box, j);

            Particles.N_Particle[j]->Calculate_Axis();
            Particles.N_Particle[j]->Calculate_Patch_Position();

            e_pair = Calculate_Potential(id_j, j, Particles, Box);

            q_factor = exp(beta * (e_pair - e_single));
            b_factor = 1 - exp(beta * (e_pair - e_single));

            cout << "e_pair " << e_pair << endl;
            cout << "e_single " << e_single << endl;
            cout << "b_factor " << b_factor << endl;

            b_factor = GSL_MAX(0, b_factor);
            q_factor = GSL_MIN(1, q_factor);

            XI = gsl_rng_uniform(r01);

            Reset_Positions(Particles, id_j);
            Particles.N_Particle[id_j]->Calculate_Axis();
            Particles.N_Particle[id_j]->Calculate_Patch_Position();

            Reset_Positions(Particles, r_id);
            Particles.N_Particle[r_id]->Calculate_Axis();
            Particles.N_Particle[r_id]->Calculate_Patch_Position();

            if (XI >= b_factor) {

                Reset_Positions(Particles, j);
                Particles.N_Particle[j]->Calculate_Axis();
                Particles.N_Particle[j]->Calculate_Patch_Position();

                if (is_element[j] == true) {

                    fl = N_failed_links;

                    q_f[fl] = q_factor;

                    Rot_Move_Map(Particles, r_id, r_id, Box, R_Rot_mat);
                    Particles.Check_Periodic_CM(r_id, Box);
                    Update_Periodic_Positions(Particles, Box, r_id);

                    Particles.N_Particle[r_id]->Calculate_Axis();
                    Particles.N_Particle[r_id]->Calculate_Patch_Position();

                    Rot_Move_Map(Particles, r_id, id_j, Box, R_Rot_mat);
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

                    Reset_Positions(Particles, r_id);
                    Particles.N_Particle[r_id]->Calculate_Axis();
                    Particles.N_Particle[r_id]->Calculate_Patch_Position();

                    N_failed_links = N_failed_links + 1;
                }
            }

            if (XI < b_factor) {

                // Need to distinguish number of bonds and number of particles
                // that go into list

                cn = N_Bonds;
                p_f[cn] = b_factor;
                //////cout<<"b_factor "<<b_factor<<endl;
                // VIRTUAL opposite move

                // go with linker particle in oppositie direction

                // Rot_Update_Positions(Particles, id_j, R_Rot_mat);

                Rot_Move_Map(Particles, r_id, r_id, Box, R_Rot_mat);
                Particles.Check_Periodic_CM(r_id, Box);
                Update_Periodic_Positions(Particles, Box, r_id);

                Particles.N_Particle[r_id]->Calculate_Axis();
                Particles.N_Particle[r_id]->Calculate_Patch_Position();

                Rot_Move_Map(Particles, r_id, id_j, Box, R_Rot_mat);
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

                cout << " e_pair " << e_pair << endl;
                cout << " e_single " << e_single << endl;

                b_factor = 1 - exp(beta * (e_pair - e_single));
                b_factor = GSL_MAX(0, b_factor);
                p_r[cn] = b_factor;

                cout << "p_r" << p_r[cn] << endl;

                //////cout<<" b_factor reverse "<<b_factor<<endl;
                // move linker back to new move

                Reset_Positions(Particles, r_id);
                Particles.N_Particle[r_id]->Calculate_Axis();
                Particles.N_Particle[r_id]->Calculate_Patch_Position();

                Reset_Positions(Particles, id_j);
                Particles.N_Particle[id_j]->Calculate_Axis();
                Particles.N_Particle[id_j]->Calculate_Patch_Position();

                N_Bonds = N_Bonds + 1;
                cout << "bond formed! b_factor: " << p_f[cn] << endl;
                // move linkee back to new move

                Rot_Pseudocluster_Recursion(j, cn, fl, Particles, Box);
            }
        }
    }
}

void pmove::Rot_Cluster_Move(particles &Particles, box *Box, fileio &Fileio,
                             int mc_time) {

    cout << "move" << endl;

    N_List = 0;
    N_Bonds = 0;
    r_id = -1;
    double phi_t;
    double phi_r;

    // choose particle id

    id = gsl_rng_uniform_int(r, Box->N);
    //////cout<<"id "<<id<<endl;

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

    phi_r = -phi_t;

    R_Rot_mat[0] = cos(phi_r);
    R_Rot_mat[1] = -sin(phi_r);
    R_Rot_mat[2] = 0;

    R_Rot_mat[3] = sin(phi_r);
    R_Rot_mat[4] = cos(phi_r);
    R_Rot_mat[5] = 0;

    R_Rot_mat[6] = 0;
    R_Rot_mat[7] = 0;
    R_Rot_mat[8] = 1;

    int links, failed_links;
    links = 0;
    failed_links = 0;

    Rot_Pseudocluster_Recursion(id, links, failed_links, Particles, Box);

    // end up with List[k] list of elements in Pseudocluster
    // linking and reverse lining probabilities

    // cout<<"N_List "<<N_List<<endl;
    // cout<<"N_Bonds "<<N_Bonds<<endl;

    // Collision Test

    for (int k = 0; k < N_List; k++) {
        int j;
        j = List[k];

        Reset_Positions(Particles, j);
        Rot_Move_Map(Particles, id, j, Box, Rot_mat);
        Particles.Check_Periodic_CM(j, Box);
        Update_Periodic_Positions(Particles, Box, j);
        Particles.N_Particle[j]->Calculate_Axis();
        Particles.N_Particle[j]->Calculate_Patch_Position();
    }

    int k = 0;
    //////cout<<"N_List   "<<N_List<<endl;
    //////cout<<"N_Bonds  "<<N_Bonds<<endl;

    exit_status = 0;
    col_count = 0;

    do {
        Particles.Collision_List[List[k]].Calculate_OP(
            Box, List[k], Particles.N_Particle,
            Particles.N_Particle[0]->cut_off, Particles.MAX_coll_p);
        // Collision Test
        Collision_Test(Particles, Box, List[k], Particles.Collision_List);
        k++;

    } while ((exit_status == 0) && (k < N_List));

    if (exit_status > 0) {

        //////cout<<"collsion!"<<endl;
        // Reset all positions!
        for (int k = 0; k < N_List; k++) {
            int j;
            j = List[k];
            Reset_Positions(Particles, j);
            Particles.N_Particle[j]->Calculate_Axis();
            Particles.N_Particle[j]->Calculate_Patch_Position();
            Set_Positions(Particles, j);
        }
    }

    if (exit_status == 0) {

        cout << "exit status 0 " << endl;

        int list_j;
        int list_k;
        int m;
        b_factor = 1;
        bool overlapped;

        //////cout<<"no collision"<<endl;
        // Calculate new factors that are interacting in the new step

        cout << "N_List " << N_List << endl;
        cout << "N_Bonds " << N_Bonds << endl;

        for (int k = 0; k < N_List; k++) {

            m = List[k];

            for (int cp1 = 0; cp1 < Particles.Collision_List[m].Nm; cp1++) {
                overlapped = false;
                list_j = Particles.Collision_List[m].Elements[cp1].nl_id;

                for (int cp2 = 0; cp2 < Particles.Collision_List_old[m].Nm;
                     cp2++) {
                    list_k =
                        Particles.Collision_List_old[m].Elements[cp2].nl_id;

                    if (list_j == list_k) {
                        overlapped = true;
                    }
                }

                if (overlapped == false) {
                    // Calculate boltzmann factor

                    // only new potential plays role

                    e_new = Calculate_Potential(m, list_j, Particles, Box);
                    b_factor = b_factor * exp(-beta * e_new);
                }
            }
        }

        //////cout<<"b_factor after new add "<<b_factor<<endl;

        // Calculate factors for particels that were interacting in the last
        // time step

        for (int k = 0; k < N_List; k++) {

            m = List[k];

            for (int cp2 = 0; cp2 < Particles.Collision_List_old[m].Nm;
                 cp2++) {
                overlapped = false;
                list_k = Particles.Collision_List_old[m].Elements[cp2].nl_id;
                for (int cp1 = 0; cp1 < Particles.Collision_List[m].Nm;
                     cp1++) {
                    list_j = Particles.Collision_List[m].Elements[cp1].nl_id;

                    if (list_j == list_k) {
                        overlapped = true;
                    }
                }
                if (overlapped == false) {
                    // Calculate boltzmann factor

                    // only old potential plays a role
                    Reset_Positions(Particles, m);
                    Particles.N_Particle[m]->Calculate_Axis();
                    Particles.N_Particle[m]->Calculate_Patch_Position();

                    e_old = Calculate_Potential(m, list_k, Particles, Box);
                    b_factor = b_factor * exp(beta * e_old);

                    // move back

                    Rot_Move_Map(Particles, id, m, Box, Rot_mat);
                    Particles.Check_Periodic_CM(m, Box);
                    Update_Periodic_Positions(Particles, Box, m);
                    Particles.N_Particle[m]->Calculate_Axis();
                    Particles.N_Particle[m]->Calculate_Patch_Position();
                }
            }
        }

        // Add up all the factors
        double p_factor = 1;
        double s_factor = 1;

        double p_fr;
        double q_fr;
        int l = 0;
        int k = 0;

        do {
            cout << " p_r[k] " << p_r[k] << endl;
            cout << " p_f[k]" << p_f[k] << endl;

            if (p_r[k] < 1e-10) {
                p_fr = 0;

            }

            else {
                p_fr = double(p_r[k]) / double(p_f[k]);
            }

            p_factor = p_factor * p_fr;

            k++;

        } while (k < N_Bonds);

        // cout<<"Boltzmann factor for Pseudocluster "<<p_factor<<endl;

        do {

            if (q_r[l] < 1e-10) {
                q_fr = 0;
            }

            else {
                q_fr = double(q_r[l]) / double(q_f[l]);
            }

            s_factor = s_factor * q_fr;
            l++;

        } while (l < N_failed_links);

        //////cout<<"p_factor "<<p_factor<<endl;
        cout << "b_factor" << b_factor << endl;
        cout << "s_factor" << s_factor << endl;
        cout << "p_factor" << p_factor << endl;

        b_factor = GSL_MIN(1, b_factor * s_factor * p_factor);
        XI = gsl_rng_uniform(r01);

        // Cluster Move Acceptance criterium

        if (b_factor < XI) {

            // Reset all positions
            //////cout<<"reset cluster move"<<endl;
            //////cout<<"Number of not moved particles "<<N_List<<endl;
            for (int k = 0; k < N_List; k++) {
                int j;
                j = List[k];
                Reset_Positions(Particles, j);
                Particles.N_Particle[j]->Calculate_Axis();
                Particles.N_Particle[j]->Calculate_Patch_Position();
            }
        }

        if (b_factor >= XI) {

            cout << "accepted" << endl;

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

            Calculate_Pair_Potential(Particles, Box);
        }
    }

    Reset_Pseudo_Cluster(Box);
}
