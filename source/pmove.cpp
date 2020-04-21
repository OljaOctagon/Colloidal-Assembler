
#include "pmove.h"

pmove::pmove() {}

double pmove::minimum(double a, double b) {

    if (a <= b) {
        return_value = a;
    }

    if (b < a) {
        return_value = b;
    }

    return return_value;
}

pmove::pmove(int size, int edge_length, double delta_tmax, double delta_rmax,
             double delta_Vmax, double Temperature,
             string is_translation_ON_in, string is_rotation_ON_in,
             string is_volumemove_ON_in, string is_grand_canonical_ON_in,
             string is_cluster_move_ON_in, int is_2D_in) {

    List = new int[size];
    p_f = new double[size];
    p_r = new double[size];

    q_r = new double[size];
    q_f = new double[size];

    pseudo_cluster_info = new int *[size];
    for (int j = 0; j < size; j++) {
        pseudo_cluster_info[j] = new int[size];
    }

    Patch_Energies = new double *[size];
    for (int j = 0; j < size; j++) {
        Patch_Energies[j] = new double[size];
    }

    Patch_Energies_old = new double *[size];
    for (int j = 0; j < size; j++) {
        Patch_Energies_old[j] = new double[size];
    }

    Total_Energy = 0;
    Total_Energy_old = 0;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            Patch_Energies[i][j] = 0.0;
        }
    }

    is_element = new bool[size];
    for (int i = 0; i < size; i++) {
        is_element[i] = false;
    }

    L_axis = new m_vector[15];
    Cluster_List = new int[size];

    for (int j = 0; j < 15; j++) {

        L_axis[j].x = 0;
        L_axis[j].y = 0;
        L_axis[j].z = 0;
    }

    is_2D = is_2D_in;

    accept_translate = 0;
    accept_rotate = 0;
    accept_iso_vol = 0;
    accept_aniso_vol = 0;

    map_dist = new m_vector[8];

    for (int i = 0; i < 8; i++) {

        map_dist[i].x = 0.0;
        map_dist[i].y = 0.0;
        map_dist[i].z = 0.0;
    }

    N_trans_moves = 0;
    N_rotate_moves = 0;
    N_iso_vol_moves = 0;
    N_aniso_vol_moves = 0;

    sigma_trans = delta_tmax;
    dmax_V = delta_Vmax;
    kappa = delta_rmax;

    // cout<<"kappa"<<kappa<<endl;

    Unity[0] = 1.0;
    Unity[1] = 0.0;
    Unity[2] = 0.0;

    Unity[3] = 0.0;
    Unity[4] = 1.0;
    Unity[5] = 0.0;

    Unity[6] = 0.0;
    Unity[7] = 0.0;
    Unity[8] = 1.0;

    zero_vec.x = 0.0;
    zero_vec.y = 0.0;
    zero_vec.z = 0.0;

    phi_reverse = 0;

    col_count = 0;

    r = gsl_rng_alloc(gsl_rng_ranlxd2);
    r01 = gsl_rng_alloc(gsl_rng_ranlxd2);

    for (int j = 0; j < 9; j++) {
        Rot_mat[j] = Unity[j];
    }

    q_0 = 1.0;
    q_1 = 0.0;
    q_2 = 0.0, q_3 = 0.0;

    T = Temperature;
    // kB = 1.3806488*pow(10.0,-23.0);
    kB = 1.0;
    beta = 1. / (kB * T);
    beta_f = beta/2.;
    

    // init discrete distribution for move_type
    P_mt = new double[5];

    is_translation_ON = is_translation_ON_in;
    is_rotation_ON = is_rotation_ON_in;
    is_volumemove_ON = is_volumemove_ON_in;
    is_grand_canonical_ON = is_grand_canonical_ON_in;
    is_cluster_move_ON = is_cluster_move_ON_in;

    mt_w[0] = 0;
    mt_w[1] = 0;
    mt_w[2] = 0;
    mt_w[3] = 0;
    mt_w[4] = 0;

    mt_sum = 0;

    cout << "translation move...... " << is_translation_ON << endl;
    cout << "rotation move.........." << is_rotation_ON << endl;
    cout << "volume move......." << is_volumemove_ON << endl;
    cout << "grand canonical move..." << is_grand_canonical_ON << endl;
    cout << "cluster move..........." << is_cluster_move_ON << endl;
    cout << "is_2D..................." << is_2D << endl;

    if (is_translation_ON.compare("on") == 0) {
        mt_w[0] = size;
    }

    if (is_rotation_ON.compare("on") == 0) {
        mt_w[1] = size;
    }

    if (is_volumemove_ON.compare("on") == 0) {
        mt_w[2] = 1;
    }

    if (is_grand_canonical_ON.compare("on") == 0) {
        mt_w[3] = 1;
    }

    if (is_cluster_move_ON.compare("on") == 0) {
        mt_w[4] = 3;
    }

    mt_sum = mt_w[0] + mt_w[1] + mt_w[2] + mt_w[3] + mt_w[4];
    cout<<"mt_sum: "<<mt_sum<<endl;
    cout<<"size: "<<size<<endl;

    P_mt[0] = double(mt_w[0]) / double(mt_sum);
    P_mt[1] = double(mt_w[1]) / double(mt_sum);
    P_mt[2] = double(mt_w[2]) / double(mt_sum);
    P_mt[3] = double(mt_w[3]) / double(mt_sum);
    P_mt[4] = double(mt_w[4]) / double(mt_sum);

    sizep = 5;

    mt = gsl_ran_discrete_preproc(sizep, P_mt);
}

pmove::~pmove() {
    delete[] L_axis;
    delete[] map_dist;
    delete[] is_element;
    // delete [] Cluster_List;
    // delete [] Patch_Energies;
    // delete [] Patch_Energies_old;
    // delete [] is_element;
}

void pmove::Set_Random_State(gsl_rng *r_start, gsl_rng *r01_start) {
    random_init_r = gsl_rng_memcpy(r, r_start);
    random_init_r01 = gsl_rng_memcpy(r01, r01_start);
}

void pmove::Trans_Update_Positions(particles &Particles, int id,
                                   m_vector &trans_vec) {

    Particles.N_Particle[id]->x_center =
        Particles.N_Particle[id]->x_center + trans_vec.x;
    Particles.N_Particle[id]->y_center =
        Particles.N_Particle[id]->y_center + trans_vec.y;
    Particles.N_Particle[id]->z_center =
        Particles.N_Particle[id]->z_center + trans_vec.z;

    // update positions

    // additional for normal particles

    for (int j = 0; j < Particles.N_Particle[id]->edge_N; j++) {

        Particles.N_Particle[id]->x[j] =
            Particles.N_Particle[id]->x[j] + trans_vec.x;
        Particles.N_Particle[id]->y[j] =
            Particles.N_Particle[id]->y[j] + trans_vec.y;
        Particles.N_Particle[id]->z[j] =
            Particles.N_Particle[id]->z[j] + trans_vec.z;
    }
}

void pmove::Rot_Update_Positions(particles &Particles, int id,
                                 double (&Rot_mat)[9]) {

    // for gear particles

    Particles.N_Particle[id]->edges_from_center();
    Particles.N_Particle[id]->distance_from_center();

    for (int j = 0; j < Particles.N_Particle[id]->edge_N; j++) {

        Particles.N_Particle[id]->new_dist_x[j] =
            Rot_mat[0] * Particles.N_Particle[id]->dist_x[j] +
            Rot_mat[1] * Particles.N_Particle[id]->dist_y[j] +
            Rot_mat[2] * Particles.N_Particle[id]->dist_z[j];
        Particles.N_Particle[id]->new_dist_y[j] =
            Rot_mat[3] * Particles.N_Particle[id]->dist_x[j] +
            Rot_mat[4] * Particles.N_Particle[id]->dist_y[j] +
            Rot_mat[5] * Particles.N_Particle[id]->dist_z[j];
        Particles.N_Particle[id]->new_dist_z[j] =
            Rot_mat[6] * Particles.N_Particle[id]->dist_x[j] +
            Rot_mat[7] * Particles.N_Particle[id]->dist_y[j] +
            Rot_mat[8] * Particles.N_Particle[id]->dist_z[j];

        Particles.N_Particle[id]->x[j] =
            Particles.N_Particle[id]->x_center +
            Particles.N_Particle[id]->new_dist_x[j];
        Particles.N_Particle[id]->y[j] =
            Particles.N_Particle[id]->y_center +
            Particles.N_Particle[id]->new_dist_y[j];
        Particles.N_Particle[id]->z[j] =
            Particles.N_Particle[id]->z_center +
            Particles.N_Particle[id]->new_dist_z[j];
    }
}

void pmove::Rot_Move_Map(particles &Particles, int id1, int id2, box *Box,
                         double (&Rot_mat)[9]) {}

void pmove::Check_Periodic_Center_of_Mass(m_vector &center_mass, box *Box) {

    int left_count, right_count;
    int front_count, back_count;
    int top_count, bottom_count;
    int edge_out;

    left_count = 0;
    right_count = 0;
    front_count = 0;
    back_count = 0;
    top_count = 0;
    bottom_count = 0;

    double lx, ly, lz;

    edge_out = 0;

    if (center_mass.x > Box->x[1]) {
        right_count = 1;
        edge_out = 1;
    }

    if (center_mass.x < Box->x[0]) {
        left_count = 1;
        edge_out = 1;
    }

    if (center_mass.y > Box->y[3]) {
        back_count = 1;
        edge_out = 1;
    }
    if (center_mass.y < Box->y[0]) {
        front_count = 1;
        edge_out = 1;
    }

    if (center_mass.z > Box->z[4]) {
        top_count = 1;
        edge_out = 1;
    }

    if (center_mass.z < Box->z[0]) {
        bottom_count = 1;
        edge_out = 1;
    }

    if (edge_out > 0) {
        lx = double(left_count) * Box->Lx - double(right_count) * Box->Lx;
        ly = double(front_count) * Box->Ly - double(back_count) * Box->Ly;
        lz = double(bottom_count) * Box->Lz - double(top_count) * Box->Lz;

        center_mass.x = center_mass.x + lx;
        center_mass.y = center_mass.y + ly;
        center_mass.z = center_mass.z + lz;
    }
}

void pmove::Rot_Move_Map(particles &Particles, int id1, box *Box,
                         m_vector &center_mass, double (&Rot_mat)[9]) {

    double particle_dist_x, particle_dist_y, particle_dist_z;

    Particles.N_Particle[id1]->distance_from_center();

    particle_dist_x = Particles.N_Particle[id1]->x_center - center_mass.x;
    particle_dist_y = Particles.N_Particle[id1]->y_center - center_mass.y;
    particle_dist_z = Particles.N_Particle[id1]->z_center - center_mass.z;

    particle_dist_x =
        particle_dist_x + Box->Lx * rint(particle_dist_x / Box->Lx);
    particle_dist_y =
        particle_dist_y + Box->Ly * rint(particle_dist_y / Box->Ly);
    particle_dist_z =
        particle_dist_z + Box->Lz * rint(particle_dist_z / Box->Lz);

    map_dist_center.x = particle_dist_x;
    map_dist_center.y = particle_dist_y;
    map_dist_center.z = particle_dist_z;

    map_rot_center.x = Rot_mat[0] * map_dist_center.x +
                       Rot_mat[1] * map_dist_center.y +
                       Rot_mat[2] * map_dist_center.z;
    map_rot_center.y = Rot_mat[3] * map_dist_center.x +
                       Rot_mat[4] * map_dist_center.y +
                       Rot_mat[5] * map_dist_center.z;
    map_rot_center.z = Rot_mat[6] * map_dist_center.x +
                       Rot_mat[7] * map_dist_center.y +
                       Rot_mat[8] * map_dist_center.z;

    Particles.N_Particle[id1]->x_center = center_mass.x + map_rot_center.x;
    Particles.N_Particle[id1]->y_center = center_mass.y + map_rot_center.y;
    Particles.N_Particle[id1]->z_center = center_mass.z + map_rot_center.z;

    for (int j = 0; j < Particles.N_Particle[id1]->edge_N; j++) {

        particle_dist_x = Particles.N_Particle[id1]->x[j] - center_mass.x;
        particle_dist_y = Particles.N_Particle[id1]->y[j] - center_mass.y;
        particle_dist_z = Particles.N_Particle[id1]->z[j] - center_mass.z;

        particle_dist_x =
            particle_dist_x + Box->Lx * rint(particle_dist_x / Box->Lx);
        particle_dist_y =
            particle_dist_y + Box->Ly * rint(particle_dist_y / Box->Ly);
        particle_dist_z =
            particle_dist_z + Box->Lz * rint(particle_dist_z / Box->Lz);

        map_dist[j].x = particle_dist_x;
        map_dist[j].y = particle_dist_y;
        map_dist[j].z = particle_dist_z;
    }

    for (int j = 0; j < Particles.N_Particle[id1]->edge_N; j++) {

        Particles.N_Particle[id1]->new_dist_x[j] = Rot_mat[0] * map_dist[j].x +
                                                   Rot_mat[1] * map_dist[j].y +
                                                   Rot_mat[2] * map_dist[j].z;
        Particles.N_Particle[id1]->new_dist_y[j] = Rot_mat[3] * map_dist[j].x +
                                                   Rot_mat[4] * map_dist[j].y +
                                                   Rot_mat[5] * map_dist[j].z;
        Particles.N_Particle[id1]->new_dist_z[j] = Rot_mat[6] * map_dist[j].x +
                                                   Rot_mat[7] * map_dist[j].y +
                                                   Rot_mat[8] * map_dist[j].z;

        Particles.N_Particle[id1]->x[j] =
            center_mass.x + Particles.N_Particle[id1]->new_dist_x[j];
        Particles.N_Particle[id1]->y[j] =
            center_mass.y + Particles.N_Particle[id1]->new_dist_y[j];
        Particles.N_Particle[id1]->z[j] =
            center_mass.z + Particles.N_Particle[id1]->new_dist_z[j];
    }
}

int pmove::Random(int N) {
    id = gsl_rng_uniform_int(r, N);
    return id;
}

void pmove::Iterate(particles &Particles, box *Box, fileio &Fileio,
                    int mc_time) {

    if (is_translation_ON.compare("on") == 0) {
        mt_w[0] = Box->N;
    }

    if (is_rotation_ON.compare("on") == 0) {
        mt_w[1] = Box->N;
    }

    if (is_volumemove_ON.compare("on") == 0) {
        mt_w[2] = 1;
    }

    if (is_grand_canonical_ON.compare("on") == 0) {
        mt_w[3] = 1;
    }

    if (is_cluster_move_ON.compare("on") == 0) {
        mt_w[4] = 3;
    }

    mt_sum = mt_w[0] + mt_w[1] + mt_w[2] + mt_w[3] + mt_w[4];
    // cout<<"mt_sum: "<<mt_sum<<endl;
    // cout<<"size: "<<size<<endl;

    P_mt[0] = double(mt_w[0]) / double(mt_sum);
    P_mt[1] = double(mt_w[1]) / double(mt_sum);
    P_mt[2] = double(mt_w[2]) / double(mt_sum);
    P_mt[3] = double(mt_w[3]) / double(mt_sum);
    P_mt[4] = double(mt_w[4]) / double(mt_sum);

    mt = gsl_ran_discrete_preproc(sizep, P_mt);

    value = gsl_ran_discrete(r, mt);
    gsl_ran_discrete_free(mt);

    if (value == 0) {

        // cout<<"translate"<<endl;
        // cluster_counter=gsl_rng_uniform_int(r,cluster_size);
        // id= Cluster_List[cluster_counter];

        id = gsl_rng_uniform_int(r, Box->N);
        Translate(Particles, Box, Fileio, id, mc_time);
    }

    if (value == 1) {

        // cluster_counter=gsl_rng_uniform_int(r,cluster_size);
        // id= Cluster_List[cluster_counter];
        // cout<<"rotate"<<endl;

        id = gsl_rng_uniform_int(r, Box->N);

        if (is_2D == 1) {

            Rotate2D(Particles, Box, Fileio, id, mc_time);
        }

        else {

            Rotate(Particles, Box, Fileio, id, mc_time);
        }
    }

    if (value == 2) {

        Iso_Vol_Change(Particles, Box, Fileio, mc_time);
        // Aniso_Vol_Change (Box, Particles.N_Particle,
        // Particles.N_Particle_old, Particles, Fileio, mc_time);
    }

    if (value == 3) {

        // cout<<"grand_can"<<endl;
        double rhalf;
        rhalf = gsl_rng_uniform(r01);

        if (rhalf > 0.5) {
            Particle_Insertion(Particles, Box, Fileio, mc_time);
        }

        else {
            Particle_Deletion(Particles, Box, Fileio, mc_time);
        }

        // cout<<"end grand can"<<endl;
    }

    if (value == 4) {

        double rhalf;
        rhalf = gsl_rng_uniform(r01);

        if (rhalf < 0.5) {
            Trans_Cluster_Move(Particles, Box, Fileio, mc_time);
        }

        if (rhalf > 0.5) {

            Rot_Cluster_Move(Particles, Box, Fileio, mc_time);
        }

        // if (rhalf>2./3.){
        // Three_Particle_Move(Particles, Box, Fileio, mc_time);
        //}
    }
}

void pmove::Calculate_Acceptances(int mc_time) {

    accept_translate_procent =
        double(accept_translate) / double(N_trans_moves);
    accept_rotate_procent = double(accept_rotate) / double(N_rotate_moves);
    accept_iso_vol_procent = double(accept_iso_vol) / double(N_iso_vol_moves);
    accept_complete_procent =
        double(accept_translate_procent + accept_rotate_procent +
               accept_iso_vol_procent) /
        double(3.0);
}

void pmove::Calculate_Cluster_List(particles &Particles, box *Box) {

    cluster_radius = Box->Lx / 3.0;
    cluster_radius_square = (Box->Lx * Box->Lx) / 9.0;

    cluster_counter = 0;

    for (int id = 0; id < Box->N; id++) {
        Cluster_List[id] = -100;
    }

    for (int id = 0; id < Box->N; id++) {

        diffc_x = Particles.N_Particle[id]->x_center - Box->x_center;
        diffc_y = Particles.N_Particle[id]->y_center - Box->y_center;
        diffc_z = Particles.N_Particle[id]->z_center - Box->z_center;

        diffc_square =
            diffc_x * diffc_x + diffc_y * diffc_y + diffc_z * diffc_z;

        if (diffc_square < cluster_radius_square) {

            Cluster_List[cluster_counter] = id;
            cluster_counter = cluster_counter + 1;
        }

        cluster_size = cluster_counter;
    }
}

void pmove::Update_Periodic_Positions(particles &Particles, box *Box, int id) {

    if (Particles.N_Particle[id]->cm_out >= 1) {

        Particles.N_Particle[id]->trans_periodic[0] =
            double(Particles.N_Particle[id]->cm_left_count) * Box->Lx -
            double(Particles.N_Particle[id]->cm_right_count) * Box->Lx;
        Particles.N_Particle[id]->trans_periodic[1] =
            double(Particles.N_Particle[id]->cm_front_count) * Box->Ly -
            double(Particles.N_Particle[id]->cm_back_count) * Box->Ly;
        Particles.N_Particle[id]->trans_periodic[2] =
            double(Particles.N_Particle[id]->cm_bottom_count) * Box->Lz -
            double(Particles.N_Particle[id]->cm_top_count) * Box->Lz;

        Particles.N_Particle[id]->x_center =
            Particles.N_Particle[id]->x_center +
            Particles.N_Particle[id]->trans_periodic[0];
        Particles.N_Particle[id]->y_center =
            Particles.N_Particle[id]->y_center +
            Particles.N_Particle[id]->trans_periodic[1];
        Particles.N_Particle[id]->z_center =
            Particles.N_Particle[id]->z_center +
            Particles.N_Particle[id]->trans_periodic[2];

        // update positions

        for (int j = 0; j < Particles.N_Particle[id]->edge_N; j++) {

            Particles.N_Particle[id]->x[j] =
                Particles.N_Particle[id]->x[j] +
                Particles.N_Particle[id]->trans_periodic[0];
            Particles.N_Particle[id]->y[j] =
                Particles.N_Particle[id]->y[j] +
                Particles.N_Particle[id]->trans_periodic[1];
            Particles.N_Particle[id]->z[j] =
                Particles.N_Particle[id]->z[j] +
                Particles.N_Particle[id]->trans_periodic[2];
        }
    }
}

void pmove::Reset_Positions(particles &Particles, int id) {

    Particles.N_Particle[id]->x_center =
        Particles.N_Particle_old[id]->x_center;
    Particles.N_Particle[id]->y_center =
        Particles.N_Particle_old[id]->y_center;
    Particles.N_Particle[id]->z_center =
        Particles.N_Particle_old[id]->z_center;

    // Particles.N_Particle[id]->edges_from_center();

    // comment in for gears

    for (int j = 0; j < Particles.N_Particle[id]->edge_N; j++) {

        Particles.N_Particle[id]->x[j] = Particles.N_Particle_old[id]->x[j];
        Particles.N_Particle[id]->y[j] = Particles.N_Particle_old[id]->y[j];
        Particles.N_Particle[id]->z[j] = Particles.N_Particle_old[id]->z[j];
    }

    Particles.N_Particle[id]->phi = Particles.N_Particle[id]->phi_old;

    Particles.N_Particle[id]->q.x = Particles.N_Particle_old[id]->q.x;
    Particles.N_Particle[id]->q.y = Particles.N_Particle_old[id]->q.y;
    Particles.N_Particle[id]->q.z = Particles.N_Particle_old[id]->q.z;
    Particles.N_Particle[id]->q.w = Particles.N_Particle_old[id]->q.w;
}

void pmove::Set_Positions(particles &Particles, int id) {

    Particles.N_Particle_old[id]->x_center =
        Particles.N_Particle[id]->x_center;
    Particles.N_Particle_old[id]->y_center =
        Particles.N_Particle[id]->y_center;
    Particles.N_Particle_old[id]->z_center =
        Particles.N_Particle[id]->z_center;

    // comment for gears

    for (int j = 0; j < Particles.N_Particle[id]->edge_N; j++) {

        Particles.N_Particle_old[id]->x[j] = Particles.N_Particle[id]->x[j];
        Particles.N_Particle_old[id]->y[j] = Particles.N_Particle[id]->y[j];
        Particles.N_Particle_old[id]->z[j] = Particles.N_Particle[id]->z[j];
    }

    Particles.N_Particle[id]->phi_old = Particles.N_Particle[id]->phi;

    Particles.N_Particle[id]->q.x = 0.0;
    Particles.N_Particle[id]->q.y = 0.0;
    Particles.N_Particle[id]->q.z =
        1.0 * sin(Particles.N_Particle[id]->phi / 2.0);
    Particles.N_Particle[id]->q.w = cos(Particles.N_Particle[id]->phi / 2.0);

    Particles.N_Particle_old[id]->q.x = Particles.N_Particle[id]->q.x;
    Particles.N_Particle_old[id]->q.y = Particles.N_Particle[id]->q.y;
    Particles.N_Particle_old[id]->q.z = Particles.N_Particle[id]->q.z;
    Particles.N_Particle_old[id]->q.w = Particles.N_Particle[id]->q.w;
}

void pmove::Rot_Update_Quarternions_VON_MISES(particles &Particles, int id) {

    mu.w = Particles.N_Particle[id]->q.w;
    mu.x = Particles.N_Particle[id]->q.x;
    mu.y = Particles.N_Particle[id]->q.y;
    mu.z = Particles.N_Particle[id]->q.z;

    // for max displacement: 0.2: kappa = 25
    // for max displacement: 0.07 kappa = 204
    // for max displacement: 2.0, kappa = 0.25
    // for max displacement: 1.0 kappa = 1.0

    b = -kappa + sqrt(kappa * kappa + 1);
    x0 = (1.0 - b) / (1.0 + b);
    c = kappa * x0 + 2 * log(1 - x0 * x0);

    do {

        // Generate a random value with the beta(a, a) distribution with a
        // = 1.5.
        //(1.5 because the parameters for beta are (m - 1)/2 and m = 4.

        do {
            u = gsl_ran_flat(r01, -1.0, 1.0);
            v = gsl_rng_uniform(r01);

            s = u * u + v * v;

        } while (s > 1.0);

        z = 0.5 + u * v * sqrt(1 - s) / double(s);

        u = gsl_rng_uniform(r01);
        w = (1.0 - (1.0 + b) * z) / (1.0 - (1.0 - b) * z);
        t = kappa * w + 2.0 * log(1.0 - x0 * w) - c;

    } while (t < log(u));

    do {
        u0 = gsl_ran_flat(r01, -1.0, 1.0);
        u1 = gsl_ran_flat(r01, -1.0, 1.0);

        sq_u01 = u0 * u0 + u1 * u1;

    } while (sq_u01 > 1.0);

    sp_x = 2. * u0 * sqrt(1.0 - sq_u01);
    sp_y = 2. * u1 * sqrt(1.0 - sq_u01);
    sp_z = 1.0 - 2. * sq_u01;

    q_t.w = w;
    q_t.x = sp_x * sqrt(1.0 - w * w);
    q_t.y = sp_y * sqrt(1.0 - w * w);
    q_t.z = sp_z * sqrt(1.0 - w * w);

    q_s.w = q_t.w * mu.w - q_t.x * mu.x - q_t.y * mu.y - q_t.z * mu.z;
    q_s.x = q_t.w * mu.x + q_t.x * mu.w + q_t.y * mu.z - q_t.z * mu.y;
    q_s.y = q_t.w * mu.y - q_t.x * mu.z + q_t.y * mu.w + q_t.z * mu.x;
    q_s.z = q_t.w * mu.z + q_t.x * mu.y - q_t.y * mu.x + q_t.z * mu.w;

    Particles.N_Particle[id]->q.w = q_s.w;
    Particles.N_Particle[id]->q.x = q_s.x;
    Particles.N_Particle[id]->q.y = q_s.y;
    Particles.N_Particle[id]->q.z = q_s.z;

    double quat_norm;

    quat_norm = q_s.w * q_s.w + q_s.x * q_s.x + q_s.y * q_s.y + q_s.z * q_s.z;
    quat_norm = sqrt(quat_norm);

    /*
    ofstream orient_out("quarternion_von_mises.dat",  ios::out | ios::app);
    orient_out<<Particles.N_Particle[id]->q.w<<"
    "<<Particles.N_Particle[id]->q.x<<"  "<<Particles.N_Particle[id]->q.y<<"
    "<<Particles.N_Particle[id]->q.z<<"  "<<quat_norm<<endl;
    orient_out.close();
    */
}

void pmove::Rot_Update_Quarternions_RANDOM(particles &Particles, int id) {

    double scalar_q;

    // do {

    do {
        u0 = gsl_ran_flat(r01, -1.0, 1.0);
        u1 = gsl_ran_flat(r01, -1.0, 1.0);

        sq_u01 = u0 * u0 + u1 * u1;

    } while (sq_u01 > 1);

    do {
        u2 = gsl_ran_flat(r01, -1.0, 1.0);
        u3 = gsl_ran_flat(r01, -1.0, 1.0);

        sq_u23 = u2 * u2 + u3 * u3;

    } while (sq_u23 > 1);

    square_term = sqrt((1 - sq_u01) / sq_u23);

    q_t.w = u0;
    q_t.x = u1;
    q_t.y = u2 * square_term;
    q_t.z = u3 * square_term;

    // scalar_q = q_t.x*Particles.N_Particle[id]->q.x + q_t.y *
    // Particles.N_Particle[id]->q.y + q_t.z*Particles.N_Particle[id]->q.z +
    // q_t.w * Particles.N_Particle[id]->q.w;

    //}while(scalar_q < 0.95);

    Particles.N_Particle[id]->q.w = q_t.w;
    Particles.N_Particle[id]->q.x = q_t.x;
    Particles.N_Particle[id]->q.y = q_t.y;
    Particles.N_Particle[id]->q.z = q_t.z;
}
