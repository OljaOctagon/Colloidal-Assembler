
#include "pmove.h"

void pmove::Rotate(particles &Particles, box *Box, fileio &Fileio, int id,
                   int mc_time) {

    Particles.N_Particle[id]->trans_periodic[0] = 0.0;
    Particles.N_Particle[id]->trans_periodic[1] = 0.0;
    Particles.N_Particle[id]->trans_periodic[2] = 0.0;

    // Rot_Update_Quarternions_RANDOM(Particles, id);
    Rot_Update_Quarternions_VON_MISES(Particles, id);

    Rot_mat[0] =
        1.0 -
        2.0 * Particles.N_Particle[id]->q.y * Particles.N_Particle[id]->q.y -
        2.0 * Particles.N_Particle[id]->q.z * Particles.N_Particle[id]->q.z;
    Rot_mat[1] =
        2.0 * Particles.N_Particle[id]->q.x * Particles.N_Particle[id]->q.y -
        2.0 * Particles.N_Particle[id]->q.w * Particles.N_Particle[id]->q.z;
    Rot_mat[2] =
        2.0 * Particles.N_Particle[id]->q.x * Particles.N_Particle[id]->q.z +
        2.0 * Particles.N_Particle[id]->q.w * Particles.N_Particle[id]->q.y;

    Rot_mat[3] =
        2.0 * Particles.N_Particle[id]->q.x * Particles.N_Particle[id]->q.y +
        2.0 * Particles.N_Particle[id]->q.w * Particles.N_Particle[id]->q.z;
    Rot_mat[4] =
        1.0 -
        2.0 * Particles.N_Particle[id]->q.x * Particles.N_Particle[id]->q.x -
        2.0 * Particles.N_Particle[id]->q.z * Particles.N_Particle[id]->q.z;
    Rot_mat[5] =
        2.0 * Particles.N_Particle[id]->q.y * Particles.N_Particle[id]->q.z -
        2.0 * Particles.N_Particle[id]->q.w * Particles.N_Particle[id]->q.x;

    Rot_mat[6] =
        2.0 * Particles.N_Particle[id]->q.x * Particles.N_Particle[id]->q.z -
        2.0 * Particles.N_Particle[id]->q.w * Particles.N_Particle[id]->q.y;
    Rot_mat[7] =
        2.0 * Particles.N_Particle[id]->q.y * Particles.N_Particle[id]->q.z +
        2.0 * Particles.N_Particle[id]->q.w * Particles.N_Particle[id]->q.x;
    Rot_mat[8] =
        1.0 -
        2.0 * Particles.N_Particle[id]->q.x * Particles.N_Particle[id]->q.x -
        2.0 * Particles.N_Particle[id]->q.y * Particles.N_Particle[id]->q.y;

    // update positions
    Rot_Update_Positions(Particles, id, Rot_mat);
    // Which edges are outside the box
    Particles.Check_Periodic_CM(id, Box);

    // is particle center outside the box?
    Update_Periodic_Positions(Particles, Box, id);

    // update axis
    Particles.N_Particle[id]->Calculate_Axis();

    // update cell list
    Particles.Update_Cell_List(id, Box);

    // Calculate Collision Partners
    Particles.Collision_List[id].Calculate(
        Box, id, Particles.Id_Cell_List, Particles.Cell_List, Particles.Cell,
        Particles.N_Particle, Particles.N_Particle[0]->cut_off,
        Particles.MAX_coll_p);
    exit_status = 0;
    col_count = 0;

    // Collision Test
    Collision_Test(Particles, Box, id, Particles.Collision_List);

    if (exit_status >= 1) {

        Reset_Positions(Particles, id);

        Particles.N_Particle[id]->q.x = Particles.N_Particle_old[id]->q.x;
        Particles.N_Particle[id]->q.y = Particles.N_Particle_old[id]->q.y;
        Particles.N_Particle[id]->q.z = Particles.N_Particle_old[id]->q.z;
        Particles.N_Particle[id]->q.w = Particles.N_Particle_old[id]->q.w;

        Particles.N_Particle[id]->Calculate_Axis();

        Particles.Reset_Cell_List(Box, id, Particles.id_num);
    }

    if (exit_status == 0) {

        Particles.N_Particle_old[id]->q.x = Particles.N_Particle[id]->q.x;
        Particles.N_Particle_old[id]->q.y = Particles.N_Particle[id]->q.y;
        Particles.N_Particle_old[id]->q.z = Particles.N_Particle[id]->q.z;
        Particles.N_Particle_old[id]->q.w = Particles.N_Particle[id]->q.w;

        Set_Positions(Particles, id);

        accept_rotate = accept_rotate + 1;
    }

    Particles.N_Particle[id]->trans_periodic[0] = 0.0;
    Particles.N_Particle[id]->trans_periodic[1] = 0.0;
    Particles.N_Particle[id]->trans_periodic[2] = 0.0;

    N_rotate_moves = N_rotate_moves + 1;
}
