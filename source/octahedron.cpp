

#include "octahedron.h"

octahedron::octahedron() {

    edge_N = 6;
    N_independent_faces = 4;
    N_cross_edges = 6;
    edge_N_vec = 10;

    x = new double[edge_N];
    y = new double[edge_N];
    z = new double[edge_N];

    dist_x = new double[edge_N];
    dist_y = new double[edge_N];
    dist_z = new double[edge_N];

    new_dist_x = new double[edge_N];
    new_dist_y = new double[edge_N];
    new_dist_z = new double[edge_N];

    facenormal = new m_vector[N_independent_faces];
    edges = new m_vector[edge_N_vec];

    // trans_init = new double[3];

    edge_out = new int[edge_N];

    trans_periodic = new double[3];
    trans_old = new double[3];
    Rot_old = new double[9];

    copy_count = 0;
}

octahedron::~octahedron() {

    delete[] x;
    delete[] y;
    delete[] z;

    delete[] dist_x;
    delete[] dist_y;
    delete[] dist_z;

    delete[] edge_out;

    delete[] new_dist_x;
    delete[] new_dist_y;
    delete[] new_dist_z;

    delete[] trans_old;
    delete[] Rot_old;
}

void octahedron::write_positions(ofstream &fout) {}

void octahedron::edges_from_center() { // for octahedron aligned with
                                       // coordinate ax

    x[0] = x_center - 0.5 * Lx;
    y[0] = y_center - 0.5 * Lx;
    z[0] = z_center;

    x[1] = x_center + 0.5 * Lx;
    y[1] = y_center - 0.5 * Lx;
    z[1] = z_center;

    x[2] = x_center + 0.5 * Lx;
    y[2] = y_center + 0.5 * Lx;
    z[2] = z_center;

    x[3] = x_center - 0.5 * Lx;
    y[3] = y_center + 0.5 * Lx;
    z[3] = z_center;

    x[4] = x_center;
    y[4] = y_center;
    z[4] = z_center + height;

    x[5] = x_center;
    y[5] = y_center;
    z[5] = z_center - height;
}

void octahedron::distance_from_center() {

    for (int j = 0; j < edge_N; j++) {

        dist_x[j] = x[j] - x_center;
        dist_y[j] = y[j] - y_center;
        dist_z[j] = z[j] - z_center;
    }
}

void octahedron::Calculate_Axis() {

    int norm_ax;

    ax_1.x = x[1] - x[0];
    ax_1.y = y[1] - y[0];
    ax_1.z = z[1] - z[0];

    norm_ax = ax_1.norm();

    ax_1.x = ax_1.x / norm_ax;
    ax_1.y = ax_1.y / norm_ax;
    ax_1.z = ax_1.z / norm_ax;

    ax_2.x = x[3] - x[0];
    ax_2.y = y[3] - y[0];
    ax_2.z = z[3] - z[0];

    norm_ax = ax_2.norm();

    ax_2.x = ax_2.x / norm_ax;
    ax_2.y = ax_2.y / norm_ax;
    ax_2.z = ax_2.z / norm_ax;

    ax_3.x = x[4] - x_center;
    ax_3.y = y[4] - y_center;
    ax_3.z = z[4] - z_center;

    norm_ax = ax_3.norm();

    ax_3.x = ax_3.x / norm_ax;
    ax_3.y = ax_3.y / norm_ax;
    ax_3.z = ax_3.z / norm_ax;
}

void octahedron::Set_Axis() {

    ax_1_old.x = ax_1.x;
    ax_1_old.y = ax_1.y;
    ax_1_old.z = ax_1.z;

    ax_2_old.x = ax_2.x;
    ax_2_old.y = ax_2.y;
    ax_2_old.z = ax_2.z;

    ax_3_old.x = ax_3.x;
    ax_3_old.y = ax_3.y;
    ax_3_old.z = ax_3.z;
}

void octahedron::Set_Lengths() {

    // Lx = 1.0;
    Lx = sqrt(2.0);
    height = 1.0;

    cut_off = 2.0 * height;

    cut_off_squared = cut_off * cut_off;
    // V = (sqrt(2.0)/3.0)*Lx*Lx*Lx;
    V = 4.0 / 3.0;
}

void octahedron::Set_Start_Lattice_Position(int id, double box_Lx, int N_box) {

    // for cubic lattice

    int N_sitesp1, N_sitesp2, N_sitesp3;
    double N_sitesp1_float, N_sitesp2_float, N_sitesp3_float;
    double l_distp, l_distpz;

    N_sitesp1 = floor(box_Lx / Lx);
    N_sitesp2 = N_sitesp1;
    N_sitesp3 = floor(box_Lx / (2.0 * height));

    N_sitesp1_float = double(N_sitesp1);
    N_sitesp3_float = double(N_sitesp3);

    l_distp = (box_Lx - N_sitesp1_float * Lx) / N_sitesp1_float;
    l_distpz = (box_Lx - N_sitesp3_float * 2.0 * height) / N_sitesp3_float;

    x_center = (Lx + l_distp) / 2.0 + double(id % N_sitesp1) * (Lx + l_distp);
    y_center = (Lx + l_distp) / 2.0 +
               double((id / N_sitesp1) % N_sitesp2) * (Lx + l_distp);
    z_center =
        (height + l_distpz / 2.0) +
        double(id / int(N_sitesp1 * N_sitesp2)) * (2.0 * height + l_distpz);

    edges_from_center();
}

void octahedron::Calculate_Face_Normals() {

    edges[0].x = x[1] - x[0];
    edges[0].y = y[1] - y[0];
    edges[0].z = z[1] - z[0];

    edges[1].x = x[3] - x[0];
    edges[1].y = y[3] - y[0];
    edges[1].z = z[3] - z[0];

    edges[2].x = x[4] - x[0];
    edges[2].y = y[4] - y[0];
    edges[2].z = z[4] - z[0];

    edges[3].x = x[4] - x[1];
    edges[3].y = y[4] - y[1];
    edges[3].z = z[4] - z[1];

    edges[4].x = x[4] - x[2];
    edges[4].y = y[4] - y[2];
    edges[4].z = z[4] - z[2];

    edges[5].x = x[4] - x[3];
    edges[5].y = y[4] - y[3];
    edges[5].z = z[4] - z[3];

    edges[6].x = x[1] - x[2];
    edges[6].y = y[1] - y[2];
    edges[6].z = z[1] - z[2];

    edges[7].x = x[3] - x[2];
    edges[7].y = y[3] - y[2];
    edges[7].z = z[3] - z[2];

    edges[8].x = x[5] - x[0];
    edges[8].y = y[5] - y[0];
    edges[8].z = z[5] - z[0];

    edges[9].x = x[5] - x[2];
    edges[9].y = y[5] - y[2];
    edges[9].z = z[5] - z[2];

    facenormal[0].x = edges[0].y * edges[2].z - edges[0].z * edges[2].y;
    facenormal[0].y = edges[0].z * edges[2].x - edges[0].x * edges[2].z;
    facenormal[0].z = edges[0].x * edges[2].y - edges[0].y * edges[2].x;

    facenormal[1].x = edges[1].y * edges[2].z - edges[1].z * edges[2].y;
    facenormal[1].y = edges[1].z * edges[2].x - edges[1].x * edges[2].z;
    facenormal[1].z = edges[1].x * edges[2].y - edges[1].y * edges[2].x;

    facenormal[2].x = edges[6].y * edges[4].z - edges[6].z * edges[4].y;
    facenormal[2].y = edges[6].z * edges[4].x - edges[6].x * edges[4].z;
    facenormal[2].z = edges[6].x * edges[4].y - edges[6].y * edges[4].x;

    facenormal[3].x = edges[7].y * edges[4].z - edges[7].z * edges[4].y;
    facenormal[3].y = edges[7].z * edges[4].x - edges[7].x * edges[4].z;
    facenormal[3].z = edges[7].x * edges[4].y - edges[7].y * edges[4].x;

    /*

    facenormal[4].x = edges[0].y*edges[7].z - edges[0].z*edges[7].y;
    facenormal[4].y = edges[0].z*edges[7].x - edges[0].x*edges[7].z;
    facenormal[4].z = edges[0].x*edges[7].y - edges[0].y*edges[7].x;

    facenormal[5].x = edges[1].y*edges[7].z - edges[1].z*edges[7].y;
    facenormal[5].y = edges[1].z*edges[7].x - edges[1].x*edges[7].z;
    facenormal[5].z = edges[1].x*edges[7].y - edges[1].y*edges[7].x;

    facenormal[6].x = edges[6].y*edges[8].z - edges[6].z*edges[8].y;
    facenormal[6].y = edges[6].z*edges[8].x - edges[6].x*edges[8].z;
    facenormal[6].z = edges[6].x*edges[8].y - edges[6].y*edges[8].x;

    facenormal[7].x = edges[4].y*edges[8].z - edges[4].z*edges[8].y;
    facenormal[7].y = edges[4].z*edges[8].x - edges[4].x*edges[8].z;
    facenormal[7].z = edges[4].x*edges[8].y - edges[4].y*edges[8].x;
    */

    double f_norm;

    for (int k = 0; k < N_independent_faces; k++) {

        f_norm = facenormal[k].norm();

        facenormal[k].x = facenormal[k].x / f_norm;
        facenormal[k].y = facenormal[k].y / f_norm;
        facenormal[k].z = facenormal[k].z / f_norm;
    }
}

double octahedron::Calculate_Projection_to_Separating_Axis(m_vector laxis) {

    double scp_oc;

    double norm_ax;
    /*
    ax_1.x = edges[0].x;
    ax_1.y = edges[0].y;
    ax_1.z = edges[0].z;

    norm_ax = ax_1.norm();

    ax_1.x = ax_1.x/norm_ax;
    ax_1.y = ax_1.y/norm_ax;
    ax_1.z = ax_1.z/norm_ax;

    ax_2.x = edges[1].x;
    ax_2.y = edges[1].y;
    ax_2.z = edges[1].z;

    norm_ax = ax_2.norm();

    ax_2.x = ax_2.x/norm_ax;
    ax_2.y = ax_2.y/norm_ax;
    ax_2.z = ax_2.z/norm_ax;


    ax_3.x = x[4] - x_center;
    ax_3.y = y[4] - y_center;
    ax_3.z = z[4] - z_center;

    norm_ax = ax_3.norm();

    ax_3.x = ax_3.x/norm_ax;
    ax_3.y = ax_3.y/norm_ax;
    ax_3.z = ax_3.z/norm_ax;
    */

    distance_from_center();

    rmin = dist_x[0] * laxis.x + dist_y[0] * laxis.y + dist_z[0] * laxis.z;
    rmax = rmin;

    for (int j = 1; j < edge_N; j++) {

        scp_oc =
            dist_x[j] * laxis.x + dist_y[j] * laxis.y + dist_z[j] * laxis.z;

        if (scp_oc < rmin) {
            rmin = scp_oc;
        } else if (scp_oc > rmax) {
            rmax = scp_oc;
        }
    }

    Rp = rmax - rmin;

    // Rp = (fabs(ax_1.x*laxis.x + ax_1.y*laxis.y + ax_1.z*laxis.z) +
    // fabs(ax_2.x*laxis.x + ax_2.y*laxis.y + ax_2.z*laxis.z))*(Lx/2.0)
    //					+ fabs(ax_3.x*laxis.x + ax_3.y*laxis.y +
    //ax_3.z*laxis.z)*(height/2.0);

    return Rp;
}
