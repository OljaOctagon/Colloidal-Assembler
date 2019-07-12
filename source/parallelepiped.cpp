#include "parallelepiped.h"

parallelepiped::parallelepiped() {

    edge_N = 8;
    N_independent_faces = 3;
    N_cross_edges = 3;
    edge_N_vec = 3;

    x = new double[edge_N];
    y = new double[edge_N];
    z = new double[edge_N];

    dist_x = new double[edge_N];
    dist_y = new double[edge_N];
    dist_z = new double[edge_N];

    new_dist_x = new double[edge_N];
    new_dist_y = new double[edge_N];
    new_dist_z = new double[edge_N];

    edges = new m_vector[edge_N_vec];
    facenormal = new m_vector[N_independent_faces];

    edge_out = new int[edge_N];

    copy_count = 0;

    trans_periodic = new double[3];

    trans_old = new double[3];
    Rot_old = new double[9];
}

parallelepiped::~parallelepiped() {

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

    delete[] trans_periodic;

    delete[] trans_old;
    delete[] Rot_old;
}

void parallelepiped::write_positions(ofstream &fout) {}

void parallelepiped::edges_from_center() { // only valid if cubes are aligned
                                           // with coordinate axes!

    x[0] = x_center + (-v1.x - v2.x - v3.x) / 2.0;
    y[0] = y_center + (-v1.y - v2.y - v3.y) / 2.0;
    z[0] = z_center + (-v1.z - v2.z - v3.z) / 2.0;

    x[1] = x_center + (v1.x - v2.x - v3.x) / 2.0;
    y[1] = y_center + (v1.y - v2.y - v3.y) / 2.0;
    z[1] = z_center + (v1.z - v2.z - v3.z) / 2.0;

    x[2] = x_center + (v1.x + v2.x - v3.x) / 2.0;
    y[2] = y_center + (v1.y + v2.y - v3.y) / 2.0;
    z[2] = z_center + (v1.z + v2.z - v3.z) / 2.0;

    x[3] = x_center + (-v1.x + v2.x - v3.x) / 2.0;
    y[3] = y_center + (-v1.y + v2.y - v3.y) / 2.0;
    z[3] = z_center + (-v1.z + v2.z - v3.z) / 2.0;

    x[4] = x_center + (-v1.x - v2.x + v3.x) / 2.0;
    y[4] = y_center + (-v1.y - v2.y + v3.y) / 2.0;
    z[4] = z_center + (-v1.z - v2.z + v3.z) / 2.0;

    x[5] = x_center + (v1.x - v2.x + v3.x) / 2.0;
    y[5] = y_center + (v1.y - v2.y + v3.y) / 2.0;
    z[5] = z_center + (v1.z - v2.z + v3.z) / 2.0;

    x[6] = x_center + (v1.x + v2.x + v3.x) / 2.0;
    y[6] = y_center + (v1.y + v2.y + v3.y) / 2.0;
    z[6] = z_center + (v1.z + v2.z + v3.z) / 2.0;

    x[7] = x_center + (-v1.x + v2.x + v3.x) / 2.0;
    y[7] = y_center + (-v1.y + v2.y + v3.y) / 2.0;
    z[7] = z_center + (-v1.z + v2.z + v3.z) / 2.0;
}

void parallelepiped::distance_from_center() {

    for (int j = 0; j < edge_N; j++) {

        dist_x[j] = x[j] - x_center;
        dist_y[j] = y[j] - y_center;
        dist_z[j] = z[j] - z_center;
    }

    // cout<<"Distance from center"<<endl;
}

void parallelepiped::Calculate_Axis() {

    double norm_ax;

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

    ax_3.x = x[4] - x[0];
    ax_3.y = y[4] - y[0];
    ax_3.z = z[4] - z[0];

    norm_ax = ax_3.norm();

    ax_3.x = ax_3.x / norm_ax;
    ax_3.y = ax_3.y / norm_ax;
    ax_3.z = ax_3.z / norm_ax;
}

void parallelepiped::Set_Axis() {

    ax_1_old.x = ax_1.x;
    ax_1_old.y = ax_1.y;
    ax_1_old.z = ax_1.z;

    ax_2_old.x = ax_2.x;
    ax_2_old.y = ax_2.y;
    ax_2_old.z = ax_2.z;

    ax_3_old.x = ax_3.x;
    ax_3_old.y = ax_3.y;
    ax_3_old.z = ax_3.z;

    // cout<<"Set Axis"<<endl;
}

void parallelepiped::Set_Lengths() {

    /*
    Lx = 1.0;
    Ly = 1.0;
    Lz = 1.0;
    */
    // alpha=80.0*M_PI/180.0;
    // beta=20.0*M_PI/180;
    // gamma = 40.0*M_PI/180.0;

    // alpha_2 = M_PI - alpha;
    // beta_2  = M_PI - beta;
    // gamma_2 = M_PI - gamma;

    v1.x = 1.0;
    v1.y = 0.0;
    v1.z = 0.0;

    v2.x = 0.0;
    v2.y = 1.3;
    v2.z = 0.2;

    v3.x = 0.3;
    v3.y = 0.2;
    v3.z = 1.3;

    Lx = v1.norm();
    Ly = v2.norm();
    Lz = v3.norm();

    alpha = acos((v1.x * v3.x + v1.y * v3.y + v1.z * v3.z) / (Lx * Lz));
    beta = acos((v2.x * v3.x + v2.y * v3.y + v2.z * v3.z) / (Ly * Lz));
    gamma = acos((v1.x * v2.x + v1.y * v2.y + v1.z * v2.z) / (Lx * Ly));

    vc.x = v1.x + v2.x + v2.x;
    vc.y = v1.y + v2.y + v2.y;
    vc.z = v1.z + v2.z + v2.z;

    cut_off = vc.norm();
    cut_off_squared = cut_off * cut_off;

    cross_p.x = v2.y * v3.z - v2.z * v3.y;
    cross_p.y = v2.z * v3.x - v2.x * v3.z;
    cross_p.z = v2.x * v3.y - v2.y * v3.x;

    V = v1.x * cross_p.x + v1.y * cross_p.y + v1.z * cross_p.z;
}

void parallelepiped::Set_Start_Lattice_Position(int id, double box_Lx,
                                                int N_box) {

    // for cubic lattice
    int N_sitesp;
    double N_sitesp_float;
    double l_distp;

    N_sitesp = rint(pow(double(N_box), 1. / 3.));
    N_sitesp_float = double(N_sitesp);

    l_distp = (box_Lx - N_sitesp_float) / N_sitesp_float;

    x_center = 0.5 + l_distp / 2.0 +
               double(id % N_sitesp + l_distp * (id % N_sitesp));
    y_center = 0.5 + l_distp / 2.0 +
               double((id / N_sitesp) % N_sitesp +
                      l_distp * ((id / N_sitesp) % N_sitesp));
    z_center = 0.5 + l_distp / 2.0 +
               double(id / int(pow((double)N_sitesp, 2)) +
                      l_distp * (id / int((double)pow(N_sitesp, 2))));

    edges_from_center();

    /*
    int N_sitesp1, N_sitesp2, N_sitesp3;
    double N_sitesp1_float, N_sitesp2_float, N_sitesp3_float;
    double l_distpx, l_distpy, l_distpz;

    double a,b;
    double a_q, b_q, c_q;

    a_q=vc.x;
    b_q=vc.y;
    c_q=vc.z;


    a=b_q/a_q;
    b=c_q/a_q;


    N_sitesp1=floor(pow(double(N_box)/double(a*b), 1./3.));
    N_sitesp2=floor(N_sitesp1*a);
    N_sitesp3=ceil(double(N_box)/double(N_sitesp1*N_sitesp2));


    N_sitesp1_float = double(N_sitesp1);
    N_sitesp2_float = double(N_sitesp2);
    N_sitesp3_float = double(N_sitesp3);


    l_distpx = (box_Lx - N_sitesp1_float*a_q)/N_sitesp1_float;
    l_distpy = (box_Lx - N_sitesp2_float*b_q)/N_sitesp2_float;
    l_distpz = (box_Lx - N_sitesp3_float*c_q)/N_sitesp3_float;

    cout<<"l1 "<<l_distpx<<endl;
    cout<<"l2 "<<l_distpy<<endl;
    cout<<"l3 "<<l_distpz<<endl;

    x_center = vc.x/2.0 + l_distpx/2.0       + double(id%N_sitesp1)*(vc.x +
    l_distpx ); y_center = vc.y/2.0 + l_distpy/2.0     +
    double((id/N_sitesp1)%N_sitesp2)*(vc.y + l_distpy ); z_center = vc.z/2.0 +
    l_distpz/2.0 + double(id/int(N_sitesp1*N_sitesp2))*(vc.z + l_distpz );


    edges_from_center();
    */
}

void parallelepiped::Calculate_Face_Normals() {

    double face_ax;

    // face normals are the same as axis in the cube:

    facenormal[0].x = ax_1.y * ax_2.z - ax_1.z * ax_2.y;
    facenormal[0].y = ax_1.z * ax_2.x - ax_1.x * ax_2.z;
    facenormal[0].z = ax_1.x * ax_2.y - ax_1.y * ax_2.x;

    face_ax = facenormal[0].norm();

    facenormal[0].x = facenormal[0].x / face_ax;
    facenormal[0].y = facenormal[0].y / face_ax;
    facenormal[0].z = facenormal[0].z / face_ax;

    facenormal[1].x = ax_1.y * ax_3.z - ax_1.z * ax_3.y;
    facenormal[1].y = ax_1.z * ax_3.x - ax_1.x * ax_3.z;
    facenormal[1].z = ax_1.x * ax_3.y - ax_1.y * ax_3.x;

    face_ax = facenormal[1].norm();

    facenormal[1].x = facenormal[1].x / face_ax;
    facenormal[1].y = facenormal[1].y / face_ax;
    facenormal[1].z = facenormal[1].z / face_ax;

    facenormal[2].x = ax_2.y * ax_3.z - ax_2.z * ax_3.y;
    facenormal[2].y = ax_2.z * ax_3.x - ax_2.x * ax_3.z;
    facenormal[2].z = ax_2.x * ax_3.y - ax_2.y * ax_3.x;

    face_ax = facenormal[2].norm();

    facenormal[2].x = facenormal[2].x / face_ax;
    facenormal[2].y = facenormal[2].y / face_ax;
    facenormal[2].z = facenormal[2].z / face_ax;

    edges[0].x = ax_1.x;
    edges[0].y = ax_1.y;
    edges[0].z = ax_1.z;

    edges[1].x = ax_2.x;
    edges[1].y = ax_2.y;
    edges[1].z = ax_2.z;

    edges[2].x = ax_3.x;
    edges[2].y = ax_3.y;
    edges[2].z = ax_3.z;

    // cout<<"Calculate Face Normals"<<endl;
}

double
parallelepiped::Calculate_Projection_to_Separating_Axis(m_vector laxis) {

    double Rp;
    double rmax;
    double rmin;
    double scp_oc;

    double norm_ax;

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

    Rp = fabs((rmax - rmin) / 2.0);
    // Rp=rmax-rmin;

    return Rp;
}
