
#include "rectangle.h"
#include <iostream>

rectangle::rectangle() {

    edge_N = 8;
    N_independent_faces = 3;
    N_cross_edges = 3;
    edge_N_vec = 3;

    N_patches = 4;
    N_patch_types = 4;

    x_patch = new double[N_patches];
    y_patch = new double[N_patches];
    z_patch = new double[N_patches];

    r_patch = new double[N_patches];
    patch_cutoff = new double[N_patches];
    patch_cutoff_squared = new double[N_patches];

    patch_energy = new double *[N_patch_types];
    for (int j = 0; j < N_patch_types; j++) {
        patch_energy[j] = new double[N_patch_types];
    }

    patch_type = new int[N_patches];

    double p1, p2, p3, p4;
    double T;
    double patch_size;
    double delta_energy;
    double level;

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("para.ini", pt);

    level = pt.get<double>("Rhombus.Energy_Level");
    T  = pt.get<double>("System.Temperature");

    delta_energy = pt.get<double>("Rhombus.Energy_Difference");

    double T_start = 0.1;

    p1 = level * T_start;
    p2 = - level * T_start;
    p3 = 0.;
    p4 = level * T_start * delta_energy;

    // patch_energy(patchtype, patchtype)

    // patches types (see simulation_roadmap.png for reference)
    // patch_type 0: neutral
    // patch_type 1: yellow
    // patch_type 2: green
    // patch_type 3: red


    // Note : stage 3 loop formers have two types of particles
    // type one is (green, yellow,neutral,neutral) = (2,1,0,0)
    // type two is (yellow, red, neutral, neutral) = (1,3,0,0)
    // type three is (neutral, yellow, yellow, neutral) = (0,1,1,0)

    // O: Neutral: alsways zeros (p3)
    patch_energy[0][0] = p3;
    // 0-Yellow: 
    patch_energy[0][1] = p3;
    // 0-Green
    patch_energy[0][2] = p3;
    // 0-Green
    patch_energy[0][3] = p3;

    // 1: Yellow:
    // Yellow-Yellow: attracitve but less  (p4)
    patch_energy[1][1] = p4;
    // Yellow-Green: repulsive (p2)
    patch_energy[1][2] = p2;
    // Yellow-Red: repulsive (p2)
    patch_energy[1][3] = p2;

    // 2: Green:
    // Green-Green: repulsive (p2)
    patch_energy[2][2] = p2;
    // Green-Red: attractive (p1)
    patch_energy[2][3] = p1;


    // 3: Red :
    // Red-Red: repulsive (p2)
    patch_energy[3][3] = p2;


    // symmetry:
    for ( int t1=0; t1<N_patch_types; t1++){
      for (int t2=0;t2<N_patch_types; t2++ ){
        if (t1>t2){
          patch_energy[t1][t2] = patch_energy[t2][t1];
        }

      }
    }

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

rectangle::~rectangle() {

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

    delete[] x_patch;
    delete[] y_patch;
    delete[] z_patch;

    delete[] r_patch;
    delete[] patch_cutoff;
    delete[] patch_cutoff_squared;
}

void rectangle::edges_from_center() { // only valid if cubes are aligned
                                         // with coordinate axes!

    x[0] = x_center + (-ra.x - rb.x - rc.x) / 2.;
    y[0] = y_center + (-ra.y - rb.y - rc.y) / 2.;
    z[0] = z_center + (-ra.z - rb.z - rc.z) / 2.;

    x[1] = x_center + (+ra.x - rb.x - rc.x) / 2.;
    y[1] = y_center + (+ra.y - rb.y - rc.y) / 2.;
    z[1] = z_center + (+ra.z - rb.z - rc.z) / 2.;

    x[2] = x_center + (+ra.x + rb.x - rc.x) / 2.;
    y[2] = y_center + (+ra.y + rb.y - rc.y) / 2.;
    z[2] = z_center + (+ra.z + rb.z - rc.z) / 2.;

    x[3] = x_center + (-ra.x + rb.x - rc.x) / 2.;
    y[3] = y_center + (-ra.y + rb.y - rc.y) / 2.;
    z[3] = z_center + (-ra.z + rb.z - rc.z) / 2.;

    x[4] = x_center + (-ra.x - rb.x + rc.x) / 2.;
    y[4] = y_center + (-ra.y - rb.y + rc.y) / 2.;
    z[4] = z_center + (-ra.z - rb.z + rc.z) / 2.;

    x[5] = x_center + (+ra.x - rb.x + rc.x) / 2.;
    y[5] = y_center + (+ra.y - rb.y + rc.y) / 2.;
    z[5] = z_center + (+ra.z - rb.z + rc.z) / 2.;

    x[6] = x_center + (+ra.x + rb.x + rc.x) / 2.;
    y[6] = y_center + (+ra.y + rb.y + rc.y) / 2.;
    z[6] = z_center + (+ra.z + rb.z + rc.z) / 2.;

    x[7] = x_center + (-ra.x + rb.x + rc.x) / 2.;
    y[7] = y_center + (-ra.y + rb.y + rc.y) / 2.;
    z[7] = z_center + (-ra.z + rb.z + rc.z) / 2.;
}

void rectangle::distance_from_center() {

    for (int j = 0; j < edge_N; j++) {

        dist_x[j] = x[j] - x_center;
        dist_y[j] = y[j] - y_center;
        dist_z[j] = z[j] - z_center;
    }
}

void rectangle::Calculate_Axis() {

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

void rectangle::Calculate_Long_Axis() {

    double norm_ax;

    long_axis.x = x[2] - x[0];
    long_axis.y = y[2] - y[0];
    long_axis.z = z[2] - z[0];

    norm_ax = long_axis.norm();

    long_axis.x = long_axis.x / norm_ax;
    long_axis.y = long_axis.y / norm_ax;
    long_axis.z = long_axis.z / norm_ax;
}

void rectangle::Set_Axis() {

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

void rectangle::Set_Lengths() {

    alpha = (90 * M_PI) / 180.0;
    beta = M_PI - alpha;
    
    Lx = 1.0;
    Ly = 2.0;
    Lz = 0.1 * Lx;

    ra.x = Lx;
    ra.y = 0.0;
    ra.z = 0.0;

    rb.x = 0;
    rb.y = Ly;
    rb.z = 0;

    rc.x = 0.0;
    rc.y = 0.0;
    rc.z = Lz;

    diag2_short = sqrt((Lx*Lx + Ly*Ly));
    diag2_long = diag2_short;

    diag3_short = sqrt(diag2_short*diag2_short + Lz*Lz);
    diag3_long = diag3_short;
   

    // cut_off = diag3_long;

    vc.x = ra.x + rb.x + rc.x;
    vc.y = ra.y + rb.y + rc.y;
    vc.z = ra.z + rb.z + rc.z;

    cut_off = vc.norm();

    cut_off_squared = cut_off * cut_off;

    cross_p.x = rb.y * rc.z - rb.z * rc.y;
    cross_p.y = rb.z * rc.x - rb.x * rc.z;
    cross_p.z = rb.x * rc.y - rb.y * rc.x;

    V = fabs(ra.x * cross_p.x + ra.y * cross_p.y + ra.z * cross_p.z);
    A = Lx*Ly;


    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("para.ini", pt);

    patch_size = pt.get<double>("Rhombus.patch_size");

    r_patch[0] = patch_size;
    r_patch[1] = patch_size;
    r_patch[2] = patch_size;
    r_patch[3] = patch_size;

    patch_cutoff[0] = r_patch[0] * 2;
    patch_cutoff_squared[0] = patch_cutoff[0] * patch_cutoff[0];

    patch_cutoff[1] = r_patch[1] * 2;
    patch_cutoff_squared[1] = patch_cutoff[1] * patch_cutoff[1];

    patch_cutoff[2] = r_patch[2] * 2;
    patch_cutoff_squared[2] = patch_cutoff[2] * patch_cutoff[2];

    patch_cutoff[3] = r_patch[3] * 2;
    patch_cutoff_squared[3] = patch_cutoff[3] * patch_cutoff[3];

    patch_delta = pt.get<double>("Rhombus.patch_delta");
    rhombus_type = pt.get<string>("Rhombus.rhombus_type");

    // available types:
    // 2_patch_2_type_enclosing, 



    patch_x = 0.5;
    d0 = patch_x;
    d1 = patch_x;
    d2 = patch_x;
    d3 = patch_x;

    patch_type[0] = 0;
    patch_type[1] = 0;
    patch_type[2] = 0;
    patch_type[3] = 0;

    if (rhombus_type.compare("2_patch_2_type_enclosing")==0) {
        d0 = patch_delta;
        d1 = patch_x;
        d2 = patch_x;
        d3 = patch_x;

        // type one is (green, yellow,neutral,neutral) = (2,1,0,0)
        // type two is (yellow, red, neutral, neutral) = (1,3,0,0)
        // type three is (neutral, yellow, yellow, neutral) = (0,1,1,0)

        patch_type[0] = 2;
        patch_type[1] = 1;
        patch_type[2] = 0;
        patch_type[3] = 0;
    }


    halo_energy = 0;
    halo_cutoff = (2 * Lx) / 10.0;

    cut_off_squared = cut_off * cut_off + Lx;
}

void rectangle::Set_Lengths(int e0, int e1, int e2, int e3) {

    patch_type[0] = e0;
    patch_type[1] = e1;
    patch_type[2] = e2;
    patch_type[3] = e3;
}

void rectangle::Set_Start_Lattice_Position(int id, double box_Lx,
                                              int N_box) {
    // for cubic lattice
    int N_sitesp_x;
    int N_sitesp_yz;

    double N_sitesp_x_float;
    double N_sitesp_yz_float;

    double l_distp_x;
    double l_distp_yz;

    N_sitesp_x = rint(pow(double(N_box * sin(alpha) * sin(alpha)), 1. / 3.));
    N_sitesp_x_float = double(N_sitesp_x);

    N_sitesp_yz = ceil(N_sitesp_x_float / sin(alpha));
    N_sitesp_yz_float = double(N_sitesp_yz);

    l_distp_x = (box_Lx - N_sitesp_x_float * Lx) / N_sitesp_x_float;
    l_distp_yz = (box_Lx - N_sitesp_yz_float * h) / N_sitesp_yz_float;

    x_center = a_x / 2.0 + l_distp_x / 2.0 + double(id % N_sitesp_x) * Lx +
               double(id % N_sitesp_x) * l_distp_x;
    y_center = h_2 + l_distp_yz / 2.0 +
               double((id / N_sitesp_x) % N_sitesp_yz) * h +
               double((id / N_sitesp_x) % N_sitesp_yz) * l_distp_yz;
    z_center = h_2 + l_distp_yz / 2.0 +
               double(id / int(N_sitesp_x * N_sitesp_yz)) * h +
               double((id / N_sitesp_x) % N_sitesp_yz) * l_distp_yz;

    edges_from_center();
}

void rectangle::Set_Start_Lattice_Position(int id, double box_Lx,
                                              double box_Ly, double box_Lz,
                                              int N_box) {

    // for cubic lattice in 2D or anisotropic box shape

    int N_sitespx, N_sitespy, N_sitespz;
    double N_sitespx_float, N_sitespy_float, N_sitespz_float;
    double l_distpx, l_distpy, l_distpz;

    N_sitespx = rint(sqrt(double(N_box / sin(alpha))));
    N_sitespx_float = double(N_sitespx);
    N_sitespy = rint(N_sitespx * sin(alpha));
    N_sitespy_float = double(N_sitespy);

    l_distpx = (box_Lx - N_sitespx_float) / N_sitespx_float;
    l_distpy = (box_Ly - N_sitespy_float) / N_sitespy_float;

    l_distpx = l_distpx*0.7;
    l_distpy = l_distpy*0.7;

    x_center = Lx / 2.0 + l_distpx / 2.0 +
               double(id % N_sitespx + l_distpx * (id % N_sitespx));
    y_center = h / 2.0 + l_distpy / 2.0 +
               double((id / N_sitespx) % N_sitespy +
                      l_distpx * ((id / N_sitespy) % N_sitespx));
    z_center = Lz / 2.0;

    edges_from_center();
}

void rectangle::Calculate_Patch_Position() {

    x_patch[0] = x[0] + d0 * (x[7] - x[0]);
    y_patch[0] = y[0] + d0 * (y[7] - y[0]);
    z_patch[0] = z[0] + patch_x * (z[7] - z[0]);

    x_patch[1] = x[2] + d1 * (x[7] - x[2]);
    y_patch[1] = y[2] + d1 * (y[7] - y[2]);
    z_patch[1] = z[2] + patch_x * (z[7] - z[2]);

    x_patch[2] = x[0] + d2 * (x[5] - x[0]);
    y_patch[2] = y[0] + d2 * (y[5] - y[0]);
    z_patch[2] = z[0] + patch_x * (z[5] - z[0]);

    x_patch[3] = x[6] + d3 * (x[1] - x[6]);
    y_patch[3] = y[6] + d3 * (y[1] - y[6]);
    z_patch[3] = z[6] + patch_x * (z[1] - z[6]);
}

void rectangle::Calculate_Face_Normals() {

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
}

double rectangle::Calculate_Projection_to_Separating_Axis(m_vector laxis) {
    double Rp;
    double rmax;
    double rmin;
    double scp_oc;
    double norm_ax;

    distance_from_center();

    rmin =  dist_x[0]*laxis.x + dist_y[0]*laxis.y + dist_z[0]*laxis.z;
    rmax = rmin;

    for (int j=1;j<edge_N;j++){
            scp_oc = dist_x[j]*laxis.x + dist_y[j]*laxis.y + dist_z[j]*laxis.z;
            if (scp_oc < rmin) {
                    rmin = scp_oc;
            }
            else if (scp_oc > rmax) {
                    rmax = scp_oc;
            }
    }

    Rp = fabs((rmax-rmin)/2.0);
    return Rp;

}
