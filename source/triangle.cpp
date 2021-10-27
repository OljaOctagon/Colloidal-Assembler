#include "triangle.h"
#include <iostream>

triangle::triangle() {

    edge_N = 6;
    N_independent_faces = 4;
    N_cross_edges = 4;
    edge_N_vec = 4;

    N_patches = 6;
    N_patch_types = 3;

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

    level = pt.get<double>("Triangle.Energy_Level");
    T  = pt.get<double>("System.Temperature");

    delta_energy = pt.get<double>("Triangle.Energy_Difference");

    p1 = level;
    p2 = - level;
    p3 = 0.;
    p4 = level * delta_energy;

    // patch_energy(patchtype, patchtype)

    // patch_energy(0,0): same type attractive
    // patch_energy(0,1): different types: 0

    // A-A
    patch_energy[0][0] = p1;
    // A-B
    patch_energy[0][1] = p3;
    // A-C
    patch_energy[0][2] = p2;


    // B-A
    patch_energy[1][0] = p3;
    // B-B
    patch_energy[1][1] = p3;
    // B-C
    patch_energy[1][2] = p3;


    // Rhombi  C-A
    // ------------------------
    // C-A
    patch_energy[2][0] = p2;
    // C-B
    patch_energy[2][1] = p3;
    // C-C
    patch_energy[2][2] = p1;
    // -----------------------

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

triangle::~triangle() {

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

void triangle::edges_from_center(){

    sinus = Lx * sin(alpha) * 1. / 3.;
    L_2 = Lx / 2.;
    H_2 = Lz / 2.;

    x[0] = x_center - L_2;
    y[0] = y_center - sinus;
    z[0] = z_center - H_2;

    x[1] = x_center + L_2;
    y[1] = y_center - sinus;
    z[1] = z_center - H_2;

    x[2] = x_center;
    y[2] = y_center + 2. * sinus;
    z[2] = z_center - H_2;

    x[3] = x_center - L_2;
    y[3] = y_center - sinus;
    z[3] = z_center + H_2;

    x[4] = x_center + L_2;
    y[4] = y_center - sinus;
    z[4] = z_center + H_2;

    x[5] = x_center;
    y[5] = y_center + 2. * sinus;
    z[5] = z_center + H_2;

}

void triangle::distance_from_center() {

    for (int j = 0; j < edge_N; j++) {

            dist_x[j] = x[j] - x_center;
            dist_y[j] = y[j] - y_center;
            dist_z[j] = z[j] - z_center;
        }

}

void triangle::Calculate_Axis() {

    double norm_ax;

    ax_1.x = x[1] - x[0];
    ax_1.y = y[1] - y[0];
    ax_1.z = z[1] - z[0];

    norm_ax = ax_1.norm();

    ax_1.x = ax_1.x / norm_ax;
    ax_1.y = ax_1.y / norm_ax;
    ax_1.z = ax_1.z / norm_ax;

    ax_2.x = x[2] - x[0];
    ax_2.y = y[2] - y[0];
    ax_2.z = z[2] - z[0];

    norm_ax = ax_2.norm();

    ax_2.x = ax_2.x / norm_ax;
    ax_2.y = ax_2.y / norm_ax;
    ax_2.z = ax_2.z / norm_ax;

    ax_3.x = x[3] - x[0];
    ax_3.y = y[3] - y[0];
    ax_3.z = z[3] - z[0];


    norm_ax = ax_3.norm();

    ax_3.x = ax_3.x / norm_ax;
    ax_3.y = ax_3.y / norm_ax;
    ax_3.z = ax_3.z / norm_ax;

}

void triangle::Calculate_Long_Axis() {

    double norm_ax;

    long_axis.x = (x[1]+x[2]) - x[0];
    long_axis.y = (y[1]+y[2]) - y[0];
    long_axis.z = (z[1]+z[2]) - z[0];

    norm_ax = long_axis.norm();

    long_axis.x = long_axis.x / norm_ax;
    long_axis.y = long_axis.y / norm_ax;
    long_axis.z = long_axis.z / norm_ax;

}

void triangle::Set_Axis() {

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

void triangle::Set_Lengths() {

    // Triangle
    alpha = (60 * M_PI) / 180.0;

    Lx = 1.0;
    Ly = Lx;
    Lz = 0.1 * Lx;
    H_2 = Ly / 2.;

    h = Lx * sin(alpha);
    h_2 = double(h) / 2.0;
    a_x = Lx * cos(alpha);


    //diag2_short = sqrt((L*L + h * h);
    //diag2_long = sqrt((Lx + a_x) * (Lx + a_x) + h * h);

    //diag3_short = sqrt(diag2_short * diag2_short + Lz * Lz);
    //diag3_long = sqrt(diag2_long * diag2_long + Lz * Lz);

    //cut_off = vc.norm();

    A = Lx * Lx * sqrt(3.) / 4.;
    V = A * Lz;

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("para.ini", pt);

    patch_size = pt.get<double>("Triangle.patch_size");

    cut_off = (2.0 * Lx / sqrt(3.)) + 2. * patch_size + 0.5;
    cut_off_squared = cut_off * cut_off;

    r_patch[0] = patch_size;
    r_patch[1] = patch_size;
    r_patch[2] = patch_size;
    r_patch[3] = patch_size;
    r_patch[4] = patch_size;
    r_patch[5] = patch_size;

    patch_cutoff[0] = r_patch[0] * 2;
    patch_cutoff_squared[0] = patch_cutoff[0] * patch_cutoff[0];

    patch_cutoff[1] = r_patch[1] * 2;
    patch_cutoff_squared[1] = patch_cutoff[1] * patch_cutoff[1];

    patch_cutoff[2] = r_patch[2] * 2;
    patch_cutoff_squared[2] = patch_cutoff[2] * patch_cutoff[2];

    patch_cutoff[3] = r_patch[3] * 2;
    patch_cutoff_squared[3] = patch_cutoff[3] * patch_cutoff[3];

    patch_cutoff[4] = r_patch[4] * 2;
    patch_cutoff_squared[4] = patch_cutoff[4] * patch_cutoff[4];

    patch_cutoff[5] = r_patch[5] * 2;
    patch_cutoff_squared[5] = patch_cutoff[5] * patch_cutoff[5];

    patch_delta = pt.get<double>("Triangle.patch_delta");
    triangle_type = pt.get<string>("Triangle.triangle_type");

    //available types:
    //six_patch, three_asymm, two_neighbour_fixedcorner, two_opposite_fixedcorner

    patch_x = 0.5;
    d0 = patch_x;
    d1 = patch_x;
    d2 = patch_x;
    d3 = patch_x;
    d4 = patch_x;
    d5 = patch_x;

    patch_type[0] = 0;
    patch_type[1] = 0;
    patch_type[2] = 0;
    patch_type[3] = 0;
    patch_type[4] = 0;
    patch_type[5] = 0;

    if (triangle_type.compare("six_patch") == 0) {
        d0 = patch_delta;
        d1 = patch_delta;
        d2 = patch_delta;
        d3 = patch_delta;
        d4 = patch_delta;
        d5 = patch_delta;

        patch_type[0] = 0;
        patch_type[1] = 0;
        patch_type[2] = 0;
        patch_type[3] = 0;
        patch_type[4] = 0;
        patch_type[5] = 0;
    }

    if (triangle_type.compare("three_asymm") == 0) {
        d0 = patch_delta;
        d1 = patch_x;
        d2 = patch_delta;
        d3 = patch_x;
        d4 = patch_delta;
        d5 = patch_x;

        patch_type[0] = 0;
        patch_type[1] = 1;
        patch_type[2] = 0;
        patch_type[3] = 1;
        patch_type[4] = 0;
        patch_type[5] = 1;
    }

    if (triangle_type.compare("two_neighbour_fixedcorner") == 0) {
        d0 = 0;
        d1 = 1 - patch_delta;
        d2 = patch_x;
        d3 = patch_x;
        d4 = patch_x;
        d5 = patch_x;

        patch_type[0] = 0;
        patch_type[1] = 0;
        patch_type[2] = 1;
        patch_type[3] = 1;
        patch_type[4] = 1;
        patch_type[5] = 1;
    }

    if (triangle_type.compare("two_opposite_fixedcorner") == 0) {
        d0 = 0;
        d1 = patch_x;
        d2 = patch_delta;
        d3 = patch_x;
        d4 = patch_x;
        d5 = patch_x;

        patch_type[0] = 0;
        patch_type[1] = 1;
        patch_type[2] = 0;
        patch_type[3] = 1;
        patch_type[4] = 1;
        patch_type[5] = 1;
    }

}

void triangle::Set_Lengths(int e0, int e1, int e2, int e3, int e4, int e5) {

    patch_type[0] = e0;
    patch_type[1] = e1;
    patch_type[2] = e2;
    patch_type[3] = e3;
    patch_type[4] = e4;
    patch_type[5] = e5;

}

void triangle::Set_Start_Lattice_Position(int id, double box_Lx, 
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

    l_distp_x = (box_Lx - N_sitesp_x_float * L) / N_sitesp_x_float;
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

void triangle::Set_Start_Lattice_Position(int id, double box_Lx, 
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

void triangle::Calculate_Patch_Position() {

    x_patch[0] = x[0] + d0 * (x[1] - x[0]);
    y_patch[0] = y[0] + d0 * (y[1] - y[0]);
    z_patch[0] = z_center;
    x_patch[1] = x[1] + d1 * (x[0] - x[1]);
    y_patch[1] = y[1] + d1 * (y[0] - y[1]);
    z_patch[1] = z_center;

    x_patch[2] = x[1] + d2 * (x[2] - x[1]);
    y_patch[2] = y[1] + d2 * (y[2] - y[1]);
    z_patch[2] = z_center;
    x_patch[3] = x[2] + d3 * (x[1] - x[2]);
    y_patch[3] = y[2] + d3 * (y[1] - y[2]);
    z_patch[3] = z_center;

    y_patch[4] = y[2] + d4 * (y[0] - y[2]);
    x_patch[4] = x[2] + d4 * (x[0] - x[2]);
    z_patch[4] = z_center;
    y_patch[5] = y[0] + d5 * (y[2] - y[0]);
    x_patch[5] = x[0] + d5 * (x[2] - x[0]);
    z_patch[5] = z_center;

}

void triangle::Calculate_Face_Normals() {

    double face_ax;

    edges[0].x =  x[1] - x[0];
    edges[0].y =  y[1] - y[0];
    edges[0].z =  z[1] - z[0];

    edges[1].x =  x[0] - x[2];
    edges[1].y =  y[0] - y[2];
    edges[1].z =  z[0] - z[2];

    edges[2].x =  x[2] - x[1];
    edges[2].y =  y[2] - y[1];
    edges[2].z =  z[2] - z[1];

    edges[3].x =  x[3] - x[0];
    edges[3].y =  y[3] - y[0];
    edges[3].z =  z[3] - z[0];


    facenormal[0].x = edges[0].y * edges[1].z - edges[0].z * edges[1].y;
    facenormal[0].y = edges[0].z * edges[1].x - edges[0].x * edges[1].z;
    facenormal[0].z = edges[0].x * edges[1].y - edges[0].y * edges[1].x;

    facenormal[1].x = edges[0].y * edges[3].z - edges[0].z * edges[3].y;
    facenormal[1].y = edges[0].z * edges[3].x - edges[0].x * edges[3].z;
    facenormal[1].z = edges[0].x * edges[3].y - edges[0].y * edges[3].x;

    facenormal[2].x = edges[1].y * edges[3].z - edges[1].z * edges[3].y;
    facenormal[2].y = edges[1].z * edges[3].x - edges[1].x * edges[3].z;
    facenormal[2].z = edges[1].x * edges[3].y - edges[1].y * edges[3].x;

    facenormal[3].x = edges[2].y * edges[3].z - edges[2].z * edges[3].y;
    facenormal[3].y = edges[2].z * edges[3].x - edges[2].x * edges[3].z;
    facenormal[3].z = edges[2].x * edges[3].y - edges[2].y * edges[3].x;


    double f_norm;

    for (int k = 0; k < N_independent_faces; k++) {

        f_norm = facenormal[k].norm();

        facenormal[k].x = facenormal[k].x / f_norm;
        facenormal[k].y = facenormal[k].y / f_norm;
        facenormal[k].z = facenormal[k].z / f_norm;
    }


    /*
    facenormal[0].x = (x[0] + x[1] + x[4] + x[3])/4. - x_center;
    facenormal[0].y = (y[0] + y[1] + y[4] + y[3])/4. - y_center;
    facenormal[0].z = (z[0] + z[1] + z[4] + z[3])/4. - z_center;

    face_ax = facenormal[0].norm();

    facenormal[0].x = facenormal[0].x / face_ax;
    facenormal[0].y = facenormal[0].y / face_ax;
    facenormal[0].z = facenormal[0].z / face_ax;

    facenormal[1].x = (x[1] + x[2] + x[5] + x[4])/4. - x_center;
    facenormal[1].y = (y[1] + y[2] + y[5] + y[4])/4. - y_center;
    facenormal[1].z = (z[1] + z[2] + z[5] + z[4])/4. - z_center;

    face_ax = facenormal[1].norm();

    facenormal[1].x = facenormal[1].x / face_ax;
    facenormal[1].y = facenormal[1].y / face_ax;
    facenormal[1].z = facenormal[1].z / face_ax;

    facenormal[2].x = (x[2] + x[0] + x[3] + x[5])/4. - x_center;
    facenormal[2].y = (y[2] + y[0] + y[3] + y[5])/4. - y_center;
    facenormal[2].z = (z[2] + z[0] + z[3] + z[5])/4. - z_center;

    face_ax = facenormal[2].norm();

    facenormal[2].x = facenormal[2].x / face_ax;
    facenormal[2].y = facenormal[2].y / face_ax;
    facenormal[2].z = facenormal[2].z / face_ax;


    facenormal[3].x = (x[0] + x[1] + x[2])/3. - x_center;
    facenormal[3].y = (y[0] + y[1] + y[2])/3. - y_center;
    facenormal[3].z = (z[0] + z[1] + z[2])/3. - z_center;

    face_ax = facenormal[3].norm();

    facenormal[3].x = facenormal[3].x / face_ax;
    facenormal[3].y = facenormal[3].y / face_ax;
    facenormal[3].z = facenormal[3].z / face_ax;

    */

    

}

double triangle::Calculate_Projection_to_Separating_Axis(m_vector laxis) {
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
    //Rp=rmax-rmin;
    return Rp;

}