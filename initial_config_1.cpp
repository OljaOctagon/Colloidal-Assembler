#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

int nmol;
double** r;
double rho;
double vol;
double L[3];
double rcut;


static void write_configuration(void){
	int		i;
	FILE	*fh_xyz;
	
	fh_xyz = fopen("initial.xyz", "w");
	
	fprintf(fh_xyz, "%d\n\n", nmol);
	for (i=0; i<nmol; i++)
		fprintf(fh_xyz, "I %f %f %f\n", r[i][0], r[i][1], r[i][2]);
		
	fclose(fh_xyz);
}

static void build_fcc(void){
	int	i, j, k, z;
	int	cells;
	double	a;
	double	test;

	test = cbrt( (double)nmol/4.0 );
	cells = (int)( cbrt( (double)nmol/4.0 ) + 1.0e-14 );
	if (gsl_fcmp((double)cells, test, 1.0e-12)){
		fprintf(stderr, "Wrong number of molecules. %d instead of %d or %d\n", nmol, 4*cells*cells*cells, 4*(cells+1)*(cells+1)*(cells+1));
		exit(EXIT_FAILURE);
	}
	a = cbrt( 4.0/rho );

	vol = (double)nmol/rho;

	L[0] = cells*a;
	L[1] = cells*a;
	L[2] = cells*a;

	z = 0;
	for (i=0; i<cells; i++)
	for (j=0; j<cells; j++)
	for (k=0; k<cells; k++){
		r[z][0] = i*a		 - 0.5*L[0];
		r[z][1] = j*a		 - 0.5*L[1];
		r[z][2] = k*a		 - 0.5*L[2];
		z++;
		r[z][0] = i*a+0.5*a	 - 0.5*L[0];
		r[z][1] = j*a+0.5*a	 - 0.5*L[1];
		r[z][2] = k*a		 - 0.5*L[2];
		z++;
		r[z][0] = i*a		 - 0.5*L[0];
		r[z][1] = j*a+0.5*a	 - 0.5*L[1];
		r[z][2] = k*a+0.5*a	 - 0.5*L[2];
		z++;
		r[z][0] = i*a+0.5*a	 - 0.5*L[0];
		r[z][1] = j*a		 - 0.5*L[1];
		r[z][2] = k*a+0.5*a	 - 0.5*L[2];
		z++;
	}
}

static void build_hcp(void){
	int	i, j, k, z;
	int		cells[3];
	double	a;
	int		npa[23][23][23];
	bool	found;
	double	s[3];
	double	b1[3] = {1.0, 0.0,					0.0}; // first basis vector (with unit lattice constant)
	double	b2[3] = {0.5, 0.866025403784439,	0.0}; // second basis vector
	double	b3[3] = {0.0, 0.0,	   1.63299316185545}; // third basis vector
	double	a1[3] = {0.0, 0.0, 0.0};			// first atom vector
	double	a2[3] = {1.0/3.0, 1.0/3.0, 0.5}; 	// second atom vector
	
	/* possible number of cells and corresponding number of particles */
	for (i=0; i<23; i++)
	for (j=0; j<23; j++)
	for (k=0; k<23; k++)
		npa[i][j][k] = 0;
	
	npa[5][6][3] = 180;	
	npa[7][8][4] = 448;
	npa[8][10][5] = 800;
	npa[10][12][6] = 1440;
	npa[12][14][7] = 2352;
	npa[13][16][8] = 3328;
	npa[14][16][9] = 2016;
	npa[16][18][10] = 5760;
	npa[17][20][10] = 6800;
	npa[19][22][12] = 10032;						
	
	found = false;
	for (i=0; i<22; i++)
	for (j=0; j<22; j++)
	for (k=0; k<22; k++)
	if (npa[i][j][k] == nmol){
		cells[0] = i;
		cells[1] = j;
		cells[2] = k;
		found = true;
	}
	
	if (!found){
		fprintf(stderr, "Wrong number of molecules. %d instead of:\n	", nmol);
		for (i=0; i<23; i++)
		for (j=0; j<23; j++)
		for (k=0; k<23; k++)
		if (npa[i][j][k] > 0)
			fprintf(stderr, "%d, ", npa[i][j][k]);
		fprintf(stderr, "...\n");
		
		exit(EXIT_FAILURE);
	} 

	vol = (double)nmol/rho;
	a = cbrt( vol/((double)cells[0]*cells[1]*b2[1]*cells[2]*b3[2]) );

	L[0] = cells[0]*a;
	L[1] = cells[1]*b2[1]*a;
	L[2] = cells[2]*b3[2]*a;

	z = 0;
	for (i=0; i<cells[0]; i++)
	for (j=0; j<cells[1]; j++)
	for (k=0; k<cells[2]; k++){
		s[0] = i*a*b1[0] + j*a*b2[0] + k*a*b3[0];
		s[1] = i*a*b1[1] + j*a*b2[1] + k*a*b3[1];
		s[2] = i*a*b1[2] + j*a*b2[2] + k*a*b3[2];
		
		r[z][0] = s[0] + a1[0]*a*b1[0] + a1[1]*a*b2[0] + a1[2]*a*b3[0];
		r[z][1] = s[1] + a1[0]*a*b1[1] + a1[1]*a*b2[1] + a1[2]*a*b3[1];
		r[z][2] = s[2] + a1[0]*a*b1[2] + a1[1]*a*b2[2] + a1[2]*a*b3[2];
		z++;
		
		r[z][0] = s[0] + a2[0]*a*b1[0] + a2[1]*a*b2[0] + a2[2]*a*b3[0];
		r[z][1] = s[1] + a2[0]*a*b1[1] + a2[1]*a*b2[1] + a2[2]*a*b3[1];
		r[z][2] = s[2] + a2[0]*a*b1[2] + a2[1]*a*b2[2] + a2[2]*a*b3[2];
		z++;
	}
	
	for (i=0; i<nmol; i++)
	for (j=0; j<3; j++)
		r[i][j] = r[i][j] - rint(r[i][j]/L[j])*L[j];
}

static void build_bcc(void){
	int	i, j, k, z;
	int	cells;
	double	a;
	double	test;

	test = cbrt( (double)nmol/2.0 );
	cells = (int)( cbrt( (double)nmol/2.0 ) + 1.0e-14 );
	if (gsl_fcmp((double)cells, test, 1.0e-12)){
		fprintf(stderr, "Wrong number of molecules. %d instead of %d or %d\n", nmol, 2*cells*cells*cells, 2*(cells+1)*(cells+1)*(cells+1));
		exit(EXIT_FAILURE);
	}
	a = cbrt( 2.0/rho );
	fprintf(stderr, "a=%lf\n", a);

	vol = (double)nmol/rho;

	L[0] = cells*a;
	L[1] = cells*a;
	L[2] = cells*a;

	z = 0;
	for (i=0; i<cells; i++)
	for (j=0; j<cells; j++)
	for (k=0; k<cells; k++){
		r[z][0] = i*a;
		r[z][1] = j*a;
		r[z][2] = k*a;
		z++;
		r[z][0] = i*a+0.5*a;
		r[z][1] = j*a+0.5*a;
		r[z][2] = k*a+0.5*a;
		z++;
	}
	
	for (i=0; i<nmol; i++)
	for (j=0; j<3; j++)
		r[i][j] = r[i][j] - rint(r[i][j]/L[j])*L[j];
}

static void build_sc(void){
		int	i, j, k, z;
	int	cells;
	double	a;
	double	test;

	test = cbrt( (double)nmol );
	cells = (int)( cbrt( (double)nmol ) + 1.0e-14 );
	if (gsl_fcmp((double)cells, test, 1.0e-12)){
		fprintf(stderr, "Wrong number of molecules. %d instead of %d or %d\n", nmol, cells*cells*cells, (cells+1)*(cells+1)*(cells+1));
		exit(EXIT_FAILURE);
	}
	a = cbrt( 1.0/rho );

	vol = (double)nmol/rho;

	L[0] = cells*a;
	L[1] = cells*a;
	L[2] = cells*a;

	z = 0;
	for (i=0; i<cells; i++)
	for (j=0; j<cells; j++)
	for (k=0; k<cells; k++){
		r[z][0] = i*a;
		r[z][1] = j*a;
		r[z][2] = k*a;
		z++;
	}
	
	for (i=0; i<nmol; i++)
	for (j=0; j<3; j++)
		r[i][j] = r[i][j] - rint(r[i][j]/L[j])*L[j];
}

int main(int argc, char** argv){
	int lattice;
	
	if (argc != 4){
		fprintf(stderr, "USAGE:  %s #particle T rho lattice\n", argv[0]);
		fprintf(stderr, "\t  lattice:\t0 ... liquid \n");
		fprintf(stderr, "\t\t\t1 ... fcc \n");
		fprintf(stderr, "\t\t\t2 ... hcp \n");
		fprintf(stderr, "\t\t\t3 ... bcc \n");
		fprintf(stderr, "\t\t\t4 ... sc \n");
		fprintf(stderr, "\t\t\t5 ... fcc/liq slab \n");
		fprintf(stderr, "\t\t\t6 ... I-43d\n");
		exit(EXIT_FAILURE);
	}
	
	nmol = atoi(argv[1]);
	//Tset = atof(argv[2]);./
	rho = atof(argv[2]);
	lattice = atoi(argv[3]);
	
	r = new double*[nmol];
	for(int j=0;j<nmol;j++){
		 r[j] = new double[3];	
		}			
	// set default cutoff
	rcut = 2.5;
	
	
	if (lattice == 0)
		build_fcc();
	else if (lattice == 1)
		build_hcp();
	else if (lattice == 2)
		build_bcc();
	else if (lattice == 3)
		build_sc();
	else{
		fprintf(stderr, "error: lattice %d unknown.\n", lattice);
		exit(EXIT_FAILURE);
	}
			
	
	//checkpoint_set("initial.cp");
	write_configuration();
		
	
	
	return EXIT_SUCCESS;
}

