/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/**
 * @file work_lda.c
 * @brief This file is to be included in LDA functionals. As often these
 *        functionals are written as a function of rs and zeta, this
 *        routine performs the necessary conversions between this and a functional
 *        of rho.
 */

#define INPUT_COUNT 1
#define HIDDEN1_COUNT 1000
#define HIDDEN2_COUNT 500
#define OUTPUT_COUNT 1


#ifndef XC_DIMENSIONS
#define XC_DIMENSIONS 3
#endif

//float fc1_W[1000][2], fc2_W[500][1001], fcout_W[1][501];

float ** dot_product(float ** array1, float ** array2, int m, int n, int q)
{
	float ** multiply = malloc(m*sizeof(float *));
	int i, c, d, k;
	for(i=0; i< m; i++) multiply[i] = malloc(q*sizeof(float));
	for (c = 0; c < m; c++) {
	  for (d = 0; d < q; d++) {
		multiply[c][d] = 0;
		for (k = 0; k < n; k++) {
		  multiply[c][d] += array1[c][k]*array2[k][d];
		}
	  }
	}
	return multiply;
}

float ** relu(float ** array, int m, int n)
{
	float ** result = malloc(m*sizeof(float *));
	int i, j;
	for(i=0; i< m; i++) result[i] = malloc(n*sizeof(float));
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			if (array[i][j] > 0) {
				result[i][j] = array[i][j];
			}
			else {
				result[i][j] = 0;
			}
		}
	}
	return result;
}

float ** transpose(float ** array, int m, int n)
{
	float ** result = malloc(n*sizeof(float *));
	int i, c, d;
	for(i=0; i< n; i++) result[i] = malloc(m*sizeof(float));
	for (c = 0; c < m; c++) {
		for( d = 0 ; d < n ; d++ ) {
			result[d][c] = array[c][d];
		}
	}
	return result;
}

float ** preprocess_input(float ** array, int m)
{
	float ** result = malloc((m+1)*sizeof(float *));
	int i, c, d;
	for(i=0; i < (m+1); i++) result[i] = (float *)malloc(sizeof(float));
	for (c = 1; c < (m+1); c++) {
		for( d = 0 ; d < 1 ; d++ ) {
			result[c][d] = array[c-1][d];
		}
	}
	result[0][0] = 1;
	return result;	
}
				


float ** lda_function(float ** inp, float ** W1, float ** W2, float ** Wout)
{
	// W1, W2, Wout are already transposed for speed up calculations
	float ** input_preprocessed = preprocess_input(inp, 1);
	float ** part1 = dot_product(W1, input_preprocessed, 1000, 2, 1);
	float ** part1_act = relu(part1, 1000, 1);
	float ** part1_preprocessed = preprocess_input(part1_act, 1000);
	float ** part2 = dot_product(W2, part1_preprocessed, 500, 1001, 1);
	float ** part2_act = relu(part2, 500, 1);
	float ** part2_preprocessed = preprocess_input(part2_act, 500);
	float ** out_part = dot_product(Wout, part2_preprocessed, 1, 501, 1);
	return out_part;
}


float ** read_weights(char * filename, int m, int n)
{
	FILE* f;
	int ii, jj, i;
	float ** result = malloc(m*sizeof(float *));
	for(i=0; i< m; i++) result[i] = malloc(n*sizeof(float));
	f = fopen(filename, "r");
	for(jj=0; jj<m; jj++) {
		for(ii=0; ii<n; ii++) {
			fscanf(f, "%e", &result[jj][ii]);
		}
	}
	fclose(f);
	return result;
}
















/**
 * @param[in,out] func_type: pointer to pspdata structure to be initialized
 */
static void 
work_lda_nn(const xc_func_type *p, int np, const double *rho, 
	 double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  xc_lda_work_t r;
  int is, ip;
  double dens, drs, d2rs, d3rs;

  /* Wigner radius */
# if   XC_DIMENSIONS == 1
  const double cnst_rs = 0.5;
# elif XC_DIMENSIONS == 2
  const double cnst_rs = 1.0/M_SQRTPI;
# else /* three dimensions */
  const double cnst_rs = RS_FACTOR;
# endif


  /* Initialize memory */
  memset(&r, 0, sizeof(r));
  
  // load weights matrices
  
  float ** fc1_W = read_weights("fc1_W.txt", 1000, 2);
  float ** fc2_W = read_weights("fc2_W.txt", 500, 1001);
  float ** fcout_W = read_weights("fcout_W.txt", 1, 501);

  r.order = -1;
  if(zk     != NULL) r.order = 0;
  if(vrho   != NULL) r.order = 1;
  if(v2rho2 != NULL) r.order = 2;
  if(v3rho3 != NULL) r.order = 3;
  if(r.order < 0) return;
/*
  for(ip = 0; ip < np; ip++){
    xc_rho2dzeta(p->nspin, rho, &dens, &r.z);

    if(dens < p->dens_threshold) goto end_ip_loop;

    r.rs = cnst_rs*pow(dens, -1.0/XC_DIMENSIONS);

    func(p, &r);

    if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
      *zk = r.f;

    if(r.order < 1) goto end_ip_loop;

    drs = -r.rs/(XC_DIMENSIONS*dens);
    
    if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC)){
      vrho[0] = r.f + dens*r.dfdrs*drs;

      if(p->nspin == XC_POLARIZED){
	vrho[1] = vrho[0] - (r.z + 1.0)*r.dfdz;
	vrho[0] = vrho[0] - (r.z - 1.0)*r.dfdz;
      }
    }
*/

	for(ip = 0; ip < np; ip++){
		xc_rho2dzeta(p->nspin, rho, &dens, &r.z);
		
		if(dens < p->dens_threshold) goto end_ip_loop;
		
		r.rs = cnst_rs*pow(dens, -1.0/XC_DIMENSIONS);
		
		func(p, &r);
		
		if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
		*zk = r.f;
		
		if(r.order < 1) goto end_ip_loop;
		
		drs = -r.rs/(XC_DIMENSIONS*dens);
		
		if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC)){
			int i;
			float ** input = malloc(1*sizeof(float *));
			for(i=0; i< 1; i++) input[i] = malloc(1*sizeof(float));
			input[0][0] = (float)log10(rho[0]);
			float ** output = lda_function(input, fc1_W, fc2_W, fcout_W);
			vrho[0] = (double)output[0][0];
			
			if(p->nspin == XC_POLARIZED){
				vrho[1] = vrho[0] - (r.z + 1.0)*r.dfdz;
				vrho[0] = vrho[0] - (r.z - 1.0)*r.dfdz;
			}
		}
  
    if(r.order < 2) goto end_ip_loop;
    
    d2rs = -drs*(1.0 + XC_DIMENSIONS)/(XC_DIMENSIONS*dens);
    
    if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC)){
      v2rho2[0] = r.dfdrs*(2.0*drs + dens*d2rs) + dens*r.d2fdrs2*drs*drs;
      
      if(p->nspin == XC_POLARIZED){
	double sign[3][2] = {{-1.0, -1.0}, {-1.0, +1.0}, {+1.0, +1.0}};
	
	for(is=2; is>=0; is--){
	  v2rho2[is] = v2rho2[0] - r.d2fdrsz*(2.0*r.z + sign[is][0] + sign[is][1])*drs
	    + (r.z + sign[is][0])*(r.z + sign[is][1])*r.d2fdz2/dens;
	}
      }
    }
    
    if(r.order < 3) goto end_ip_loop;

    d3rs = -d2rs*(1.0 + 2.0*XC_DIMENSIONS)/(XC_DIMENSIONS*dens);
    
    if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC)){
      v3rho3[0] = r.dfdrs*(3.0*d2rs + dens*d3rs) + 
	3.0*r.d2fdrs2*drs*(drs + dens*d2rs) + r.d3fdrs3*dens*drs*drs*drs;
      
      if(p->nspin == XC_POLARIZED){
	double sign[4][3] = {{-1.0, -1.0, -1.0}, {-1.0, -1.0, +1.0}, {-1.0, +1.0, +1.0}, {+1.0, +1.0, +1.0}};
	
	for(is=3; is>=0; is--){
	  double ff;
	  
	  v3rho3[is]  = v3rho3[0] - (2.0*r.z  + sign[is][0] + sign[is][1])*(d2rs*r.d2fdrsz + drs*drs*r.d3fdrs2z);
	  v3rho3[is] += (r.z + sign[is][0])*(r.z + sign[is][1])*(-r.d2fdz2/dens + r.d3fdrsz2*drs)/dens;
	  
	  ff  = r.d2fdrsz*(2.0*drs + dens*d2rs) + dens*r.d3fdrs2z*drs*drs;
	  ff += -2.0*r.d2fdrsz*drs - r.d3fdrsz2*(2.0*r.z + sign[is][0] + sign[is][1])*drs;
	  ff += (r.z + sign[is][0])*(r.z + sign[is][1])*r.d3fdz3/dens;
	  ff += (2.0*r.z  + sign[is][0] + sign[is][1])*r.d2fdz2/dens;
	  
	  v3rho3[is] += -ff*(r.z + sign[is][2])/dens;
	}
      }
    }

  end_ip_loop:
    rho += p->n_rho;

    if(zk != NULL)
      zk += p->n_zk;
    
    if(vrho != NULL)
      vrho += p->n_vrho;

    if(v2rho2 != NULL)
      v2rho2 += p->n_v2rho2;

    if(v3rho3 != NULL)
      v3rho3 += p->n_v3rho3;

  } /* for(ip) */
}
