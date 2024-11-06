
#ifndef __faulttools_h__
#define __faulttools_h__

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "ft_constants.h"

#define NQP 1

int zheevj3(double complex A[3][3], double complex Q[3][3], double w[3]);
void compute_eigV(long int n_points, double E[], double eig_vec[], double eig_val[]);
void compute_eigV_sorted(long int n_points, double E[], double eig_vec[], double eig_val[]);
void compute_covariance_matrix_3d(long int npoints, double data[], double cov_matrix[][3]);

void ft_compute_cov_eig_vectors(float radius, double point_coor[], long int npoints, double vectors[]);

#endif
