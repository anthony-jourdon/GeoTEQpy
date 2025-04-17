/*
====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  GeoTEQpy
  filename: faulttools.h

  This file is part of GeoTEQpy.

  GeoTEQpy is free software: you can redistribute it and/or modify it under the terms 
  of the GNU General Public License as published by the Free Software Foundation, either 
  version 3 of the License, or any later version.

  GeoTEQpy is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with GeoTEQpy. 
  If not, see <https://www.gnu.org/licenses/>.
====================================================================================================
*/

#ifndef __faulttools_h__
#define __faulttools_h__

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "ft_constants.h"

#define NQP 1

void sort_ascendant(double w[], int idx[]);
int zheevj3(double complex A[3][3], double complex Q[3][3], double w[3]);
void compute_eigV(long int n_points, double E[], double eig_vec[], double eig_val[]);
void compute_eigV_sorted(long int n_points, double E[], double eig_vec[], double eig_val[]);
void compute_covariance_matrix_3d(long int npoints, double data[], double cov_matrix[][3]);
void compute_covariance_eigenvectors(long int npoints, double point_data[], double eig_val[], double eig_vec[][3]);

void ft_compute_cov_eig_vectors(float radius, double point_coor[], long int npoints, double vectors[]);

#endif
