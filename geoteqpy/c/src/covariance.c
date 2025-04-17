/*
====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  GeoTEQpy
  filename: covariance.c

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

#include <stdio.h>
#include "faulttools.h"

static void get_mean(long int npoints, int ndim, double data[], double x_bar[])
{
  long int n;
  int      d;

  /* Initialize x_bar to 0.0 */
  memset(x_bar, 0.0, ndim*sizeof(double));
  /* Sum for each dimension (or type of variable) */
  for (n=0; n<npoints; n++) {
    for (d=0; d<ndim; d++) {
      x_bar[d] += data[ndim*n + d];
    }
  }
  /* Average */
  for (d=0; d<ndim; d++) {
    x_bar[d] = x_bar[d] / npoints;
  }
}

void compute_covariance_matrix_3d(long int npoints, double data[], double cov_matrix[][3])
{
  int   i,j,k;
  int   nvariables;
  double x_bar[3];

  /* 3 variables: x,y,z spatial directions */
  nvariables = 3;

  /* get the mean of each variable */
  get_mean(npoints,nvariables,data,x_bar);

  memset(cov_matrix, 0.0, (sizeof(double)) * 9 );
  /* compute covariance matrix */
  for (i=0; i<nvariables; i++) {
    for (j=0; j<nvariables; j++) {
      for (k=0; k<npoints; k++) {
        cov_matrix[i][j] += (data[nvariables*k + i] - x_bar[i])*(data[nvariables*k + j] - x_bar[j]) / npoints;
      }
    }
  }
}

void compute_covariance_eigenvectors(
  long int npoints,
  double point_data[],
  double eig_val[],
  double eig_vec[][3]
)
{
  int i,j,ierr;
  int sorted_idx[3];
  double cov_matrix[3][3],w[3];
  double complex Q[3][3],A[3][3];

  /* compute the covariance matrix of this set of points */
  compute_covariance_matrix_3d(npoints,point_data,cov_matrix);
  /* typecast it for the eigen computation */
  for (i=0; i<NSD_3D; i++) {
    for (j=0; j<NSD_3D; j++) {
      A[i][j] = (double complex)cov_matrix[i][j];
    }
  }
  /* compute the eigen vectors and eigen values of the covariance matrix */
  ierr = zheevj3(A, Q, w);

  /* Determine index of the sorted eigen values in ascending order such that
      w[idx[0]] <= w[idx[1]] <= w[idx[2]]
  */
  sort_ascendant(w,sorted_idx);
  /* Sort eigen values and vectors */
  for (i=0; i<NSD_3D; i++) {
    eig_val[i] = w[sorted_idx[i]];
    for (j=0; j<NSD_3D; j++) {
      eig_vec[i][j] = (double)Q[sorted_idx[i]][j];
    }
  }
  return;
}

void ft_compute_covariance_matrix_3d(long int npoints, double data[], double cov_v[])
{
  int   i,j,cnt;
  double cov_matrix[3][3];

  compute_covariance_matrix_3d(npoints,data,cov_matrix);
  cnt = 0;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      cov_v[cnt] = cov_matrix[i][j];
      cnt++;
    }
  }
}
