/*
====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  ptatin3d-extract-faults-tools
  filename: cov_eigen_v.c

  This file is part of ptatin3d-extract-faults-tools.

  ptatin3d-extract-faults-tools is free software: you can redistribute it and/or modify it under the terms 
  of the GNU General Public License as published by the Free Software Foundation, either 
  version 3 of the License, or any later version.

  ptatin3d-extract-faults-tools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with ptatin3d-extract-faults-tools. 
  If not, see <https://www.gnu.org/licenses/>.
====================================================================================================
*/
#include <stdio.h>
#include "faulttools.h"

#define LOG_DEBUG 0

static inline void sort_ascendant(double w[], int idx[])
{
  if (w[0] <= w[1] && w[0] <= w[2]) {
    if (w[1] <= w[2]) {
      idx[0] = 0;
      idx[1] = 1;
      idx[2] = 2;
    } else {
      idx[0] = 0;
      idx[1] = 2;
      idx[2] = 1;
    }
  } else if (w[1] <= w[0] && w[1] <= w[2]) {
    if (w[0] <= w[2]) {
      idx[0] = 1;
      idx[1] = 0;
      idx[2] = 2;
    } else {
      idx[0] = 1;
      idx[1] = 2;
      idx[2] = 0;
    }
  } else {
    if (w[0] <= w[1]) {
      idx[0] = 2;
      idx[1] = 0;
      idx[2] = 1;
    } else {
      idx[0] = 2;
      idx[1] = 1;
      idx[2] = 0;
    }
  }
}

static inline double distance2(double pt_A[], double pt_B[])
{
  int    dim;
  double dist;

  dist = 0.0;
  for (dim=0; dim<NSD_3D; dim++) {
    dist += (pt_A[dim] - pt_B[dim]) * (pt_A[dim] - pt_B[dim]);
  }
  return dist;
}

static void _count(double centre[], double r2, double p_coor[],long int npoints, int *in_cnt)
{
  long int p;
  int cnt;

  cnt = 0;
  for (p=0; p<npoints; p++) {
    double point[3];
    double dist;

    point[0] = p_coor[3*p    ];
    point[1] = p_coor[3*p + 1];
    point[2] = p_coor[3*p + 2];
    dist  = distance2(centre,point);

    if (dist <= r2) {
      cnt++;
    }
  }
  *in_cnt = cnt;
}

static void is_in_sphere(double centre[], double r2, double p_coor[],long int npoints, int *indices)
{
  long int p;
  int      in_cnt;

  in_cnt = 0;
  for (p=0; p<npoints; p++) {
    double point[3];
    double dist;

    point[0] = p_coor[3*p    ];
    point[1] = p_coor[3*p + 1];
    point[2] = p_coor[3*p + 2];
    dist  = distance2(centre,point);

    if (dist <= r2) {
      indices[in_cnt] = p;
      in_cnt++;
    }
  }
}

static void compute_pointset_eigv(double centre[], 
                                  double r2, 
                                  double point_coor[], 
                                  long int npoints, 
                                  double eig_val[], 
                                  double eig_vec[][3])
{
  int            *indices,idx[3];
  int            in_cnt,np,i,j,d,ierr;
  double         *data;
  double         cov_matrix[3][3];
  double complex Q[3][3],A[3][3];
  double         w[3];

  /* count how many points are in the sphere */
  _count(centre,r2,point_coor,npoints,&in_cnt);
  /* allocate memory to store indices of these points */
  indices = (int *)malloc(in_cnt * sizeof(int));

  /* get their indices */
  is_in_sphere(centre,r2,point_coor,npoints,indices);
  /* allocate memory to store data */
  data = (double *)malloc(3 * in_cnt * sizeof(double));

#if (LOG_DEBUG > 0)
  printf("Number of points in sphere: %d\n",in_cnt);
#endif

  for (np=0; np<in_cnt; np++) {
    /* get coordinates of the points inside the sphere */
    for (d=0; d<NSD_3D; d++) {
      data[3*np + d] = point_coor[3*indices[np] + d];
    }
#if (LOG_DEBUG > 0)
  printf("data[%d]: point index = %d, coords = [%f, %f, %f]\n",np,indices[np],data[3*np + 0],data[3*np + 1],data[3*np + 2]);
#endif
  }
  
  /* compute the covariance matrix of this set of points */
  compute_covariance_matrix_3d(in_cnt,data,cov_matrix);
  for (i=0; i<NSD_3D; i++) {
    for (j=0; j<NSD_3D; j++) {
      A[i][j] = (double complex)cov_matrix[i][j];
#if (LOG_DEBUG > 0)
  printf("Cov matrix[%d][%d] = %f\n",i,j,cov_matrix[i][j]);
#endif
    }
  }

  /* compute the eigen vectors and eigen values of the covariance matrix */
  ierr = zheevj3(A, Q, w);

  /* Determine index of the sorted eigen values in ascending order such that
       w[idx[0]] <= w[idx[1]] <= w[idx[2]]
  */
  sort_ascendant(w,idx);
  /* Sort eigen values and vectors */
  for (i=0; i<NSD_3D; i++) {
    eig_val[i] = w[idx[i]];
    for (j=0; j<NSD_3D; j++) {
      eig_vec[i][j] = (double)Q[idx[i]][j];
#if (LOG_DEBUG > 0)
  printf("eig vec[%d][%d] = %f\n",i,j,eig_vec[i][j]);
#endif
    }
#if (LOG_DEBUG > 0)
  printf("eig val[%d] = %f\n",i,eig_val[i]);
#endif
  }
  free(indices);
  free(data);
}

static void compute_cov_eig_vectors(float radius, double point_coor[], long int npoints, double vectors[])
{
  int    p,d,i,j,cnt;
  double r2;

  r2 = radius*radius;
  for (p=0; p<npoints; p++) {
    double eig_val[3],eig_vec[3][3],centre[3];

    /* define the centre at current point */
    for (d=0; d<NSD_3D; d++) {
      centre[d] = point_coor[3*p + d];
    }
    /* 
    search points in the radius of that point, 
    compute its covariance matrix 
    and the eigen vectors of that matrix 
    */
#if (LOG_DEBUG > 0)
  printf("****** Point %d ******\n",p);
#endif
    compute_pointset_eigv(centre,r2,point_coor,npoints,eig_val,eig_vec);    
    
    cnt = 0;
    for (i=0; i<NSD_3D; i++) {
      for (j=0; j<NSD_3D; j++) {
        vectors[9*p + cnt] = eig_vec[i][j];
        cnt++;
      }
    }
  }
}

void ft_compute_cov_eig_vectors(float radius, double point_coor[], long int npoints, double vectors[])
{
  compute_cov_eig_vectors(radius,point_coor,npoints,vectors);
}
