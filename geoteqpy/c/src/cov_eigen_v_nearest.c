/*
====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  GeoTEQpy
  filename: cov_eigen_v_nearest.c

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
#include <stdbool.h>

#include "faulttools.h"
#include "kdtree.h"

#define LOG_DEBUG 0

static void ComputeCovEigVectors(
  KDTree kdtree, 
  long int k_nearest, 
  long int npoints, 
  double points[], 
  double e_vectors[]
)
{
  maxheap  *heap;
  long int np;
  double   *nearest_points;
  double   eig_val[3],eig_vec[3][3];

  MaxHeapCreate(k_nearest,&heap);
  nearest_points = (double *)malloc(sizeof(double)*NSD_3D*k_nearest);

  for (np=0; np<npoints; np++) {
    int     i,j,k,cnt;
    double  *point_coor = &points[NSD_3D*np];
    /* Find the k nearest neighbours of that point */
    KDTreeFindKNearest(kdtree,point_coor,heap);

    /* get the nearest points */
    for (k=0; k<k_nearest; k++) {
      heap_node *node = &heap->nodes[k];
      int       idx = node->global_idx;
      /* Copy coordinates of nearest points in a dedicated array */
      memcpy(&nearest_points[NSD_3D*k],&points[NSD_3D*idx],sizeof(double)*NSD_3D);
    }

    /* Compute the covariance matrix eigenvectors and eigenvalues of these k nearest neighbour points */
    compute_covariance_eigenvectors(k_nearest,nearest_points,eig_val,eig_vec);

    cnt = 0;
    for (i=0; i<NSD_3D; i++) {
      for (j=0; j<NSD_3D; j++) {
        e_vectors[9*np + cnt] = eig_vec[i][j];
        cnt++;
      }
    }
    /* Reset the max heap for next point */
    MaxHeapReset(heap);
  }
  MaxHeapDestroy(heap);
  free(nearest_points);
  return;
}

void cov_eig_vectors_knearest(
  long int npoints,
  long int k_nearest,
  double p_coords[],
  double e_vectors[]
)
{
  KDTree   kdtree;
  kd_node  nodes;
  int      d;
  long int np;

  /* allocate kdtree data structure */
  KDTreeCreate(3,&kdtree);
  /* set points */
  KDTreeSetPoints(kdtree,npoints);
  /* fill kdtree points */
  KDTreeGetPoints(kdtree,NULL,&nodes);
  for (np=0; np<npoints; np++) {
    /* assign coordinates */
    for (d=0; d<NSD_3D; d++) {
      nodes[np].x[d] = p_coords[NSD_3D*np + d];
    }
    /* assign indices */
    nodes[np].index = np;
  }
  KDTreeSetup(kdtree);
  /* Compute the covariance matrix eigenvectors for each point */
  ComputeCovEigVectors(kdtree,k_nearest,npoints,p_coords,e_vectors);
  /* free kdtree */
  KDTreeDestroy(&kdtree);
}