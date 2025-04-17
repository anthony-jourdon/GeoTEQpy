/*
====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  GeoTEQpy
  filename: pca.c

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
#include "SortPCtx.h"

void sphere_cov_eig_vectors(long int npoints, 
  int ncells, 
  int msize[], 
  double p_coords[], 
  int p_cellidx[], 
  float radius, 
  double e_vectors[])
{
Points *plist;
int    *pcell_list;

pcell_list = (int*)malloc( (ncells+1)*sizeof(int) );
/* Allocate the data structure that will contain the points information */
plist      = (Points*)malloc( (npoints)*sizeof(Points) );

compute_covariance_eigenvectors_sphere(radius,npoints,ncells,p_cellidx,pcell_list,msize,pcell_list,p_coords,plist,1,e_vectors);

free(plist);
free(pcell_list);
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
  compute_covariance_eigenvectors_k_nearest(kdtree,k_nearest,npoints,p_coords,e_vectors);
  /* free kdtree */
  KDTreeDestroy(&kdtree);
}
