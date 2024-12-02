/*
====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  GeoTEQpy
  filename: cov_eigen_v_mesh.c

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


#define LOG_DEBUG 0

static int sort_ComparePoints(const void *dataA, const void *dataB)
{
  Points *pointA;
  Points *pointB;

  pointA = (Points*)dataA;
  pointB = (Points*)dataB;

  /* 
  Compare 2 elements of the array that is passed to qsort()
  Here we indicate that the sorting is done on the cell index 
  in ascending order (from min(cell index) to max(cell index))
  such that qsort() will produce an array sorted by the cell index
  */

  if (pointA->cell_index < pointB->cell_index) {
    return -1;
  } else if (pointA->cell_index > pointB->cell_index) {
    return 1;
  } else {
    return 0;
  }
}

static void sort_Points(const long int npoints, Points list[])
{
  size_t np;

  np = (size_t)npoints;
  /* sort the array list[] with the condition provided in sort_ComparePoints() */
  qsort( list, np, sizeof(Points), sort_ComparePoints );
}

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

  /* Compute the distance squared between 2 points */
  dist = 0.0;
  for (dim=0; dim<NSD_3D; dim++) {
    dist += (pt_A[dim] - pt_B[dim]) * (pt_A[dim] - pt_B[dim]);
  }
  return dist;
}

static inline bool is_in_sphere(double centre[], double point[], double r2)
{
  double dist2;

  dist2 = distance2(centre,point);
  return (dist2 <= r2);
}

static void GetChunckBoundsIndices(int cell_idx, int msize[], int chunck_bounds[], int patch_extent)
{
  int cidx2d,cell_index_k,cell_index_j,cell_index_i;
  int mx,my,mz;
  /* 
  i,j,k indices in which we will search points 
  Because cell size are approx. equal to the radius of the sphere
  we know that the points must be contained in neighbouring cells.
  However, that can be controlled by patch_extent > 1 
  to search in more than the direct neighbours cells
  */

  mx = msize[0];
  my = msize[1];
  mz = msize[2];

  /* get the i,j,k indices of the current cell */
  cell_index_k = cell_idx / (mx*my);
  cidx2d = cell_idx - cell_index_k*(mx*my);
  cell_index_j = cidx2d / mx;
  cell_index_i = cidx2d - cell_index_j * mx;

  /* define start and end i,j,k indices of the cells we want to search in */
  chunck_bounds[0] = cell_index_i - patch_extent;   //istart 
  chunck_bounds[1] = cell_index_i + patch_extent+1; //iend

  chunck_bounds[2] = cell_index_j - patch_extent;   //jstart 
  chunck_bounds[3] = cell_index_j + patch_extent+1; //jend 

  chunck_bounds[4] = cell_index_k - patch_extent;   //kstart
  chunck_bounds[5] = cell_index_k + patch_extent+1; //kend

  /* ensure we are in the mesh bounds */
  if (chunck_bounds[0] < 0) { chunck_bounds[0] = 0; } //istart
  if (chunck_bounds[2] < 0) { chunck_bounds[2] = 0; } //jstart 
  if (chunck_bounds[4] < 0) { chunck_bounds[4] = 0; } //kstart

  if (chunck_bounds[1] > mx) { chunck_bounds[1] = mx; } //iend
  if (chunck_bounds[3] > my) { chunck_bounds[3] = my; } //jend
  if (chunck_bounds[5] > mz) { chunck_bounds[5] = mz; } //kend

#if (LOG_DEBUG > 0)
  printf("cell[%d]: start[i,j,k] = [%d, %d, %d], end[i,j,k] = [%d, %d, %d]\n",cell_idx,chunck_bounds[0],chunck_bounds[2],chunck_bounds[4],chunck_bounds[1],chunck_bounds[3],chunck_bounds[5]);
#endif

}

static int GetMaxPointsPerChunk(int chunck_bounds[], int msize[], int offset[])
{
  int istart,jstart,kstart;
  int iend,jend,kend;
  int max_points,nppc;
  int i,j,k,cell_idx,mx,my;

  mx = msize[0]; // number of cells in x direction
  my = msize[1]; // number of cells in y direction

  /* i,j,k bounds of the cells chunck in which we count */
  istart = chunck_bounds[0];
  iend   = chunck_bounds[1];
  jstart = chunck_bounds[2];
  jend   = chunck_bounds[3];
  kstart = chunck_bounds[4];
  kend   = chunck_bounds[5];

  max_points = 0;
  for (k=kstart; k<kend; k++) {
    for (j=jstart; j<jend; j++) {
      for (i=istart; i<iend; i++) {
        /* cell index */
        cell_idx = i + j*mx + k*mx*my;
        /* number of point in the cell */
        nppc = offset[cell_idx+1] - offset[cell_idx];
#if (LOG_DEBUG > 0)
  printf("cell[%d]: points in cell = %d\n",cell_idx,nppc);
#endif
        /* increase count */
        max_points += nppc;
      }
    }
  }
  return max_points;
}

static void PointsSearch_NeighbourCells_Sphere(int chunck_bounds[], 
                                               int msize[], 
                                               int offset[],
                                               double p_coords[], 
                                               double centre[],
                                               double r2,
                                               Points *plist,
                                               int *in_sphere_idx, 
                                               int *npoints)
{
  int point_count,i,j,k;
  int mx,my;
  int istart,iend;
  int jstart,jend;
  int kstart,kend;

  mx = msize[0]; // number of cells in x direction
  my = msize[1]; // number of cells in y direction

  /* i,j,k bounds of the cells chunck in which we count */
  istart = chunck_bounds[0];
  iend   = chunck_bounds[1];
  jstart = chunck_bounds[2];
  jend   = chunck_bounds[3];
  kstart = chunck_bounds[4];
  kend   = chunck_bounds[5];

  /* Initialize counter */
  point_count = 0;
  /* Iterate over cell chunck */
  for (k=kstart; k<kend; k++) {
    for (j=jstart; j<jend; j++) {
      for (i=istart; i<iend; i++) {
        int chunck_cidx,nppc,pc;

        /* cell index */
        chunck_cidx = i + j*mx + k*mx*my;
        /* points per cell in this cell */
        nppc = offset[chunck_cidx+1] - offset[chunck_cidx];

#if (LOG_DEBUG > 0)
  printf("Cell[%d]: points in cell = %d\n",chunck_cidx,nppc);
#endif

        /* if there is no point move to next cell */
        if (nppc == 0) { continue; }

        /* iterate over points that are in this cell */
        for (pc=0; pc<nppc; pc++) {
          int    d,cell2point,pidx;
          double point_coor[3];

          /* Get the point index of the point contained in that cell */
          cell2point = offset[chunck_cidx] + pc;
          /* unsorted point index */
          pidx = plist[ cell2point ].point_index;
          /* Get coordinates of the point */
          for (d=0; d<NSD_3D; d++) {
            point_coor[d] = p_coords[3*pidx + d];
          }

#if (LOG_DEBUG > 0)
  printf("current point = %d, offset[%d] = %d, cell to points index = %d\n",pc,chunck_cidx,offset[chunck_cidx],cell2point);
#endif
          /* Test if we are in the sphere */
          if ( is_in_sphere(centre,point_coor,r2) ) {
            /* Store the index of the point in the sphere */
            in_sphere_idx[point_count] = cell2point;
#if (LOG_DEBUG > 0)
  printf("in_sphere_idx[%d] = %d\n",point_count,cell2point);
#endif
            /* Increase counter */
            point_count++;
          }
        }
      }
    }
  }
  /* return the number of points */
  *npoints = point_count;
}

static void ComputePointsetEigv(int cell_idx, 
                                int msize[], 
                                int offset[], 
                                double p_coords[],
                                double centre[], 
                                double r2,
                                int patch_extent,
                                Points *plist,
                                double eig_val[], 
                                double eig_vec[][3])
{
  int    i,j,d,nps,idx[3],ierr;
  int    max_point_per_chunck,point_count,chunck_bounds[6];
  int    *in_sphere_idx;
  double *data;
  double cov_matrix[3][3],w[3];
  double complex Q[3][3],A[3][3];

  /* i,j,k indices of the cells chunk in which we will search */
  GetChunckBoundsIndices(cell_idx,msize,chunck_bounds,patch_extent);

  /* Count how much point there are in the cell chunck we look in */
  max_point_per_chunck = GetMaxPointsPerChunk(chunck_bounds,msize,offset);

#if (LOG_DEBUG > 0)
  printf("Cell[%d]: max_point_per_chunck = %d\n",cell_idx,max_point_per_chunck);
#endif

  /* Allocate an array to store the indices of the points that are inside the sphere */
  in_sphere_idx = (int*)malloc( (max_point_per_chunck)*sizeof(int) ); 
  /* Search which points are in the sphere and fill the array containing their indices */
  PointsSearch_NeighbourCells_Sphere(chunck_bounds,msize,offset,p_coords,centre,r2,plist,in_sphere_idx,&point_count);

#if (LOG_DEBUG > 0)
  printf("Number of points in the sphere: %d\n",point_count);
#endif

  /* Allocate the array containing the coords of the points inside the sphere */
  data = (double*)malloc( (3*point_count)*sizeof(double) );
  for (nps=0; nps<point_count; nps++) {
    int pidx;
    /* unsorted point index (original list) */
    pidx = plist[ in_sphere_idx[ nps ] ].point_index;
    /* fill the data array on which we will compute the covariance matrix */
    for (d=0; d<NSD_3D; d++) {
      data[3*nps + d] = p_coords[3*pidx + d];
    }
#if (LOG_DEBUG > 0)
  printf("datapoint[%d]:  point[%d]: coords = [%f, %f, %f]\n",nps,in_sphere_idx[ nps ],data[3*nps + 0],data[3*nps + 1],data[3*nps + 2]);
#endif
  }
  /* Free the array containing the indices, we don't need it anymore */
  free(in_sphere_idx);

  /* compute the covariance matrix of this set of points */
  compute_covariance_matrix_3d(point_count,data,cov_matrix);
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
  sort_ascendant(w,idx);
  /* Sort eigen values and vectors */
  for (i=0; i<NSD_3D; i++) {
    eig_val[i] = w[idx[i]];
    for (j=0; j<NSD_3D; j++) {
      eig_vec[i][j] = (double)Q[idx[i]][j];
    }
  }
  free(data);
}

static void ComputeCovEigVectors(float radius, 
                                 long int npoints, 
                                 int msize[], 
                                 int offset[],
                                 double p_coords[], 
                                 Points *plist, 
                                 int patch_extent,
                                 double e_vectors[])
{
  int    i,j,p,d,cnt;
  double r2;

  if (patch_extent < 0) { 
    printf("patch_extent < 0 is not possible, the search must be at least in the current cell\n");
    exit(1);
  }

  r2 = radius*radius;

  for (p=0; p<npoints; p++) {
    int    cell_idx,pidx;
    double centre[3],eig_val[3];
    double eig_vec[3][3];
    
    /* cell in which the point is contained */
    cell_idx = plist[p].cell_index;
    /* unsorted point index (original list) */
    pidx = plist[p].point_index;

    /* centre of the sphere we are looking in */
    for (d=0; d<NSD_3D; d++) {
      centre[d] = p_coords[3*pidx + d];
    }
    
#if (LOG_DEBUG > 0)
  printf("****** [[ Point %d, Cell %d ]] ******\n",p,cell_idx);
#endif
    ComputePointsetEigv(cell_idx,msize,offset,p_coords,centre,r2,patch_extent,plist,eig_val,eig_vec);

    cnt = 0;
    for (i=0; i<NSD_3D; i++) {
      for (j=0; j<NSD_3D; j++) {
        e_vectors[9*pidx + cnt] = eig_vec[i][j];
        cnt++;
      }
    }
  }
}

static void CreatePoint2CellConnectivity(long int npoints, int ncells, double p_coords[], int p_cellidx[], int pcell_list[], Points *plist)
{
  int p,c,count,tmp;

  for (p=0; p<npoints; p++) {
    plist[p].point_index     = p;
    plist[p].cell_index      = p_cellidx[p]; 
  }

  /* Sort the points in cell index ascending order */
  sort_Points(npoints,plist);

  /* sum points per cell */
  memset( pcell_list, 0, (ncells+1)*sizeof(int) ); // Initialize to 0
  for (p=0; p<npoints; p++) {
    pcell_list[ plist[p].cell_index ]++;

#if (LOG_DEBUG > 0)
    printf("plist[%d]: cell_index = %d, point_index = %d, coords = [%f, %f, %f]\n",p,plist[p].cell_index,plist[p].point_index,p_coords[3*p],p_coords[3*p+1],p_coords[3*p+2]);
#endif
  }

  /* 
  Create offset list 
  This list contains the sum of the points from cell index 0 to cell index c
  */
  count = 0;
  for (c=0; c<ncells; c++) {
    tmp = pcell_list[c];
    pcell_list[c] = count;
    count += tmp;
  }
  pcell_list[c] = count;

#if (LOG_DEBUG > 0)
  printf("****** offset list ******\n");
  for (c=0; c<ncells; c++) {
    printf("cell[%d]: offset = %d, points in cell = %d\n",c,pcell_list[c],pcell_list[c+1]-pcell_list[c]);
  }
#endif

}

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

  CreatePoint2CellConnectivity(npoints,ncells,p_coords,p_cellidx,pcell_list,plist);

  ComputeCovEigVectors(radius,npoints,msize,pcell_list,p_coords,plist,1,e_vectors);

#if (LOG_DEBUG > 0)
  printf("done\n");
#endif

  free(plist);
  free(pcell_list);
}
