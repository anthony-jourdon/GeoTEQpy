/*
====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  ptatin3d-extract-faults-tools
  filename: compute_medial_axis.c

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
#include <stdbool.h>
#include <math.h>
#include "kdtree.h"
#include "ft_constants.h"

#define FLOAT64_EQUALITY(a, b) ( (fabs((a)-(b)) < 1.0e-12) ? 1 : 0 )

typedef struct {
  double coor[NSD_3D];
} vector;

/* compute dot product between 2 vectors */
double dot(vector u, vector v) 
{
  int    d;
  double result = 0.0;
  
  for (d=0; d<NSD_3D; d++) {
    result += u.coor[d]*v.coor[d];
  }
  
  return result;
}

/* compute the norm of a vector */
double norm(vector u)
{
  double result = 0.0;
  
  result = sqrt(dot(u,u));
  
  return result;
}

double cosine_vec(vector p, vector q)
{
  double result;

  result = dot(p,q) / ( norm(p) * norm(q) );

  if (result > 1.0) { 
    return 1.0; 
  } else if (result < -1.0) {
    return -1.0;
  }
  
  return result;  
}

double compute_radius(vector p, vector q, vector n)
{
  int    d;
  double dist,cos_theta,radius;
  vector v;

  /* compute the vector linking points p and q */
  for (d=0; d<NSD_3D; d++) {
    v.coor[d] = p.coor[d] - q.coor[d];
  }
  /* distance between p and q */
  dist = norm(v);
  /* cosine of the angle between the normal and the vector v */
  cos_theta = dot(n,v) / dist;
  /* trigonometry to obtain the radius from a chord of the circle */
  radius = fabs(dist / (2.0 * cos_theta));

  return radius;
}

/* 
Check equality between vectors, to be equal they must:
  - have the same norm
  - have the same components in the same order
*/
bool vector_equality(vector u, vector v)
{
  if ( FLOAT64_EQUALITY( norm(u), norm(v) ) ) {
    if ( FLOAT64_EQUALITY( u.coor[0], v.coor[0] ) ) {
      if ( FLOAT64_EQUALITY( u.coor[1], v.coor[1] ) ) {
        if ( FLOAT64_EQUALITY( u.coor[2], v.coor[2] ) ) {
          return 1;
        } else {
          return 0;
        }
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

bool noise_handling(vector p, vector q, vector ball_centre, int it, double angle_planar, double angle_preserve, double radius)
{
  int    d;
  double separation_angle;
  vector qc,pc,pq;

  for (d=0; d<NSD_3D; d++) {
    qc.coor[d] = q.coor[d] - ball_centre.coor[d];
    pc.coor[d] = p.coor[d] - ball_centre.coor[d];
    pq.coor[d] = p.coor[d] - q.coor[d];
  }

  separation_angle = acos( cosine_vec(pc, qc) );

  if (it == 0 && angle_planar > separation_angle) {
    return 1;
  }

  if (it > 0 && angle_preserve > separation_angle && radius > norm(pq)) {
    return 1;
  }

  return 0;
}

vector compute_medial_axis(KDTree kdtree,
                           vector p,
                           vector n,
                           int max_it,
                           double radius_init,
                           double Tol,
                           double angle_planar,
                           double angle_preserve)
{
  int    d,it;
  double radius,radius_prev;
  vector ma,ball_centre;

  /* initialize the radius of the ball */
  radius_prev = radius_init;

  it = 0;
  while (it < max_it) {
    kd_node  nearest;
    double   nearest_dist;
    vector   q;

    /* compute the position of the ball centre */
    for (d=0; d<NSD_3D; d++) {
      ball_centre.coor[d] = p.coor[d] - n.coor[d] * radius_prev;
    }
    /* assign ball centre coordinates to find the nearest point */
    KDTreeFindNearest(kdtree,ball_centre.coor,&nearest,&nearest_dist);
    /* assign the nearest point coordinates */
    for (d=0; d<NSD_3D; d++) {
      q.coor[d] = nearest->x[d];
    }

    /* if the radius is already smaller than the distance to the closest point => stop */
    if ( nearest_dist >= radius_prev ) {
      //printf("point[%d]: nearest distance (%f) >= radius (%f)\n",nearest->index,nearest_dist,radius_prev);
      break;
    }
    /* if the point q is the same than point p => stop */
    if ( vector_equality(p, q) ) {
      break;
    }

    /* compute the new radius of the ball */
    radius = compute_radius(p, q, n);
    /* if the radius did not decrease it means that we cannot find a smaller ball => stop */
    if (fabs(radius_prev - radius) < Tol) {
      //printf("Converged! previous radius - radius < Tol: %f - %f < %f\n",radius_prev,radius,Tol);
      break;
    }

    if ( noise_handling(p,q,ball_centre,it,angle_planar,angle_preserve,radius) ) {
      break;
    }

    /* update radius */
    radius_prev = radius;
    /* increase iteration count */
    it++;
  }
  ma = ball_centre;

  return ma;
}

void compute_medial_axis_transform(double points[],
                                   double normals[],
                                   int npoints,
                                   double radius_init,
                                   double angle_planar,
                                   double angle_preserve,
                                   double medial_axis_points[])
{
  int     np,d;
  KDTree  kdtree;
  kd_node node;

  /* allocate kdtree data structure */
  KDTreeCreate(3,&kdtree);
  /* set points */
  KDTreeSetPoints(kdtree,npoints);
  /* fill kdtree points */
  KDTreeGetPoints(kdtree,NULL,&node);
  for (np=0; np<npoints; np++) {
    /* assign coordinates */
    for (d=0; d<NSD_3D; d++) {
      node[np].x[d] = points[NSD_3D*np + d];
    }
    /* assign indices */
    node[np].index = np;
  }
  KDTreeSetup(kdtree);

  /* loop over points and apply the shrinking ball algorithm */
  for (np=0; np<npoints; np++) {
    vector p,n,ma;
    
    /* current point and its normal */
    for (d=0; d<NSD_3D; d++) {
      p.coor[d] = points[  NSD_3D*np + d ];
      n.coor[d] = normals[ NSD_3D*np + d ];
    }
    /* compute the medial axis point with the shrinking ball algorithm */
    ma = compute_medial_axis(kdtree,p,n,20,radius_init,1.0e-12,angle_planar,angle_preserve);
    /* fill the output array */
    for (d=0; d<NSD_3D; d++) {
      medial_axis_points[ NSD_3D*np + d ] = ma.coor[d];
    }
  }
  KDTreeDestroy(&kdtree);
}
