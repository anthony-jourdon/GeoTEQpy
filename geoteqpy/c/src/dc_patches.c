/*
====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  ptatin3d-extract-faults-tools
  filename: dc_patches.c

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

static int get_random_int_bounds(int start, int end)
{
  int number;
  number = (rand() % (end - start + 1)) + start;
  return number;
}

static void compute_patches_radius(int npatches, float radius0, int N0, float D, int number[], float radius[])
{
  int n;

  for(n=0; n<npatches; n++) {
    radius[n] = pow(2,(npatches-n-1)) * radius0;
    number[n] = pow(pow(2,D),n) * N0;
  }
}

static void apply_dc_patches(long int npoints, int npatches, int number[], float radius[], float coor[], float Dc0[], float Dc[])
{
  int n,p,d,np,cnt;

  cnt = 0;
  for(n=0; n<npatches; n++) {
    for(p=0; p<number[n]; p++) {
      int   pidx;
      float centre[NSD_3D];

      srand(cnt);
      cnt++;
      /* generate a random point index */
      pidx = get_random_int_bounds(0,npoints);
      /* define sphere centre at coords of that point */
      for (d=0; d<NSD_3D; d++) {
        centre[d] = coor[3*pidx + d];
      }
      /* loop over points */
      for (np=0; np<npoints; np++) {
        float sep,r2;

        /* compute distance squared between the point and the centre of the sphere */
        sep = 0.0;
        for (d=0; d<NSD_3D; d++) {
          sep += (coor[3*np + d] - centre[d]) * (coor[3*np + d] - centre[d]);
        }
        r2 = radius[n] * radius[n];
        if (sep <= r2) { /* inside the sphere */
          /* set the value corresponding to the current patches */
          Dc[np] = Dc0[n];
        }
      }
    }
  }
}

void ft_apply_dc_patches(long int npoints, 
                         int npatches,
                         float coor[],
                         float radius0,
                         int N0,
                         float D,
                         float Dc0[],
                         float Dc[])
{

  float *radius;
  int   *number;

  radius = (float*)malloc(npatches * sizeof(float));
  number = (int*  )malloc(npatches * sizeof(int)  );

  compute_patches_radius(npatches,radius0,N0,D,number,radius);
  apply_dc_patches(npoints,npatches,number,radius,coor,Dc0,Dc);

}
