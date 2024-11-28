/*
====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  ptatin3d-extract-faults-tools
  filename: eigen_v.c

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

#include "faulttools.h"

#define SQR(x)      ((x)*(x))                        // x^2
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  // |x|^2

#define FLOAT64_EQUALITY(a, b) ( (fabs((a)-(b)) < 1.0e-20) ? 1 : 0 )

// ----------------------------------------------------------------------------
int zheevj3(double complex A[3][3], double complex Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a hermitian 3x3
// matrix A using the Jacobi algorithm.
// The upper triangular part of A is destroyed during the calculation,
// the diagonal elements are read but not destroyed, and the lower
// triangular elements are not referenced at all.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
{
  const int n = 3;
  double sd, so;                  // Sums of diagonal resp. off-diagonal elements
  double complex s, t;            // sin(phi), tan(phi) and temporary storage
  double c;                       // cos(phi)
  double g, h, z;                 // Temporary storage
  double thresh;
  int i,nIter,p,q,r;
  
  // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
  for (i=0; i < n; i++) {
    Q[i][i] = 1.0;
    for (int j=0; j < i; j++)
    Q[i][j] = Q[j][i] = 0.0;
  }
#endif
  
  // Initialize w to diag(A)
  for (i=0; i < n; i++)
  w[i] = creal(A[i][i]);
  
  // Calculate SQR(tr(A))
  sd = 0.0;
  for (i=0; i < n; i++)
  sd += fabs(w[i]);
  sd = SQR(sd);
  
  // Main iteration loop
  for (nIter=0; nIter < 50; nIter++) {
    // Test for convergence
    so = 0.0;
    for (p=0; p < n; p++)
    for (q=p+1; q < n; q++)
    so += fabs(creal(A[p][q])) + fabs(cimag(A[p][q]));
    if (so == 0.0)
    return(0);
    
    if (nIter < 4)
    thresh = 0.2 * so / SQR(n);
    else
    thresh = 0.0;
    
    // Do sweep
    for (p=0; p < n; p++)
    for (q=p+1; q < n; q++) {
      g = 100.0 * (fabs(creal(A[p][q])) + fabs(cimag(A[p][q])));
      if (nIter > 4  &&  fabs(w[p]) + g == fabs(w[p])
          &&  fabs(w[q]) + g == fabs(w[q])) {
        A[p][q] = 0.0;
      }
      else if (fabs(creal(A[p][q])) + fabs(cimag(A[p][q])) > thresh) {
        // Calculate Jacobi transformation
        h = w[q] - w[p];
        if (fabs(h) + g == fabs(h))
        t = A[p][q] / h;
        else {
          if (h < 0.0)
          t = -2.0 * A[p][q] / (sqrt(SQR(h) + 4.0*SQR_ABS(A[p][q])) - h);
          else if (h == 0.0)
          t = A[p][q] * (1.0 / cabs(A[p][q]));  // A[p][q]/fabs(A[p][q]) could cause overflows
          else
          t = 2.0 * A[p][q] / (sqrt(SQR(h) + 4.0*SQR_ABS(A[p][q])) + h);
        }
        c = 1.0/sqrt(1.0 + SQR_ABS(t));
        s = t * c;
        z = creal(t * conj(A[p][q]));
        
        // Apply Jacobi transformation
        A[p][q] = 0.0;
        w[p] -= z;
        w[q] += z;
        for (r=0; r < p; r++) {
          t = A[r][p];
          A[r][p] = c*t - conj(s)*A[r][q];
          A[r][q] = s*t + c*A[r][q];
        }
        for (r=p+1; r < q; r++) {
          t = A[p][r];
          A[p][r] = c*t - s*conj(A[r][q]);
          A[r][q] = s*conj(t) + c*A[r][q];
        }
        for (r=q+1; r < n; r++) {
          t = A[p][r];
          A[p][r] = c*t - s*A[q][r];
          A[q][r] = conj(s)*t + c*A[q][r];
        }
        
        // Update eigenvectors
#ifndef EVALS_ONLY
        for (r=0; r < n; r++) {
          t = Q[r][p];
          Q[r][p] = c*t - conj(s)*Q[r][q];
          Q[r][q] = s*t + c*Q[r][q];
        }
#endif
      }
    }
  }
  
  return(-1);
}

void compute_eigV(long int n_points, double E[], double eig_vec[], double eig_val[])
{
  long int       n;
  int            i,j,cnt,ierr;
  double complex A[3][3],Q[3][3];
  double         w[3];
  
  for (n=0; n<n_points; n++){
    /* Construct the 3x3 tensor from 1D array */
    cnt = 0;
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        A[i][j] = E[9*n + cnt];
        cnt++;
      }
    }
    
    /* Compute eigen values and vectors */
    ierr = zheevj3(A, Q, w);
    
    cnt = 0;
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        eig_vec[9*n + cnt] = (double)Q[i][j];
        cnt++;
      }
    }
    
    cnt = 0;
    for (i=0; i<3; i++){
      eig_val[3*n + cnt] = w[i];
      cnt++;
    }
  }
}

void compute_eigV_sorted(long int n_points, double E[], double eig_vec[], double eig_val[])
{
  long int       n,idx[3];
  int            i,j,cnt,ierr;
  double complex A[3][3],Q[3][3];
  double         val[3],vec[3][3];
  double         w[3];
  
  for (n=0; n<n_points; n++){
    /* Construct the 3x3 tensor from 1D array */
    cnt = 0;
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        A[i][j] = E[9*n + cnt];
        cnt++;
      }
    }
    
    /* Compute eigen values and vectors */
    ierr = zheevj3(A, Q, w);

    /* Determine index of the sorted eigen values in ascending order such that
       w[idx[0]] <= w[idx[1]] <= w[idx[2]]
    */
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

    /* Sort eigen values and vectors */
    for (i=0; i<3; i++) {
      val[i] = w[idx[i]];
      for (j=0; j<3; j++) {
        vec[i][j] = (double)Q[idx[i]][j];
      }
    }
    
    cnt = 0;
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        eig_vec[9*n + cnt] = vec[i][j];
        cnt++;
      }
    }
    
    cnt = 0;
    for (i=0; i<3; i++){
      eig_val[3*n + cnt] = val[i];
      cnt++;
    }
  }
}

void ft_compute_eigV(long int n_points, double E[], double eig_vec[], double eig_val[])
{
  compute_eigV(n_points,E,eig_vec,eig_val);
}

void ft_compute_eigV_sorted(long int n_points, double E[], double eig_vec[], double eig_val[])
{
  compute_eigV_sorted(n_points,E,eig_vec,eig_val);
}
