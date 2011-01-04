/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file eigen.h
   Eigenvalue computation and matrix diagonalization.  CLAPACK is used.
   @ingroup base
*/

#ifndef EIG_H
#define EIG_H

#include <matrix.h>
#include <complex_vector.h>
#include <complex_matrix.h>

/** Diagonalize a square, real, non-symmetric matrix.
  @param[in,out] M Input matrix to diagonalize (n x n)
  @param[out] eval Eigen values computed from diagonalization preallocate dimension n
  @param[out] revect Normalized matrix of right eigen vectors from diagonalization preallocate dimension (n x n)
  @param[out] levect Normalized matrix of left eigen vectors from diagonalization preallocate dimension (n x n)
  @result 0 on success, otherwise failure
*/
int mat_diagonalize(Matrix *M, Zvector *eval, Zmatrix *revect, Zmatrix *levect);

/** Compute eigenvalues only of square, real non-symmetric matrix.
  @param[in,out] M Input matrix to find eigen values from dimension (n x n)
  @param[out] eval Eigen values, preallocate dimension n
  @result 0 on success, otherwise failure
*/
int mat_eigenvals(Matrix *M, Zvector *evals);

#endif
