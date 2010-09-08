/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell
 * University All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: eigen.c,v 1.5 2008-11-16 02:32:54 acs Exp $*/

/** \file eigen.c
    Eigenvalue computation and matrix diagonalization.  CLAPACK is used.
    \ingroup base
*/

#include <stdlib.h>
#include <eigen.h>
#include <math.h>
#include <external_libs.h>
#include <misc.h>


/* Diagonalize a square, real, nonsymmetric matrix.  Computes vector
   of eigenvalues and matrices of right and left eigenvectors,
   normalized so that they are inverses.  Returns 0 on success, 1 on
   failure. */ 
int mat_diagonalize(Matrix *M, /* input matrix (n x n) */
                    Zvector *eval, 
                                /* computed vector of eigenvectors --
                                   preallocate dim. n */
                    Zmatrix *revect, 
                                /* computed matrix of right
                                   eigenvectors (columns) --
                                   preallocate n x n */
                    Zmatrix *levect
                                /* computed matrix of left
                                   eigenvectors (rows) -- preallocate
                                   n x n  */
		    ) { 

#ifdef SKIP_LAPACK
  die("ERROR: LAPACK required for matrix diagonalization.\n");
#else
  char jobvl = 'V', jobvr = 'V';
  LAPACK_INT n = (LAPACK_INT)M->nrows, lwork = (LAPACK_INT)(4*M->nrows), info;
  LAPACK_DOUBLE tmp[n*n], wr[n], wi[n], vl[n * n], vr[n * n], work[4 * n];
  int i, j, k;
  enum {REALVAL, CONJ1, CONJ2} eval_type=REALVAL;
  int check = 0;
  
  if (n != M->ncols)
    die("ERROR in mat_diagonalize: M->nrows (%i) != M->ncols (%i)\n",
	M->nrows, M->ncols);

  /* convert matrix to representation used by LAPACK (column-major) */
  for (j = 0; j < n; j++) 
    for (i = 0; i < n; i++) 
      tmp[j*n + i] = (LAPACK_DOUBLE)mat_get(M, i, j);

#ifdef R_LAPACK
  F77_CALL(dgeev)(&jobvl, &jobvr, &n, tmp, &n, wr, wi, vl, &n,
		  vr, &n, work, &lwork, &info);
#else
  dgeev_(&jobvl, &jobvr, &n, tmp, &n, wr, wi, vl, 
         &n, vr, &n, work, &lwork, &info);
#endif

  if (info != 0) {
    fprintf(stderr, "ERROR executing the LAPACK 'dgeev' routine.\n"); 
    return 1;
  }

  for (j = 0; j < n; j++) {     /* column */
    Complex z = z_set((double)wr[j], (double)wi[j]);
    zvec_set(eval, j, z);

    /* if there are repeated eigenvalues, diagonalization may fail */
    for (k = 0; !check && k < j; k++) 
      if (z_abs(z_sub(zvec_get(eval, k), z)) < EQ_THRESHOLD)
        check = 1;

    if (wi[j] == 0)             /* real eigenvalue */
      eval_type = REALVAL;
    else if (j < n-1 && wr[j] == wr[j+1] && wi[j] == -wi[j+1]) 
      eval_type = CONJ1;        /* first in conjugate pair */
    else if (j > 0 && wr[j-1] == wr[j] && wi[j-1] == -wi[j]) 
      eval_type = CONJ2;        /* second in conjugate pair */
    else
      die("ERROR in mat_diagnoalize: complex eigenvalue does not have a conjugate pair");

    for (i = 0; i < n; i++) {   /* row */
      Complex zr, zl;
      switch (eval_type) {
      case REALVAL:
	zr = z_set((double)vr[j*n + i], 0);
	zl = z_set((double)vl[j*n + i], 0);
        break;
      case CONJ1:
	zr = z_set((double)vr[j*n + i], (double)vr[(j+1)*n + i]);
        zl = z_set((double)vl[j*n + i], -(double)vl[(j+1)*n + i]);
                                /* left eigenvectors according to
                                   dgeev are complex conjugates of
                                   those described by Press et al. --
                                   see dgeev man page */
        break;
      case CONJ2:
	zr = z_set((double)vr[(j-1)*n + i], -(double)vr[j*n + i]);
        zl = z_set((double)vl[(j-1)*n + i], (double)vl[j*n + i]);
        break;
      }

      zmat_set(revect, i, j, zr);
      zmat_set(levect, j, i, zl);
    }
  }

  /* rescale such that revect and levect are inverses */
  /* see Press et al. (Numerical Recipes) for a very clear explanation
     of the relationship between the left and right eigenvectors  */
  for (i = 0; i < n; i++) {
    /* compute dot product of row i in levect and col i in revect */
    Complex dotprod = z_set(0, 0);
    for (j = 0; j < n; j++) 
      dotprod = z_add(dotprod, z_mul(zmat_get(levect, i, j), zmat_get(revect, j, i)));
      
    /* now scale levect accordingly */
    for (j = 0; j < n; j++) {
      Complex oldval = zmat_get(levect, i, j);
      Complex scaled_val =  z_div(oldval, dotprod);
      zmat_set(levect, i, j, scaled_val);
    }
  }  

  /* verify correctness, if necessary */
  if (check) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        Complex z = z_set(0,0);
        for (k = 0; k < n; k++) 
          z = z_add(z, z_mul(zmat_get(revect, i, k), z_mul(zvec_get(eval, k), zmat_get(levect, k, j))));

        if (fabs(z.y) > EQ_THRESHOLD ||
            fabs(z.x - mat_get(M, i, j)) > EQ_THRESHOLD) {
          die("ERROR: diagonalization failed (got %e + %ei, expected %e).\n", 
		  z.x, z.y, mat_get(M, i, j));
	}
      }
    }
  }

#endif
  return 0;
}


/* Compute eigenvalues only of square, real nonsymmetric matrix.
   Returns 0 on success, 1 on failure.  
*/ 
int mat_eigenvals(Matrix *M, /* input matrix (n x n) */
		  Zvector *evals
				/* computed vector of eigenvalues
				   (preallocate n-dim) */
		  ) {
#ifdef SKIP_LAPACK
  die("ERROR: LAPACK required for eigenvalue computation.\n");
#else
  char jobvl = 'N';
  char jobvr = 'N';
  LAPACK_INT n = (LAPACK_INT)M->nrows, lwork = (LAPACK_INT)4*n, info;
  LAPACK_DOUBLE wr[n], wi[n], work[4 * n], tmp[n][n];
  int i, j;

  for (i = 0; i < n; i++) 
    for (j = 0; j < n; j++)
      tmp[i][j] = (LAPACK_DOUBLE)mat_get(M, j, i);

  if (n != M->ncols) 
    die("ERROR in mat_eigenvals: M->nrows (%i) != M->ncols (%i)\n",
	M->nrows, M->ncols);

#ifdef R_LAPACK
  F77_CALL(dgeev)(&jobvl, &jobvr, &n, (LAPACK_DOUBLE*)tmp, &n, wr, wi, NULL,
		  &n, NULL, &n, work, &lwork, &info);
#else
  dgeev_(&jobvl, &jobvr, &n, (LAPACK_DOUBLE*)tmp, &n, wr, wi, NULL, 
         &n, NULL, &n, work, &lwork, &info);
#endif

  if (info != 0) {
    fprintf(stderr, "ERROR executing the LAPACK 'dgeev' routine.\n"); 
    return 1;
  }

  for (i = 0; i < n; i++) {
    Complex z = z_set((double)wr[i], (double)wi[i]);
    zvec_set(evals, i, z);
  }

#endif
  return 0;
}                     

