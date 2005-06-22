/* $Id: eigen.c,v 1.1 2005-06-22 07:11:19 acs Exp $
   Written by Adam Siepel, Summer 2005
   Copyright 2005, Adam Siepel, University of California
*/

/** \file eigen.c
    Eigenvalue computation and matrix diagonalization.  CLAPACK is used.
    \ingroup base
*/

#include <eigen.h>
#include <assert.h>
#include <math.h>

#ifndef SKIP_CLAPACK 
#include <f2c.h>  
#include <clapack.h> 
#endif 

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

#ifdef SKIP_CLAPACK
  fprintf(stderr, "ERROR: CLAPACK required for matrix diagonalization.\n");
  assert(0);
#else
  char jobvl = 'V', jobvr = 'V';
  long int n = M->nrows, lwork = 4 * M->nrows, info; 
  double wr[n],wi[n], vr[n * n], vl[n * n], work[4 * n];
  int i, j, k;
  double tmp[n*n];
  enum {REAL, CONJ1, CONJ2} eval_type;
  int check = 0;

  assert(n == M->ncols);

  /* convert matrix to representation used by LAPACK (column-major) */
  for (j = 0; j < n; j++) 
    for (i = 0; i < n; i++) 
      tmp[j*n + i] = mat_get(M, i, j);

  dgeev_(&jobvl, &jobvr, &n, (doublereal*)tmp, &n, wr, wi, vl, 
         &n, vr, &n, work, &lwork, &info);

  if (info != 0) {
    fprintf(stderr, "ERROR executing the LAPACK 'dgeev' routine.\n"); 
    return 1;
  }

  for (j = 0; j < n; j++) {     /* column */
    Complex z = z_set(wr[j], wi[j]);
    zvec_set(eval, j, z);

    /* if there are repeated eigenvalues, diagonalization may fail */
    for (k = 0; !check && k < j; k++) 
      if (z_abs(z_sub(zvec_get(eval, k), z)) < EQ_THRESHOLD)
        check = 1;

    if (wi[j] == 0)             /* real eigenvalue */
      eval_type = REAL;
    else if (j < n-1 && wr[j] == wr[j+1] && wi[j] == -wi[j+1]) 
      eval_type = CONJ1;        /* first in conjugate pair */
    else if (j > 0 && wr[j-1] == wr[j] && wi[j-1] == -wi[j]) 
      eval_type = CONJ2;        /* second in conjugate pair */
    else
      assert(0);

    for (i = 0; i < n; i++) {   /* row */
      Complex zr, zl;
      switch (eval_type) {
      case REAL:
	zr = z_set(vr[j*n + i], 0);
	zl = z_set(vl[j*n + i], 0);
        break;
      case CONJ1:
	zr = z_set(vr[j*n + i], vr[(j+1)*n + i]);
        zl = z_set(vl[j*n + i], -vl[(j+1)*n + i]);
                                /* left eigenvectors according to
                                   dgeev are complex conjugates of
                                   those described by Press et al. --
                                   see dgeev man page */
        break;
      case CONJ2:
	zr = z_set(vr[(j-1)*n + i], -vr[j*n + i]);
        zl = z_set(vl[(j-1)*n + i], vl[j*n + i]);
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

        if (fabs(z.x) > EQ_THRESHOLD ||
            fabs(z.y - mat_get(M, i, j)) > EQ_THRESHOLD) {
          fprintf(stderr, "ERROR: diagonalization failed (got %f + %fi, expected %f).\n", 
		  z.x, z.y, mat_get(M, i, j));
	  assert(0);
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
#ifdef SKIP_CLAPACK
  fprintf(stderr, "ERROR: CLAPACK required for eigenvalue computation.\n");
  assert(0);
#else
  char jobvl = 'N';
  char jobvr = 'N';
  long int n = M->nrows;
  double wr[n];
  double wi[n];
  double work[4 * n];
  long int lwork = 4 * n;
  long int info; 
  int i, j;
  double tmp[n][n];

  for (i = 0; i < n; i++) 
    for (j = 0; j < n; j++) 
      tmp[i][j] = mat_get(M, j, i);

  assert(n == M->ncols);

  dgeev_(&jobvl, &jobvr, &n, (doublereal*)tmp, &n, wr, wi, NULL, 
         &n, NULL, &n, work, &lwork, &info);

  if (info != 0) {
    fprintf(stderr, "ERROR executing the LAPACK 'dgeev' routine.\n"); 
    return 1;
  }

  for (i = 0; i < n; i++) {
    Complex z = z_set(wr[i], wi[i]);
    zvec_set(evals, i, z);
  }

#endif
  return 0;
}                     

