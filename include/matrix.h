/* matrix - basic matrix and vector functions, mostly implemented as wrappers for CLAPACK and CBLAS routines */

/* $Id: matrix.h,v 1.1.1.1 2004-06-03 22:43:11 acs Exp $ 
   Written by Adam Siepel, Spring 2002
   Copyright 2002, Adam Siepel, University of California 

   NOTE: Function-by-function documentation appears in the header file
   only and is formatted for use with "c2man".  See generated man
   pages in section 3.
*/

#ifndef MAT_H
#define MAT_H

#include <gsl/gsl_matrix.h>

/* For use in mat_mult_complex_to_real: when the imaginary component
   of a complex number is smaller than this threshold, the number will be
   considered real */
#define MAT_REAL_EPSILON 1e-10


/* Diagonalize a square, real, nonsymmetric matrix.  Computes vector
   of eigenvalues and matrices of right and left eigenvectors,
   normalized so that they are inverses.

   Returns 0 on success, 1 on failure. */ 
int mat_diagonalize (gsl_matrix *M, /* input matrix (n x n) */
                     gsl_vector_complex *eval, 
                                /* computed vector of eigenvectors --
                                   preallocate dim. n */
                     gsl_matrix_complex *revect, 
                                /* computed matrix of right
                                   eigenvectors (columns) --
                                   preallocate n x n */
                     gsl_matrix_complex *levect); 
                                /* computed matrix of left
                                   eigenvectors (rows) -- preallocate
                                   n x n  */

/* Compute eigenvalues only of square, real nonsymmetric matrix.
   Returns 0 on success, 1 on failure.  
*/ 
int mat_eigenvals(gsl_matrix *M, /* input matrix (n x n) */
                  gsl_vector_complex *evals); /* computed vector of
                                                 eigenvalues
                                                 (preallocate
                                                 n-dim) */

/* Invert square, real, nonsymmetric matrix.  Uses LU decomposition
   (LAPACK routines dgetrf and dgetri).

   Returns 0 on success, 1 on failure. */
int mat_invert(gsl_matrix *M_inv, gsl_matrix *M);

/* Multiply two real matrices.  This is a wrapper for CBLAS's
   dgemm.  */
void mat_mult(gsl_matrix *product, /* computed product (preallocate n x n) */
              gsl_matrix *M1,   /* first input matrix (n x m) */
              gsl_matrix *M2);  /* second input matrix (m x n) */

/* Multiply two complex matrices. */
void mat_mult_complex(gsl_matrix_complex *product, gsl_matrix_complex *M1, 
                      gsl_matrix_complex *M2);

/* Multiply two complex matrices such that their product is real */
void mat_mult_complex_to_real(gsl_matrix *product, gsl_matrix_complex *M1, 
                              gsl_matrix_complex *M2);

/* Multiply a real matrix by a real vector.  This is a wrapper for
   CBLAS's dgemv. */
void mat_vect_mult(gsl_vector *product, /* computed product
                                           (preallocate n-dim) */
                   gsl_matrix *M, /* input matrix (n x m) */
                   gsl_vector *v); /* input vector (n-dim) */

/* Compute and return dot (inner) product of two n-dimensional
   real-valued vectors. */
double vect_dot_prod(gsl_vector *v1, /* first input vector (n-dim) */
                     gsl_vector *v2); /* second input vector (n-dim) */

/* Compute and return cross (outer) product of two n-dimensional
   real-valued vectors. */
void vect_cross_prod(gsl_matrix *mat, /* computed cross-product matrix
                                         (preallocate n x n) */
                     gsl_vector *v1, /* first input vector (n-dim) */
                     gsl_vector *v2); /* second input vector (n-dim) */

/* Compute and return 2-norm of vector. */
double vect_norm(gsl_vector *v);

void gsl_matrix_pretty_print(FILE *F, gsl_matrix *M); 
void gsl_matrix_complex_pretty_print(FILE *F, gsl_matrix_complex *M);

#endif
