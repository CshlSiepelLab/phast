/* matrix - basic matrix and vector functions, mostly implemented as wrappers for CLAPACK and CBLAS routines */

/* $Id: matrix.c,v 1.1.1.1 2004-06-03 22:43:11 acs Exp $ 
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California 
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <matrix.h>
#include <f2c.h>
#include <clapack.h>
#include <gsl/gsl_cblas.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <misc.h>

#define EQ_THRESHOLD 1e-10

int mat_diagonalize(gsl_matrix *M, gsl_vector_complex *eval, 
                    gsl_matrix_complex *revect, 
                    gsl_matrix_complex *levect) { 

  char jobvl = 'V', jobvr = 'V';
  long int n = M->size1, lwork = 4 * M->size1, info; 
  double wr[n],wi[n], vr[n * n], vl[n * n], work[4 * n];
  int i, j, k;
  double tmp[n*n];
  enum {REAL, CONJ1, CONJ2} eval_type;
  int check = 0;

  assert(n == M->size2);

  /* convert matrix to representation used by LAPACK (column-major) */
  for (j = 0; j < n; j++) 
    for (i = 0; i < n; i++) 
      tmp[j*n + i] = gsl_matrix_get(M, i, j);

  dgeev_(&jobvl, &jobvr, &n, (doublereal*)tmp, &n, wr, wi, vl, 
         &n, vr, &n, work, &lwork, &info);

  if (info != 0) {
    fprintf(stderr, "ERROR executing the LAPACK 'dgeev' routine.\n"); 
    return 1;
  }

  for (j = 0; j < n; j++) {     /* column */
    gsl_complex z;
    
    GSL_SET_REAL(&z, wr[j]);
    GSL_SET_IMAG(&z, wi[j]);
    gsl_vector_complex_set(eval, j, z);

    /* if there are repeated eigenvalues, diagonalization may fail */
    for (k = 0; !check && k < j; k++) 
      if (GSL_COMPLEX_EQ(gsl_vector_complex_get(eval, k), z)) 
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
      gsl_complex zr, zl;
      switch (eval_type) {
      case REAL:
        GSL_SET_REAL(&zr, vr[j*n + i]);
        GSL_SET_IMAG(&zr, 0);
        GSL_SET_REAL(&zl, vl[j*n + i]);
        GSL_SET_IMAG(&zl, 0);
        break;
      case CONJ1:
        GSL_SET_REAL(&zr, vr[j*n + i]);
        GSL_SET_IMAG(&zr, vr[(j+1)*n + i]);
        GSL_SET_REAL(&zl, vl[j*n + i]);
        GSL_SET_IMAG(&zl, -vl[(j+1)*n + i]);
                                /* left eigenvectors according to
                                   dgeev are complex conjugates of
                                   those described by Press et al. --
                                   see dgeev man page */
        break;
      case CONJ2:
        GSL_SET_REAL(&zr, vr[(j-1)*n + i]);
        GSL_SET_IMAG(&zr, -vr[j*n + i]);
        GSL_SET_REAL(&zl, vl[(j-1)*n + i]);
        GSL_SET_IMAG(&zl, vl[j*n + i]);
        break;
      }

      gsl_matrix_complex_set(revect, i, j, zr);
      gsl_matrix_complex_set(levect, j, i, zl);
    }
  }

  /* rescale such that revect and levect are inverses */
  /* see Press et al. (Numerical Recipes) for a very clear explanation
     of the relationship between the left and right eigenvectors  */
  for (i = 0; i < n; i++) {
    /* compute dot product of row i in levect and col i in revect */
    gsl_complex dotprod;
    GSL_SET_REAL(&dotprod, 0);
    GSL_SET_IMAG(&dotprod, 0);
    for (j = 0; j < n; j++) 
      dotprod = gsl_complex_add(dotprod, gsl_complex_mul(gsl_matrix_complex_get(levect, i, j), gsl_matrix_complex_get(revect, j, i)));
      
    /* now scale levect accordingly */
    for (j = 0; j < n; j++) {
      gsl_complex oldval = gsl_matrix_complex_get(levect, i, j);
      gsl_complex scaled_val =  gsl_complex_div(oldval, dotprod);
      gsl_matrix_complex_set(levect, i, j, scaled_val);
    }
  }  

  /* verify correctness, if necessary */
  if (check) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        gsl_complex z; GSL_SET_REAL(&z, 0); GSL_SET_IMAG(&z, 0);
        for (k = 0; k < n; k++) 
          z = gsl_complex_add(z, gsl_complex_mul(gsl_matrix_complex_get(revect, i, k), gsl_complex_mul(gsl_vector_complex_get(eval, k), gsl_matrix_complex_get(levect, k, j))));
        if (fabs(GSL_IMAG(z)) > EQ_THRESHOLD ||
            fabs(GSL_REAL(z) - gsl_matrix_get(M, i, j)) > EQ_THRESHOLD) {
          fprintf(stderr, "ERROR: diagonalization failed.\n");
          assert(0);
        }
      }
    }
  }

  return 0;
}

int mat_eigenvals(gsl_matrix *M, gsl_vector_complex *evals) {
  char jobvl = 'N';
  char jobvr = 'N';
  long int n = M->size1;
  double wr[n];
  double wi[n];
  double work[4 * n];
  long int lwork = 4 * n;
  long int info; 
  int i, j;
  double tmp[n][n];

  for (i = 0; i < n; i++) 
    for (j = 0; j < n; j++) 
      tmp[i][j] = gsl_matrix_get(M, j, i);

  assert(n == M->size2);

  dgeev_(&jobvl, &jobvr, &n, (doublereal*)tmp, &n, wr, wi, NULL, 
         &n, NULL, &n, work, &lwork, &info);

  if (info != 0) {
    fprintf(stderr, "ERROR executing the LAPACK 'dgeev' routine.\n"); 
    return 1;
  }

  for (i = 0; i < n; i++) {
    gsl_complex z;
    GSL_SET_REAL(&z, wr[i]);
    GSL_SET_IMAG(&z, wi[i]);
    gsl_vector_complex_set(evals, i, z);
  }

  return 0;
}                     

int mat_invert(gsl_matrix *M_inv, gsl_matrix *M) {
  int i, j;
  long int info, n = M->size1;
  long int ipiv[n];
  double tmp[n][n];
  long int lwork = n;
  double work[lwork];

  assert(M->size1 == M->size2 && M_inv->size1 == M_inv->size2 && 
         M->size1 == M_inv->size1);

  for (i = 0; i < n; i++) 
    for (j = 0; j < n; j++) 
      tmp[i][j] = gsl_matrix_get(M, j, i);

  dgetrf_(&n, &n, (doublereal*)tmp, &n, ipiv, &info);

  if (info != 0) {
    fprintf(stderr, "ERROR: unable to compute LU factorization of matrix (for matrix inversion); dgetrf returned value of %d.\n", (int)info); 
    return 1;
  }

  dgetri_(&n, (doublereal*)tmp, &n, ipiv, work, &lwork, &info);

  if (info != 0) {
    if (info > 0)
      fprintf(stderr, "ERROR: matrix is singular -- cannot invert.\n");
    else
      fprintf(stderr, "ERROR: unable to invert matrix.  Element %d had an illegal value (according to dgetri).\n", (int)info); 
    return 1;
  }

  for (i = 0; i < M->size1; i++) 
    for (j = 0; j < M->size1; j++) 
      gsl_matrix_set(M_inv, i, j, tmp[j][i]);

  return 0;
}

void mat_mult(gsl_matrix *product, gsl_matrix *M1, gsl_matrix *M2) {
  assert(M1->size2 == M2->size1 && M1->size1 == M2->size2 && 
         product->size1 == M1->size1 && product->size2 == M2->size2);

  /* below adapted from call in GSL source code (gsl_blas_dgemm) */
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, product->size1, 
               product->size2, M1->size2, 1, M1->data, M1->tda, M2->data, 
               M2->tda, 0, product->data, product->tda);
}

void mat_vect_mult(gsl_vector *product, gsl_matrix *M, gsl_vector *v) {
  assert(M->size1 == v->size && v->size == product->size);

  cblas_dgemv (CblasRowMajor, CblasNoTrans, M->size1, 
               M->size2, 1, M->data, M->tda, v->data, 
               1, 0, product->data, 1);
}

void mat_mult_complex(gsl_matrix_complex *product, gsl_matrix_complex *M1, 
                      gsl_matrix_complex *M2) {
  int i, j, k;

  assert(M1->size2 == M2->size1 && M1->size1 == M2->size2 && 
         product->size1 == M1->size1 && product->size2 == M2->size2);

  gsl_matrix_complex_set_zero(product);
  for (i = 0; i < product->size1; i++) {
    for (j = 0; j < product->size2; j++) {
      for (k = 0; k < M1->size2; k++) {
        gsl_complex newterm, newsum;
        newterm = gsl_complex_mul(gsl_matrix_complex_get(M1, i, k),
                                  gsl_matrix_complex_get(M2, k, j));
        newsum = gsl_complex_add(gsl_matrix_complex_get(product, i, j), 
                                 newterm);
        gsl_matrix_complex_set(product, i, j, newsum);
      }
    }
  }
}

void mat_mult_complex_to_real(gsl_matrix *product, gsl_matrix_complex *M1, 
                              gsl_matrix_complex *M2) {

  int i, j, k;

  assert(M1->size2 == M2->size1 && M1->size1 == M2->size2 && 
         product->size1 == M1->size1 && product->size2 == M2->size2);

  gsl_matrix_set_zero(product);
  for (i = 0; i < product->size1; i++) {
    for (j = 0; j < product->size2; j++) {
      gsl_complex z;
      GSL_SET_REAL(&z, 0);
      GSL_SET_IMAG(&z, 0);
      for (k = 0; k < M1->size2; k++) {
        z = gsl_complex_add(z, gsl_complex_mul(gsl_matrix_complex_get(M1, i, k),
                                               gsl_matrix_complex_get(M2, k, j)));
      }
      assert(GSL_IMAG(z) < MAT_REAL_EPSILON);
      gsl_matrix_set(product, i, j, GSL_REAL(z));
    }
  }
}

double vect_dot_prod(gsl_vector *vect1, gsl_vector *vect2) {
  int i;
  double retval = 0;
  assert(vect1->size == vect2->size);
  for (i = 0; i < vect1->size; i++)
    retval += (gsl_vector_get(vect1, i) * gsl_vector_get(vect2, i));
  return retval;
}

void vect_cross_prod(gsl_matrix *mat, gsl_vector *v1, gsl_vector *v2) {
  int i, j;
  assert(v1->size == v2->size && mat->size1 == v1->size && 
         mat->size2 == v1->size);
  for (i = 0; i < v1->size; i++) {
    for (j = 0; j < v1->size; j++) {
      gsl_matrix_set(mat, i, j, 
                     gsl_vector_get(v1, i) * gsl_vector_get(v2, j));
    }
  }
}

double vect_norm(gsl_vector *v) {
  double ss = 0;
  int i;
  for (i = 0; i < v->size; i++) {
    double element = gsl_vector_get(v, i);
    ss += element*element;
  }
  return sqrt(ss);
}

void gsl_matrix_pretty_print(FILE *F, gsl_matrix *M) {
  int i, j;
  char *formatstr = "%11.6f ";
  double min = INFTY;

  /* find minimum non-zero absolute value; if it is very small, then
     print with exponential notation */
  for (i = 0; i < M->size1; i++) {
    for (j = 0; j< M->size2; j++) {
      double val = fabs(gsl_matrix_get(M, i, j));
      if (val != 0 && val < min) min = val;
    }
  }
  if (min < 1e-3) formatstr = "%14.6e ";

  for (i = 0; i < M->size1; i++) {
    for (j = 0; j < M->size2; j++) {
      fprintf(F, formatstr, gsl_matrix_get(M, i, j));
    }
    fprintf(F, "\n");
  }
}

void gsl_matrix_complex_pretty_print(FILE *F, gsl_matrix_complex *M) {
  int i, j;
  for (i = 0; i < M->size1; i++) {
    for (j = 0; j < M->size2; j++) {
      gsl_complex z = gsl_matrix_complex_get(M, i, j);
      fprintf(F, "%11.6f ", GSL_REAL(z));
      fprintf(F, "%c ", GSL_IMAG(z) >= 0 ? '+' : '-');
      fprintf(F, "%8.6fi ", fabs(GSL_IMAG(z)));
    }
    fprintf(F, "\n");
  }
}

