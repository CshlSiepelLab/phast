/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: matrix.c,v 1.10 2008-11-16 02:43:24 acs Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <matrix.h>
#include <external_libs.h>
#include <math.h>
#include <misc.h>

Matrix *mat_new(int nrows, int ncols) {
  int i;
  Matrix *m = smalloc(sizeof(Matrix));
  m->data = smalloc(nrows * sizeof(void*));
  for (i = 0; i < nrows; i++)
    m->data[i] = smalloc(ncols * sizeof(double));
  m->nrows = nrows;
  m->ncols = ncols;
  return m;
}

Matrix *mat_new_from_array(double **array, int nrows, int ncols) {
  int i, j;
  Matrix *m = mat_new(nrows, ncols);
  for (i = 0; i < nrows; i++)
    for (j = 0; j < ncols; j++)
      m->data[i][j] = array[i][j];
  return m;
}

void mat_free(Matrix *m) {
  int i;
  for (i = 0; i < m->nrows; i++)
    sfree(m->data[i]);
  sfree(m->data);
  sfree(m);
}

Vector *mat_get_row(Matrix *m, int row) {
  int j;
  Vector *v = vec_new(m->ncols);
  for (j = 0; j < m->ncols; j++)
    v->data[j] = m->data[row][j];
  return v;
}

Vector *mat_get_col(Matrix *m, int col) {
  int i;
  Vector *v = vec_new(m->nrows);
  for (i = 0; i < m->nrows; i++)
    v->data[i] = m->data[i][col];
  return v;
}

void mat_set_identity(Matrix *m) {
  int i, j;
  for (i = 0; i < m->nrows; i++)
    for (j = 0; j < m->ncols; j++)
      m->data[i][j] = (i == j ? 1 : 0);
}

void mat_zero(Matrix *m) {
  int i, j;
  for (i = 0; i < m->nrows; i++)
    for (j = 0; j < m->ncols; j++)
      m->data[i][j] = 0;
}

void mat_set_all(Matrix *m, double val) {
  int i, j;
  for (i = 0; i < m->nrows; i++)
    for (j = 0; j < m->ncols; j++)
      m->data[i][j] = val;
}

void mat_copy(Matrix *dest, Matrix *src) {
  int i, j;
  if (dest->nrows != src->nrows)
    die("ERROR: mat_copy: dest->nrows (%i) != src->nrows (%i)\n",
	dest->nrows, src->nrows);
  if (dest->ncols != src->ncols)
    die("ERROR: mat_copy: dest->ncols (%i) != src->ncols (%i)\n",
	dest->ncols, src->ncols);
  for (i = 0; i < dest->nrows; i++)
    for (j = 0; j < dest->ncols; j++)
      dest->data[i][j] = src->data[i][j];
}

Matrix *mat_create_copy(Matrix *src) {
  Matrix *dest = mat_new(src->nrows, src->ncols);
  mat_copy(dest, src);
  return dest;
}

Matrix *mat_transpose(Matrix *src) {
  int i, j;
  Matrix *retval = mat_new(src->ncols, src->nrows);
  for (i = 0; i < src->nrows; i++)
    for (j = 0; j < src->ncols; j++)
      retval->data[j][i] = src->data[i][j];
  return retval;
}

void mat_scale(Matrix *m, double scale_factor) {
  int i, j;
  for (i = 0; i < m->nrows; i++)
    for (j = 0; j < m->ncols; j++)
      m->data[i][j] *= scale_factor;
}

void mat_print(Matrix *m, FILE *F) {
  int i, j;
  char *formatstr = "%11.6f ";
  double min = INFTY;

  /* find minimum non-zero absolute value; if it is very small, then
     print with exponential notation */
  for (i = 0; i < m->nrows; i++) {
    for (j = 0; j< m->ncols; j++) {
      double val = fabs(m->data[i][j]);
      if (val != 0 && val < min) min = val;
    }
  }
  if (min < 1e-3) formatstr = "%14.6e ";

  for (i = 0; i < m->nrows; i++) {
    for (j = 0; j < m->ncols; j++)
      fprintf(F, formatstr, m->data[i][j]);
    fprintf(F, "\n");
  }
}

void mat_read(Matrix *m, FILE *F) {
  int i, j;
  for (i = 0; i < m->nrows; i++)
    for (j = 0; j < m->ncols; j++)
      if (1 != fscanf(F, "%lf ", &m->data[i][j]))
	die("ERROR reading matrix");
}

Matrix *mat_new_from_file(FILE *F, int nrows, int ncols) {
  Matrix *m = mat_new(nrows, ncols);
  mat_read(m, F);
  return m;
}

void mat_mult(Matrix *prod, Matrix *m1, Matrix *m2) {
  if (!(m1->ncols == m2->nrows && m1->nrows == m2->ncols &&
	prod->nrows == m1->nrows && prod->ncols == m2->ncols))
    die("ERROR mat_mult: bad matrix dimensions\n");
  int i, j, k;
  for (i = 0; i < prod->nrows; i++) {
    for (j = 0; j < prod->ncols; j++) {
      prod->data[i][j] = 0;
      for (k = 0; k < m1->ncols; k++)
	prod->data[i][j] += m1->data[i][k] * m2->data[k][j];
    }
  }
}

void mat_vec_mult(Vector *prod, Matrix *m, Vector *v) {
  int i, j;
  if (!(m->nrows == prod->size && v->size == m->ncols))
    die("ERROR mat_vec_mult: bad dimensions\n");
  for (i = 0; i < m->nrows; i++) {
    prod->data[i] = 0;
    for (j = 0; j < m->ncols; j++) {
      prod->data[i] += m->data[i][j] * v->data[j];
    }
  }
}

void mat_plus_eq(Matrix *thism, Matrix *addm) {
  int i, j;
  if (!(thism->nrows == addm->nrows && thism->ncols == addm->ncols))
    die("mat_plus_eq: bad dimensions\n");
  for (i = 0; i < thism->nrows; i++)
    for (j = 0; j < thism->ncols; j++)
      thism->data[i][j] += addm->data[i][j];
}

void mat_minus_eq(Matrix *thism, Matrix *subm) {
  int i, j;
  if (!(thism->nrows == subm->nrows && thism->ncols == subm->ncols))
    die("ERROR mat_minus_eq: bad dimensions\n");
  for (i = 0; i < thism->nrows; i++)
    for (j = 0; j < thism->ncols; j++)
      thism->data[i][j] -= subm->data[i][j];
}

void mat_linear_comb(Matrix *dest, Matrix *src1, double coef1,
                     Matrix *src2, double coef2) {
  int i, j;
  if (!(dest->nrows == src1->nrows && dest->ncols == src1->ncols &&
	dest->nrows == src2->nrows && dest->ncols == src2->ncols))
    die("ERROR mat_linear_comb: bad dimensions\n");
  for (i = 0; i < dest->nrows; i++)
    for (j = 0; j < dest->ncols; j++)
      dest->data[i][j] = coef1*src1->data[i][j] + coef2*src2->data[i][j];
}

void mat_linear_comb_many(Matrix *dest, int n, Matrix **src, double *coef) {
  int i, j, k;
  if (n == 0) return;
  if (dest->nrows != src[0]->nrows || dest->ncols != src[0]->ncols)
    die("ERROR: mat_linear_comb_many: bad dimensions\n");
  for (i=1; i < n; i++)
    if (src[i]->nrows != src[0]->nrows ||
	src[i]->ncols != src[0]->ncols)
      die("ERROR: mat_linear_comb_many: bad dimensions in mat %i\n", i);
  for (i=0; i < dest->nrows; i++)  {
    for (j=0; j < dest->ncols; j++) {
      dest->data[i][j] = 0.0;
      for (k=0; k < n; k++)
	dest->data[i][j] += coef[k]*src[k]->data[i][j];
    }
  }
}


void mat_resize(Matrix *m, int nrows, int ncols) {
  int i;
  if (!(nrows >= 0 && ncols >= 0))
    die("ERROR mat_resize: nrows=%i ncols=%i\n", nrows, ncols);
  for (i = nrows; i < m->nrows; i++) sfree(m->data[i]);
  m->data = srealloc(m->data, nrows * sizeof(void*));
  for (i = 0; i < nrows; i++)
    m->data[i] = srealloc(m->data[i], ncols * sizeof(double));
  m->nrows = nrows;
  m->ncols = ncols;
}

/* Invert square, real, nonsymmetric matrix.  Uses LU decomposition
   (LAPACK routines dgetrf and dgetri).  Returns 0 on success, 1 on
   failure. */
int mat_invert(Matrix *M_inv, Matrix *M) {
#ifdef SKIP_LAPACK
  die("ERROR: LAPACK required for matrix inversion.\n");
#else
  int i, j;
  LAPACK_INT info, n = (LAPACK_INT)M->nrows, ipiv[n], lwork=(LAPACK_INT)n;
  LAPACK_DOUBLE tmp[n][n], work[lwork];

  if (!(M->nrows == M->ncols && M_inv->nrows == M_inv->ncols &&
	M->nrows == M_inv->nrows))
    die("ERROR mat_invert: bad dimensions\n");

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      tmp[i][j] = (LAPACK_DOUBLE)mat_get(M, j, i);

#ifdef R_LAPACK
  F77_CALL(dgetrf)(&n, &n, (LAPACK_DOUBLE*)tmp, &n, ipiv, &info);
#else
  dgetrf_(&n, &n, (LAPACK_DOUBLE*)tmp, &n, ipiv, &info);
#endif

  if (info != 0) {
    fprintf(stderr, "ERROR: unable to compute LU factorization of matrix (for matrix inversion); dgetrf returned value of %d.\n", (int)info);
    return 1;
  }
#ifdef R_LAPACK
  F77_CALL(dgetri)(&n, (LAPACK_DOUBLE*)tmp, &n, ipiv, work, &lwork, &info);
#else
  dgetri_(&n, (LAPACK_DOUBLE*)tmp, &n, ipiv, work, &lwork, &info);
#endif

  if (info != 0) {
    if (info > 0)
      fprintf(stderr, "ERROR: matrix is singular -- cannot invert.\n");
    else
      fprintf(stderr, "ERROR: unable to invert matrix.  Element %d had an illegal value (according to dgetri).\n", (int)info);
    return 1;
  }

  for (i = 0; i < M->nrows; i++)
    for (j = 0; j < M->nrows; j++)
      mat_set(M_inv, i, j, (double)tmp[j][i]);

#endif
  return 0;
}

/* Compute A = B * C * D where A, B, C, D are square matrices of the
   same dimension, and C is diagonal.  C is described by a vector
   representing its diagonal elements.  */
void mat_mult_diag(Matrix *A, Matrix *B, Vector *C, Matrix *D) {
  int i, j, k;
  for (i = 0; i < C->size; i++) {
    for (j = 0; j < C->size; j++) {
      A->data[i][j] = 0;
      for (k = 0; k < C->size; k++)
        A->data[i][j] += B->data[i][k] * C->data[k] * D->data[k][j];
    }
  }
}

int mat_equal(Matrix *A, Matrix *B) {
  int i, j;
  if (A->nrows != B->nrows || A->ncols != B->ncols) return 0;
  for (i=0; i < A->nrows; i++)
    for (j=0; j < A->ncols; j++)
      if (A->data[i][j] != B->data[i][j]) return 0;
  return 1;
}

#ifndef SKIP_LAPACK
void mat_to_lapack(Matrix *m, LAPACK_DOUBLE *arr) {
  int i, j, pos=0;
  for (j=0; j<m->ncols; j++)
    for (i=0; i< m->nrows; i++)
      arr[pos++] = (LAPACK_DOUBLE)m->data[i][j];
}

void mat_from_lapack(Matrix *m, LAPACK_DOUBLE *arr) {
  int i, j, pos=0;
  for (j=0; j < m->ncols; j++)
    for (i=0; i < m->nrows; i++)
      m->data[i][j] = (double)arr[pos++];
}

#endif
