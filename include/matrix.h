/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: matrix.h,v 1.8 2008-11-16 02:32:54 acs Exp $ */

/* \file matrix.h
   Matrices of real numbers (doubles)
   \ingroup base
*/

#ifndef MAT_H
#define MAT_H

#include <vector.h>

/* Equality threshold -- consider equal if this close */
#define EQ_THRESHOLD 1e-10

/** Matrix structure -- just a 2d array of doubles and its dimensions */
struct matrix_struct {
  double **data;
  int nrows;
  int ncols;
};
typedef struct matrix_struct Matrix;

Matrix *mat_new(int nrows, int ncols);
Matrix *mat_new_from_array(double **array, int nrows, int ncols);
void mat_free(Matrix *m);
double mat_get(Matrix *m, int row, int col);
Vector *mat_get_row(Matrix *m, int row);
Vector *mat_get_col(Matrix *m, int col);
void mat_set(Matrix *m, int row, int col, double val);
void mat_set_identity(Matrix *m);
void mat_zero(Matrix *m);
void mat_set_all(Matrix *m, double val);
void mat_copy(Matrix *dest, Matrix *src);
Matrix *mat_create_copy(Matrix *src);
Matrix *mat_transpose(Matrix *src);
void mat_scale(Matrix *m, double scale_factor);
void mat_print(Matrix *m, FILE *F);
void mat_read(Matrix *m, FILE *F);
Matrix *mat_new_from_file(FILE *F, int nrows, int ncols);
void mat_mult(Matrix *prod, Matrix *m1, Matrix *m2);
void mat_vec_mult(Vector *prod, Matrix *m, Vector *v);
void mat_plus_eq(Matrix *thism, Matrix *addm);
void mat_minus_eq(Matrix *thism, Matrix *subm);
void mat_linear_comb(Matrix *dest, Matrix *src1, double coef1, 
                     Matrix *src2, double coef2);
void mat_resize(Matrix *m, int nrows, int ncols);
int mat_invert(Matrix *M_inv, Matrix *M);
void mat_mult_diag(Matrix *A, Matrix *B, Vector *C, Matrix *D);

/***************************************************************************
 * inline functions; also defined in matrix.c 
 ***************************************************************************/

/* we'll only inline the functions likely to be used heavily in inner
   loops */  

extern inline
double mat_get(Matrix *m, int row, int col) {
  return m->data[row][col];
}

extern inline
void mat_set(Matrix *m, int row, int col, double val) {
  m->data[row][col] = val;
}

#endif
