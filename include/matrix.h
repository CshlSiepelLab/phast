/* $Id: matrix.h,v 1.3 2005-06-24 17:41:52 acs Exp $ 
   Written by Adam Siepel, 2002-2005
   Copyright 2002-2005, Adam Siepel, University of California 
*/

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
void mat_scale(Matrix *m, double scale_factor);
void mat_print(Matrix *m, FILE *F);
void mat_read(Matrix *m, FILE *F);
Matrix *mat_new_from_file(FILE *F, int nrows, int ncols);
void mat_mult(Matrix *prod, Matrix *m1, Matrix *m2);
void mat_vec_mult(Vector *prod, Matrix *m, Vector *v);
void mat_plus_eq(Matrix *thism, Matrix *addm);
void mat_minus_eq(Matrix *thism, Matrix *subm);
int mat_invert(Matrix *M_inv, Matrix *M);

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
