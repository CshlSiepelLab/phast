/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: complex_matrix.h,v 1.6 2008-11-16 02:32:54 acs Exp $ */

/** \file complex_matrix.h
    Matrices of complex numbers
    \ingroup base
*/

#ifndef ZMAT_H
#define ZMAT_H

#include <complex.h>
#include <complex_vector.h>
#include <external_libs.h>

/** Structure for matrix of complex numbers -- 2d array of Complex
    objects and its dimensions */
typedef struct {
  Complex **data;
  int nrows;
  int ncols;
} Zmatrix;

Zmatrix *zmat_new(int nrows, int ncols);
void zmat_free(Zmatrix *m);
Zvector *zmat_get_row(Zmatrix *m, int row);
Zvector *zmat_get_col(Zmatrix *m, int col);
void zmat_set_identity(Zmatrix *m);
void zmat_zero(Zmatrix *m);
void zmat_set_all(Zmatrix *m, Complex val);
void zmat_copy(Zmatrix *dest, Zmatrix *src);
Zmatrix *zmat_create_copy(Zmatrix *src);
void zmat_scale_complex(Zmatrix *m, Complex scale_factor);
void zmat_scale(Zmatrix *m, double scale_factor);
void zmat_print(Zmatrix *m, FILE *F);
void zmat_read(Zmatrix *m, FILE *F);
Zmatrix *zmat_new_from_file(FILE *F, int nrows, int ncols);
void zmat_mult(Zmatrix *prod, Zmatrix *m1, Zmatrix *m2);
void zmat_vec_mult(Zvector *prod, Zmatrix *m, Zvector *v);
void zmat_mult_real(Matrix *prod, Zmatrix *m1, Zmatrix *m2);
void zmat_plus_eq(Zmatrix *thism, Zmatrix *addm);
void zmat_minus_eq(Zmatrix *thism, Zmatrix *subm);
void zmat_mult_real_diag(Matrix *A, Zmatrix *B, Zvector *C, Zmatrix *D,
                         Zmatrix *scratch);
int zmat_as_real(Matrix *dest, Zmatrix *src, int strict);

/***************************************************************************
 * inline functions
 ***************************************************************************/

/* we'll only inline the functions likely to be used heavily in inner
   loops */  

static PHAST_INLINE
Complex zmat_get(Zmatrix *m, int row, int col) {
  return m->data[row][col];
}

static PHAST_INLINE
void zmat_set(Zmatrix *m, int row, int col, Complex val) {
  m->data[row][col] = val;
}

#endif
