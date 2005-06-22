/* $Id: complex_matrix.c,v 1.1 2005-06-22 07:11:20 acs Exp $
   Written by Adam Siepel, Summer 2005
   Copyright 2005, Adam Siepel, University of California
*/

/** \file complex_matrix.c
    Matrices of complex numbers
    \ingroup base
*/

#include <misc.h>
#include <complex_matrix.h>
#include <assert.h>

Zmatrix *zmat_new(int nrows, int ncols) {
  int i;
  Zmatrix *m = smalloc(sizeof(Zmatrix));
  m->data = smalloc(nrows * sizeof(void*));
  for (i = 0; i < nrows; i++)
    m->data[i] = smalloc(ncols * sizeof(Complex));
  m->nrows = nrows;
  m->ncols = ncols;
  return m;
}

void zmat_free(Zmatrix *m) {
  int i;
  for (i = 0; i < m->nrows; i++)
    free(m->data[i]);
  free(m->data);
  free(m);
}

Complex zmat_get(Zmatrix *m, int row, int col) {
  return m->data[row][col];
}

Zvector *zmat_get_row(Zmatrix *m, int row) {
  int j;
  Zvector *v = zvec_new(m->ncols);
  for (j = 0; j < m->ncols; j++)
    v->data[j] = m->data[row][j];
  return v;
}

Zvector *zmat_get_col(Zmatrix *m, int col) {
  int i;
  Zvector *v = zvec_new(m->nrows);
  for (i = 0; i < m->nrows; i++)
    v->data[i] = m->data[i][col];
  return v;
}

void zmat_set(Zmatrix *m, int row, int col, Complex val) {
  m->data[row][col] = val;
}

void zmat_set_identity(Zmatrix *m) {
  int i, j;
  for (i = 0; i < m->nrows; i++)
    for (j = 0; j < m->ncols; j++)
      m->data[i][j] = z_set(i == j ? 1 : 0, 0);
}

void zmat_zero(Zmatrix *m) {
  int i, j;
  for (i = 0; i < m->nrows; i++)
    for (j = 0; j < m->ncols; j++)
      m->data[i][j] = z_set(0, 0);
}

void zmat_set_all(Zmatrix *m, Complex val) {
  int i, j;
  for (i = 0; i < m->nrows; i++)
    for (j = 0; j < m->ncols; j++)
      m->data[i][j] = val;
}

void zmat_copy(Zmatrix *dest, Zmatrix *src) {
  int i, j;
  assert(dest->nrows == src->nrows && dest->ncols == src->ncols);
  for (i = 0; i < dest->nrows; i++)
    for (j = 0; j < dest->ncols; j++)
      dest->data[i][j] = src->data[i][j];
}

Zmatrix *zmat_create_copy(Zmatrix *src) {
  Zmatrix *dest = zmat_new(src->nrows, src->ncols);
  zmat_copy(dest, src);
  return dest;
}

void zmat_scale(Zmatrix *m, double scale_factor) {
  int i, j;
  for (i = 0; i < m->nrows; i++)
    for (j = 0; j < m->ncols; j++)
      m->data[i][j] = z_mul_real(m->data[i][j], scale_factor);
}

void zmat_print(Zmatrix *m, FILE *F) {
  int i, j;
  for (i = 0; i < m->nrows; i++) {
    for (j = 0; j < m->ncols; j++) 
      fprintf(F, "%14.6e %14.6e ", m->data[i][j].x, m->data[i][j].y);
    fprintf(F, "\n");
  }
}

void zmat_read(Zmatrix *m, FILE *F) {
  int i, j;
  for (i = 0; i < m->nrows; i++)
    for (j = 0; j < m->ncols; j++)
      fscanf(F, "%lf %lf ", &m->data[i][j].x, &m->data[i][j].y);
}

Zmatrix *zmat_new_from_file(FILE *F, int nrows, int ncols) {
  Zmatrix *m = zmat_new(nrows, ncols);
  zmat_read(m, F);
  return m;
}

void zmat_mult(Zmatrix *prod, Zmatrix *m1, Zmatrix *m2) {
  assert(m1->ncols == m2->nrows && m1->nrows == m2->ncols && 
         prod->nrows == m1->nrows && prod->ncols == m2->ncols);
  int i, j, k;
  zmat_zero(prod);
  for (i = 0; i < prod->nrows; i++) 
    for (j = 0; j < prod->ncols; j++) 
      for (k = 0; k < m1->ncols; k++) 
	prod->data[i][j] = z_add(prod->data[i][j],
				 z_mul(m1->data[i][k], m2->data[k][j]));
}

void zmat_vec_mult(Zvector *prod, Zmatrix *m, Zvector *v) {
  int i, j;
  assert(m->nrows == v->size && v->size == prod->size);
  for (i = 0; i < m->nrows; i++) {
    prod->data[i] = z_set(0, 0);
    for (j = 0; j < m->ncols; j++) {
      prod->data[i] = z_add(prod->data[i],
			    z_mul(m->data[i][j], v->data[j]));
    }
  }
}

/* multiply two complex matrices whose product is expected to be real */
void zmat_mult_real(Matrix *prod, Zmatrix *m1, Zmatrix *m2) {
  assert(m1->ncols == m2->nrows && m1->nrows == m2->ncols && 
         prod->nrows == m1->nrows && prod->ncols == m2->ncols);
  int i, j, k;
  mat_zero(prod);
  for (i = 0; i < prod->nrows; i++) {
    for (j = 0; j < prod->ncols; j++) {
      Complex z = z_set(0, 0);
      for (k = 0; k < m1->ncols; k++) 
	z = z_add(z, z_mul(m1->data[i][k], m2->data[k][j]));
      if (z.y > 1e-6)
	die("ERROR in zmat_mult_real: product of complex matrices not real.\n");
      mat_set(prod, i, j, z.x);
    }
  }
}

void zmat_plus_eq(Zmatrix *thism, Zmatrix *addm) {
  int i, j;
  assert(thism->nrows == addm->nrows && thism->ncols == addm->ncols);
  for (i = 0; i < thism->nrows; i++)
    for (j = 0; j < thism->ncols; j++)  
      thism->data[i][j] = z_add(thism->data[i][j], addm->data[i][j]);
}

void zmat_minus_eq(Zmatrix *thism, Zmatrix *subm) {
  int i, j;
  assert(thism->nrows == subm->nrows && thism->ncols == subm->ncols);
  for (i = 0; i < thism->nrows; i++)
    for (j = 0; j < thism->ncols; j++)  
      thism->data[i][j] = z_sub(thism->data[i][j], subm->data[i][j]);
}
