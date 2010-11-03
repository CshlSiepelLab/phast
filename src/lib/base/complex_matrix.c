/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: complex_matrix.c,v 1.4 2008-11-16 02:32:54 acs Exp $*/

/** \file complex_matrix.c
    Matrices of complex numbers
    \ingroup base
*/

#include <misc.h>
#include <complex_matrix.h>

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
    sfree(m->data[i]);
  sfree(m->data);
  sfree(m);
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
  if (!(dest->nrows == src->nrows && dest->ncols == src->ncols))
    die("ERROR zmat_copy: bad dimensions\n");
  for (i = 0; i < dest->nrows; i++)
    for (j = 0; j < dest->ncols; j++)
      dest->data[i][j] = src->data[i][j];
}

Zmatrix *zmat_create_copy(Zmatrix *src) {
  Zmatrix *dest = zmat_new(src->nrows, src->ncols);
  zmat_copy(dest, src);
  return dest;
}

void zmat_scale_complex(Zmatrix *m, Complex val) {
   int i,j;
   for (i=0; i<m->nrows; i++)
     for (j=0; j<m->ncols; j++)
	m->data[i][j] = z_mul(m->data[i][j], val);
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
      if (2 != fscanf(F, "%lf %lf ", &m->data[i][j].x, &m->data[i][j].y))
	die("ERROR reading complex matrix");
}

Zmatrix *zmat_new_from_file(FILE *F, int nrows, int ncols) {
  Zmatrix *m = zmat_new(nrows, ncols);
  zmat_read(m, F);
  return m;
}

void zmat_mult(Zmatrix *prod, Zmatrix *m1, Zmatrix *m2) {
  int i, j, k;
  if (!(m1->ncols == m2->nrows && m1->nrows == m2->ncols && 
	prod->nrows == m1->nrows && prod->ncols == m2->ncols))
    die("ERROR zmat_mult: bad dimensions\n");
  zmat_zero(prod);
  for (i = 0; i < prod->nrows; i++) 
    for (j = 0; j < prod->ncols; j++) 
      for (k = 0; k < m1->ncols; k++) 
	prod->data[i][j] = z_add(prod->data[i][j],
				 z_mul(m1->data[i][k], m2->data[k][j]));
}

void zmat_vec_mult(Zvector *prod, Zmatrix *m, Zvector *v) {
  int i, j;
  if (!(m->nrows == v->size && v->size == prod->size))
    die("ERROR zmat_vec_mult: bad dimensions\n");
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
  int i, j, k;
  if (!(m1->ncols == m2->nrows && m1->nrows == m2->ncols && 
	prod->nrows == m1->nrows && prod->ncols == m2->ncols))
    die("ERROR zmat_mult_real: bad dimensions\n");
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
  if (!(thism->nrows == addm->nrows && thism->ncols == addm->ncols))
    die("ERROR zmat_plus_eq: bad dimensions\n");
  for (i = 0; i < thism->nrows; i++)
    for (j = 0; j < thism->ncols; j++)  
      thism->data[i][j] = z_add(thism->data[i][j], addm->data[i][j]);
}

void zmat_minus_eq(Zmatrix *thism, Zmatrix *subm) {
  int i, j;
  if (!(thism->nrows == subm->nrows && thism->ncols == subm->ncols))
    die("ERROR zmat_minus_eq: bad dimensions\n");
  for (i = 0; i < thism->nrows; i++)
    for (j = 0; j < thism->ncols; j++)  
      thism->data[i][j] = z_sub(thism->data[i][j], subm->data[i][j]);
}

/* Compute A = B * C * D where A, B, C, D are square matrices of the
   same dimension, and C is diagonal.  Allow B, C, D to be complex
   valued but assume their product is real valued (as when B,C,D
   represent diagonalization of A).  C is described by a vector
   representing its diagonal elements.  A temp matrix can optionally
   be passed in to improve efficiency.  */
void zmat_mult_real_diag(Matrix *A, Zmatrix *B, Zvector *C, Zmatrix *D,
                         Zmatrix *scratch) {
  int i, j;
  int size = C->size;
  Zmatrix *tmp;

  if (!(A->nrows == A->ncols && A->nrows == B->nrows &&
	B->nrows == B->ncols && B->nrows == C->size && 
	C->size == D->nrows && D->nrows == D->ncols))
    die("ERROR zmat_mult_real_diag: bad dimensions\n");

  if (scratch == NULL) 
    tmp= zmat_new(size, size);
  else {
    if (!(scratch->nrows == size && scratch->ncols == size))
      die("ERROR zmat_mult_real_diag: scratch has wrong size\n");
    tmp = scratch;
  }

  /* first compute tmp = C*D */
  for (i = 0; i < size; i++) 
    for (j = 0; j < size; j++) 
      zmat_set(tmp, i, j, z_mul(zvec_get(C, i), zmat_get(D, i, j))); 

  /* now compute A = B*tmp */
  zmat_mult_real(A, B, tmp);

  if (scratch == NULL)
    zmat_free(tmp);
}

/* "cast" complex matrix as real, by extracting real component of each
   element.  If strict == TRUE ensure imaginary components are zero
   (or very close)  */
int zmat_as_real(Matrix *dest, Zmatrix *src, int strict) {
  int i, j, rv=0;
  if (!(dest->nrows == src->nrows && dest->ncols == src->ncols))
    die("ERROR zmat_as_real: bad dimensions\n");
  for (i = 0; i < src->nrows; i++) {
    for (j = 0; j < src->ncols; j++) {
      dest->data[i][j] = src->data[i][j].x;
      if (src->data[i][j].y >= 1e-6) {
	rv = 1;
	if (strict) 
	  die("ERROR in zmat_as_real: src matrix has imaginary component %ei",
	      src->data[i][j].y);
      }
    }
  }
  return rv;
}
