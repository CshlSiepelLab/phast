/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: vector.c,v 1.3 2008-11-12 02:07:59 acs Exp $*/

#include <vector.h>
#include <misc.h>
#include <lists.h>

Vector *vec_new(int size) {
  Vector *v = smalloc(sizeof(Vector));
  v->data = smalloc(size * sizeof(double));
  v->size = size;
  return v;
}

Vector *vec_new_from_array(double *array, int size) {
  int i;
  Vector *v = vec_new(size);
  for (i = 0; i < size; i++)
    v->data[i] = array[i];
  return v;
}

Vector *vec_new_from_list(List *l) {
  int i;
  Vector *v = vec_new(lst_size(l));
  for (i = 0; i < v->size; i++)
    v->data[i] = lst_get_dbl(l, i);
  return v;
}

void vec_free(Vector *v) {
  free(v->data);
  free(v);
}

void vec_set_all(Vector *v, double val) {
  int i;
  for (i = 0; i < v->size; i++)
    v->data[i] = val;
}

void vec_copy(Vector *dest, Vector *src) {
  int i;
  if (dest->size != src->size)
    die("ERROR vec_copy bad dimensions\n");
  for (i = 0; i < src->size; i++) 
    dest->data[i] = src->data[i];
}

Vector* vec_create_copy(Vector *src) {
  Vector *copy = vec_new(src->size);
  vec_copy(copy, src);
  return copy;
}

void vec_print(Vector *v, FILE *F) {
  int i;
  for (i = 0; i < v->size; i++)
    fprintf(F, "%f ", v->data[i]);
  fprintf(F, "\n");
}

void vec_fprintf(Vector *v, FILE *F, char *formatstr) {
  int i;
  for (i = 0; i < v->size; i++)
    fprintf(F, formatstr, v->data[i]);
  fprintf(F, "\n");
}

void vec_read(Vector *v, FILE *F) {
  int i;
  for (i = 0; i < v->size; i++)
    fscanf(F, "%lf ", &v->data[i]);
}

Vector *vec_new_from_file(FILE *F, int size) {
  Vector *v = vec_new(size);
  vec_read(v, F);
  return v;
}

void vec_zero(Vector *v) {
  int i;
  for (i = 0; i < v->size; i++) v->data[i] = 0;
}

void vec_plus_eq(Vector *thisv, Vector *addv) {
  int i;
  if (thisv->size != addv->size)
    die("ERROR vec_plus_eq: bad dimensions\n");
  for (i = 0; i < thisv->size; i++) 
    thisv->data[i] += addv->data[i];  
}

void vec_minus_eq(Vector *thisv, Vector *subv) {
  int i;
  if (thisv->size != subv->size)
    die("ERROR vec_minus_eq: bad dimensions\n");
  for (i = 0; i < thisv->size; i++) 
    thisv->data[i] -= subv->data[i];  
}

void vec_scale(Vector *v, double scale_factor) {
  int i;
  for (i = 0; i < v->size; i++) 
    v->data[i] *= scale_factor;
}

/* Compute and return inner (dot) product of two n-dimensional
   real-valued vectors. */
double vec_inner_prod(Vector *v1, /* first input vector (n-dim) */
		      Vector *v2 /* second input vector (n-dim) */
		      ) {
  int i;
  double retval = 0;
  if (v1->size != v2->size)
    die("ERROR vec_inner_prod: bad dimensions\n");
  for (i = 0; i < v1->size; i++)
    retval += (vec_get(v1, i) * vec_get(v2, i));
  return retval;
}

/* Compute and return outer (cross) product of two n-dimensional
   real-valued vectors. */
void vec_outer_prod(Matrix *mat, /* computed cross-product matrix
                                         (preallocate n x n) */
		    Vector *v1, /* first input vector (n-dim) */
		    Vector *v2 /* second input vector (n-dim) */
		    ) {
  int i, j;
  if (!(v1->size == v2->size && mat->nrows == v1->size && 
	mat->ncols == v1->size))
    die("ERROR vec_outer_prod: bad dimensions\n");
  for (i = 0; i < v1->size; i++) {
    for (j = 0; j < v1->size; j++) {
      mat_set(mat, i, j, vec_get(v1, i) * vec_get(v2, j));
    }
  }
}

/* Compute and return 2-norm of vector. */
double vec_norm(Vector *v) {
  double ss = 0;
  int i;
  for (i = 0; i < v->size; i++) 
    ss += v->data[i] * v->data[i];
  return sqrt(ss);
}


/* Force elements of vector to sum to 1 */
void vec_normalize(Vector *v) {
  double sum=0.0;
  int i;
  for (i=0; i<v->size; i++) 
    sum += v->data[i];
  for (i=0; i<v->size; i++)
    v->data[i] /= sum;
}


/* Compute pointwise average of vectors.  If counts is NULL, each
   source vector is assumed to have a count of 1 */
void vec_ave(Vector *dest_v, List *source_vs, List *counts) {
  int i, j;
  double n = 0;

  if (lst_size(source_vs) < 1)
    die("ERROR vec_ave: lst_size(source_vs)=%i\n", lst_size(source_vs));

  if (counts != NULL) {
    if (lst_size(source_vs) != lst_size(counts))
      die("ERROR vec_ave: lst_size(source_vs)=%i != lst_size(counts)=%i\n",
	  lst_size(source_vs), lst_size(counts));
    for (i = 0; i < lst_size(counts); i++) n += lst_get_int(counts, i);
  }
  else n = lst_size(source_vs);

  vec_set_all(dest_v, 0);
  
  for (i = 0; i < lst_size(source_vs); i++) {
    int count = (counts == NULL ? 1 : lst_get_int(counts, i));
    for (j = 0; j < dest_v->size; j++) 
      dest_v->data[j] += vec_get(lst_get_ptr(source_vs, i), j) * count;
  }

  vec_scale(dest_v, 1/n);
}

