/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* \file complex_vector.c
    Vectors of complex numbers
    \ingroup base
*/

#include <misc.h>
#include <complex_vector.h> 


Zvector *zvec_new(int size) {
  Zvector *v = smalloc(sizeof(Zvector));
  v->data = smalloc(size * sizeof(Complex));
  v->size = size;
  return v;
}

void zvec_free(Zvector *v) {
  sfree(v->data);
  sfree(v);
}

void zvec_set_all(Zvector *v, Complex val) {
  int i;
  for (i = 0; i < v->size; i++)
    v->data[i] = val;
}

void zvec_copy(Zvector *dest, Zvector *src) {
  int i;
  if (dest->size != src->size)
    die("ERROR zvec_copy bad dimensions\n");
  for (i = 0; i < src->size; i++) 
    dest->data[i] = src->data[i];
}

Zvector* zvec_create_copy(Zvector *src) {
  Zvector *copy = zvec_new(src->size);
  zvec_copy(copy, src);
  return copy;
}

void zvec_print(Zvector *v, FILE *F) {
  int i;
  for (i = 0; i < v->size; i++)
    fprintf(F, "%f %f ", v->data[i].x, v->data[i].y);
  fprintf(F, "\n");
}

void zvec_read(Zvector *v, FILE *F) {
  int i;
  for (i = 0; i < v->size; i++)
    if (2 != fscanf(F, "%lf %lf ", &(v->data[i].x), &(v->data[i].y)))
      die("ERROR reading complex vector");
}

Zvector *zvec_new_from_file(FILE *F, int size) {
  Zvector *v = zvec_new(size);
  zvec_read(v, F);
  return v;
}

void zvec_zero(Zvector *v) {
  int i;
  for (i = 0; i < v->size; i++) 
    v->data[i] = z_set(0, 0);
}

void zvec_plus_eq(Zvector *thisv, Zvector *addv) {
  int i;
  if (thisv->size != addv->size)
    die("zvec_plus_eq: bad dimensions\n");
  for (i = 0; i < thisv->size; i++) 
    thisv->data[i] = z_add(thisv->data[i], addv->data[i]);  
}

void zvec_minus_eq(Zvector *thisv, Zvector *subv) {
  int i;
  if (thisv->size != subv->size)
    die("ERROR zvec_minus_eq: bad dimensions\n");
  for (i = 0; i < thisv->size; i++) 
    thisv->data[i] = z_sub(thisv->data[i], subv->data[i]);  
}

void zvec_scale(Zvector *v, double scale_factor) {
  int i;
  for (i = 0; i < v->size; i++) 
    v->data[i] = z_mul_real(v->data[i], scale_factor);
}

/* compute Hadamard (pointwise) product of two vectors */
void zvec_had_prod(Zvector *dest, Zvector *src1, Zvector *src2) {
  int i;
  if (!(dest->size == src1->size && src1->size == src2->size))
    die("ERROR zvec_had_prod: bad dimensions\n");
  for (i = 0; i < dest->size; i++)
    dest->data[i] = z_mul(src1->data[i], src2->data[i]);
}

/* "cast" complex vector as real, by extracting real component of each
   element.  If strict == TRUE ensure imaginary components are zero
   (or very close)  */
int zvec_as_real(Vector *dest, Zvector *src, int strict) {
  int i, rv = 0;
  if (dest->size != src->size)
    die("ERROR zvec_as_real: bad dimensions\n");
  for (i = 0; i < src->size; i++) {
    dest->data[i] = src->data[i].x;
    if (src->data[i].y >= 1.0e-6) {
      rv = 1;
      if (strict) 
	die("Error in zvec_as_real: src vector has imaginary component %ei", src->data[i].y);
    }
  }
  return rv;
}
