/* $Id: vector.h,v 1.1 2005-06-22 07:11:20 acs Exp $
   Written by Adam Siepel, Summer 2005
   Copyright 2005, Adam Siepel, University of California
*/

/** \file vector.h
    Vectors of doubles.  Very simple implementation -- essentially
    just an array.
    \ingroup base
*/

#ifndef VEC_H
#define VEC_H

#include <stdio.h>
struct lst_struct;

/* temporary -- elim gsl */

struct matrix_struct;

/** Vector structure -- just an array of doubles and its length */
typedef struct {
  double *data;			/**< underlying array of doubles */
  int size;			/**< length of array */
} Vector;

Vector *vec_new(int size);
Vector *vec_new_from_array(double *array, int size);
Vector *vec_new_from_list(struct lst_struct *l);
void vec_free(Vector *v);
double vec_get(Vector *v, int i);
void vec_set(Vector *v, int i, double val);
void vec_set_all(Vector *v, double val);
void vec_copy(Vector *dest, Vector *src);
Vector* vec_create_copy(Vector *src);
void vec_print(Vector *v, FILE *F);
void vec_fprintf(Vector *v, FILE *F, char *formatstr);
void vec_read(Vector *v, FILE *F);
Vector *vec_new_from_file(FILE *F, int size);
void vec_zero(Vector *v);
void vec_plus_eq(Vector *thisv, Vector *addv);
void vec_minus_eq(Vector *thisv, Vector *subv);
void vec_scale(Vector *v, double scale_factor);
double vec_inner_prod(Vector *v1, Vector *v2); 
void vec_outer_prod(struct matrix_struct *mat, Vector *v1, Vector *v2); 
double vec_norm(Vector *v);

#endif
