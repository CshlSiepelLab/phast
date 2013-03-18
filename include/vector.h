/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** \file vector.h
    Functions and structures to hold, manipulate, and normalize a vector of type double.
    \ingroup base
*/

#ifndef VEC_H
#define VEC_H

#include <stdio.h>
#include <external_libs.h>

struct lst_struct;

/* temporary -- elim gsl */

struct matrix_struct;

/** Vector structure -- just an array of doubles and its length */
typedef struct {
  double *data;			/**< underlying array of doubles */
  int size;			/**< length of array */
} Vector;

/** \name Vector create/free operations. */

/** \{ */

/** Create new vector.
  @param size Number of elements.
*/
Vector *vec_new(int size);

/** Create new vector and initialized it with values from given array.
  
  @param array Initialization values.
  @param size Number of elements.
*/
Vector *vec_new_from_array(double *array, int size);

/** Create new vector and initialized it with values from given list.
  
  Vector size is the size of the supplied list.
  
  @param l Initialization values.
  
  \sa lists.h
*/
Vector *vec_new_from_list(struct lst_struct *l);

/** Allocate a new vector a initialize it from file.

  @param F File to use as input (must have at least as many values as the given
           size).
  @param size Size of the new vector.
  
  \sa vec_read
*/
Vector *vec_new_from_file(FILE *F, int size);

/** Allocate a new vector with identical contents as the given vector.

  @param src Vector to be copied.
*/
Vector* vec_create_copy(Vector *src);

/** Change the size of a vector

    @param v vector to be resized
    @param new_size New vector length
    @return v.  v->data has been reallocated and v->size updated.
    @note This behaves like regular realloc, in that it will not affect elements 0..min(old_size, new_size)-1.  If new_size > old_size, elements with indices [old_size, ..., new_size-1] will not be initialized.
*/
Vector *vec_realloc(Vector *v, int new_size);

/** Release memory used by vector.
  
  @param v Vector to be freed.
*/
void vec_free(Vector *v);

/** \} */

/** \name Vector initialization operations. */

/** \{ */

/** Set all elements to zero.

  @param v Vector to be reset to zero.
  
  \sa vec_set_all
*/
void vec_zero(Vector *v);

/** Set all elements of vector to a given value.

  @param v Vector of elements.
  @param val New value for vector elements.
*/
void vec_set_all(Vector *v, double val);

/** Read a vector from file.

  @param v Vector where values are stored.
  @param F File to use as input (must have at least as many values as the size
           of the supplied vector).

  \sa vec_new_from_file
*/
void vec_read(Vector *v, FILE *F);

/** Copy content between two vectors of equal size.

  @param src Source vector.
  @param dest Destination vector.
*/
void vec_copy(Vector *dest, Vector *src);

/** \} */

/** Retrieve value of i'th element from vector.

  @param v Vector of elements.
  @param i 0-based element index.
*/
static PHAST_INLINE
double vec_get(Vector *v, int i) {
  return v->data[i];

}

/** Set value of i'th element in vector.

  @param v Vector of elements.
  @param i 0-based element index.
  @param val New value to store.
*/
static PHAST_INLINE
void vec_set(Vector *v, int i, double val) {
  v->data[i] = val;
}

/** Print supplied vector.

  Prints the vector, with values separated by spaces, terminated with a newline.
  
  @param v Vector to be printed.
  @param F File to use as output.
  
  \sa vec_fprintf
*/
void vec_print(Vector *v, FILE *F);

/** Print supplied vector.

  Prints the vector, with values separated by spaces, terminated with a newline.
  Output is formatted using the supplied formatstr which should contain the format
  description for a single value.
  
  @param v Vector to be printed.
  @param F File to use as output.
  @param formatstr Element format description as in fprintf.
  
  \sa vec_print
*/
void vec_fprintf(Vector *v, FILE *F, char *formatstr);

/** \name Vector math operations. */

/** \{ */

/** Increment vector elements with values of supplied vector.

  \code 
  // equivalent to this, but with vectors
  thisv += addv; 
  \endcode
  
  @param thisv Vector to be incremented.
  @param addv Vector to be used as increment.
  
  \sa vec_minus_eq, vec_scale
*/
void vec_plus_eq(Vector *thisv, Vector *addv);

/** Decrement vector elements with values of supplied vector.

  \code 
  // equivalent to this, but with vectors
  thisv -= subv; 
  \endcode
  
  @param thisv Vector to be decremented.
  @param subv Vector to be used as decrement.
  
  \sa vec_plus_eq, vec_scale
*/
void vec_minus_eq(Vector *thisv, Vector *subv);

/** Scale vector elements by supplied constant.

  \code 
  // equivalent to this, but with vectors
  v *= scale_factor;
  \endcode
  
  @param v Vector to be scaled.
  @param scale_factor Scale factor.
  
  \sa vec_plus_eq, vec_minus_eq
*/
void vec_scale(Vector *v, double scale_factor);

/** Compute inner product of two n-dimensional real-valued vectors.
  @param v1 First input vector (n-dim).
  @param v2 Second input vector (n-dim).
  
  \sa vec_outer_prod, vec_norm
*/
double vec_inner_prod(Vector *v1, Vector *v2); 

/** Compute outer (cross) product of two n-dimensional real-valued vectors.

  @param mat Computed cross-product matrix (preallocate n x n)
  @param v1 First input vector (n-dim).
  @param v2 Second input vector (n-dim).
  
  \sa vec_inner_prod, vec_norm
  \sa matrix.h
*/
void vec_outer_prod(struct matrix_struct *mat, Vector *v1, Vector *v2); 

/** Compute a 2-norm of vector.

  @param v Input vector.
  \sa vec_inner_prod, vec_outer_prod
*/
double vec_norm(Vector *v);

/** Force elements of vector to sum to 1.

  @param v Input vector.
*/
void vec_normalize(Vector *v);

/** Compute point wise average of vectors. 

  If counts is NULL, each source vector is assumed to have a count of 1.
  
  @param dest_v n-dim output vector.
  @param source_vs List of n-dim source vector pointers.
  @param counts List of (integer) vector counts. If NULL, count = 1.

  \sa lists.h
*/
void vec_ave(Vector *dest_v, struct lst_struct *source_vs, 
	     struct lst_struct *counts);

/* \} */

#endif
