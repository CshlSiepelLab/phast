/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file complex_vector.h
    Vectors of complex numbers
    @ingroup base
*/

#ifndef ZVEC_H
#define ZVEC_H

#include <complex.h>
#include <stdio.h>
#include <external_libs.h>

/** Structure for vector of complex numbers.
   Array of Complex objects and its length */
typedef struct { 
  Complex *data;		/**< array of Complex objects */ 
  int size;			/**< length of array */ 
} Zvector; 

/** \name Complex Vector allocation function
 \{ */

/** Create a new complex vector that supports 'size' elements 
  @param size Amount of elements vector can contain
  @result Newly allocated complex vector
*/
Zvector *zvec_new(int size);

/** \name Complex Vector cleanup function
 \{ */

/** Free complex vector v
   @param v Complex vector to free*/
void zvec_free(Zvector *v);

/** \name Complex Vector data access function
 \{ */

/* we'll only inline the functions likely to be used heavily in inner
   loops */  

/** Get an individual complex number of a vector at location i
   @param v Vector containing complex number to get
   @param i Index of complex number in Vector 'v' to get
   @result Complex number v[i]
 */
static PHAST_INLINE
Complex zvec_get(Zvector *v, int i) { /* check */
  return v->data[i];
}

/** \name Complex Vector data modification functions
 \{ */


/* we'll only inline the functions likely to be used heavily in inner
   loops */  

/** Set an individual complex number of a vector at location i
   @param v Vector containing complex number to set
   @param i Index of complex number in Vector to set
   @param val Complex value to set v[i] to */
static PHAST_INLINE
void zvec_set(Zvector *v, int i, Complex val) {
  v->data[i] = val;
}

/** Set all elements of complex vector to complex value val
   @param v Vector to modify all elements of
   @param val Value to set all elements of v to
 */
void zvec_set_all(Zvector *v, Complex val);

/** Zero a given vector
   @param v Vector to zero
 */
void zvec_zero(Zvector *v);


/** Vector add.
For each element in thisv add the corresponding element from addv 
@code
e.g.
[3,4,5,2] // thisv (initial)
+
[1,2,2,1] // addv
=
[4,6,7,3] //thisv (end)
@endcode
  @param[in,out] thisv vector to be modified
  @param[in] addv  vector containing complex numbers
*/
void zvec_plus_eq(Zvector *thisv, Zvector *addv);

/** Vector Subtract.
For each element in X in vector thisv, subtract complex number X in subv
@code
e.g.
[3,4,5,2] //thisv (initial)
-
[1,2,2,1] //subv
=
[2,2,3,1] //thisv (end)
@endcode
  @param[in,out] thisv vector to be modified
  @param[in] subv vector containing complex numbers  */
void zvec_minus_eq(Zvector *thisv, Zvector *subv);

/** Scale each element X in vector v by scale_factor.
 @code
e.g.
[2]        //scale_factor
*
[3,4,5,2]  //v (initial)
=
[6,8,10,4] //v(end)
@endcode
  @param[in,out] v vector to have its elements scaled
  @param scale_factor how much to multiply each element in v by
 */
void zvec_scale(Zvector *v, double scale_factor);

/** Compute Hadamard (pointwise) product of two vectors */
void zvec_had_prod(Zvector *dest, Zvector *src1, Zvector *src2);


/** \name Complex Vector data copy functions
 \{ */

/** Copy data from a complex vector into another existing complex vector
   @param dest Vector data is copied to
   @param src Vector data is copied from*/
void zvec_copy(Zvector *dest, Zvector *src);

/** Copy data from a complex vector into a new vector object
   @param Vector data is copied from
   @result Newly created vector with data from src
*/
Zvector* zvec_create_copy(Zvector *src);

/** \name Complex Vector read/save as file functions
 \{ */

/** Save vector data to a file
    @param v Vector to save to file
    @param F File descriptor of file to save to
*/
void zvec_print(Zvector *v, FILE *F);

/** Read vector data from a file into a complex vector
    @param v Vector to populate with data from file
    @param F File containing vector data to read in
*/
void zvec_read(Zvector *v, FILE *F);

/** Read vector data from a file into a newly allocated complex vector
    @param F File descriptor to file containing complex vector data
    @param size Size (Amount of elements) of vector to put file data into
    @result Newly allocated vector with data from file
*/
Zvector *zvec_new_from_file(FILE *F, int size);

/** \name Complex Vector cast function
 \{ */

/** Cast complex vector as real.
   Extract real component of each element. 
   @param[out] dest Vector of casted real numbers
   @param[in] src Vector of un-casted complex numbers
   @param[in] strict If == 1 and any of the real complex numbers contain an imaginary component that is not almost zero die
   @result If non-zero at least one "casted" complex number had an imaginary component that was not almost zero
*/
int zvec_as_real(Vector *dest, Zvector *src, int strict);

/** \} */



#endif
