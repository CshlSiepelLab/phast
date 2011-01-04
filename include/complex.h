/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file complex.h
    Functions and structures to hold, manipulate, and test complex numbers.
    All functions are inline; there is no complex.c file

    @ingroup base
*/

#ifndef COMPLEX_H
#define COMPLEX_H

#include <math.h>
#include <external_libs.h>

/** Structure representing complex number */
typedef struct {
  double x;			/**< real component */
  double y;			/**< imaginary component */
} Complex;

/** \name Initialize complex number function 
 \{ */

/** Create a complex number from Real and Imaginary components
    @param x Real component
    @param y Imaginary component
*/

static PHAST_INLINE
Complex z_set(double x, double y) {
  Complex z;
  z.x = x;
  z.y = y;
  return z;
}

/** \name Basic arithmetic for complex numbers functions
 \{ */

/** Add two complex numbers and return the result.
    @param z1 First complex number to add
    @param z2 Second complex number to add
    @result = z1 + z2
 */
static PHAST_INLINE
Complex z_add(Complex z1, Complex z2) {
  Complex sum;
  sum.x = z1.x + z2.x;
  sum.y = z1.y + z2.y;
  return sum;
}

/** Subtract one complex number from another and return the result.
   @param z1 Complex number to subtract from
   @param z2 Complex number to subtract
   @result  = z1 - z2
*/
static PHAST_INLINE
Complex z_sub(Complex z1, Complex z2) {
  Complex diff;
  diff.x = z1.x - z2.x;
  diff.y = z1.y - z2.y;
  return diff;
}

/** Multiply two complex numbers and return the result.
  @param z1 Complex number to multiply
  @param z2 Complex number to multiply by
  @result = z1 * z2
 */
static PHAST_INLINE
Complex z_mul(Complex z1, Complex z2) {
  Complex prod;
  prod.x = z1.x * z2.x - z1.y * z2.y;
  prod.y = z1.x * z2.y + z1.y * z2.x;
  return prod;
}

/** Divide one complex number by another and return the result.
  @param z1 Complex number to divide
  @param z2 Complex number to divide by
  @result = z1 / z2
*/
static PHAST_INLINE
Complex z_div(Complex z1, Complex z2) {
  Complex quot;
  double denom = z2.x * z2.x + z2.y * z2.y;
  quot.x = (z1.x * z2.x + z1.y * z2.y) / denom;
  quot.y = (z1.y * z2.x - z1.x * z2.y) / denom;
  return quot;
}

/** Multiply a complex number by a real number.
  @param z1 Complex number to multiply
  @param z2 Real number to multiply by
  @result = z1 * z2
 */
static PHAST_INLINE
Complex z_mul_real(Complex z, double x) {
  Complex prod;
  prod.x = z.x * x;
  prod.y = z.y * x;
  return prod;
}

/** Compute exponential function of complex number.
  @param z Power to raise to
  @return Exponential value of z
*/
static PHAST_INLINE
Complex z_exp(Complex z) {
  Complex retval;
  double exp_zx = exp(z.x);
  retval.x = exp_zx * cos(z.y);
  retval.y = exp_zx * sin(z.y);
  return retval;
}

/** Get absolute value of complex number.
  @param z Positive or negative complex number
  @result positive version of z
*/
static PHAST_INLINE
double z_abs(Complex z) {
  return sqrt(z.x * z.x + z.y * z.y);
}

/** \name Test equality complex number function
 \{ */


/** Check if two complex numbers are equal.
  @param z1 First complex number to compare
  @param z2 Second complex number to compare
  @result If (Z1 == Z2) { return 1; } else { return 0 }
*/
static PHAST_INLINE
int z_eq(Complex z1, Complex z2) {
  return(z1.x == z2.x && z1.y == z2.y);
}
/** \} */

#endif
