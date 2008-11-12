/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: complex.c,v 1.3 2008-11-12 02:07:59 acs Exp $ */

/** \file complex.c
    Complex numbers
    \ingroup base
*/

#include <complex.h>
#include <math.h>

Complex z_set(double x, double y) {
  Complex z;
  z.x = x;
  z.y = y;
  return z;
}

Complex z_add(Complex z1, Complex z2) {
  Complex sum;
  sum.x = z1.x + z2.x;
  sum.y = z1.y + z2.y;
  return sum;
}

Complex z_sub(Complex z1, Complex z2) {
  Complex diff;
  diff.x = z1.x - z2.x;
  diff.y = z1.y - z2.y;
  return diff;
}

Complex z_mul(Complex z1, Complex z2) {
  Complex prod;
  prod.x = z1.x * z2.x - z1.y * z2.y;
  prod.y = z1.x * z2.y + z1.y * z2.x;
  return prod;
}

Complex z_div(Complex z1, Complex z2) {
  Complex quot;
  double denom = z2.x * z2.x + z2.y * z2.y;
  quot.x = (z1.x * z2.x + z1.y * z2.y) / denom;
  quot.y = (z1.y * z2.x - z1.x * z2.y) / denom;
  return quot;
}

Complex z_mul_real(Complex z, double x) {
  Complex prod;
  prod.x = z.x * x;
  prod.y = z.y * x;
  return prod;
}

Complex z_exp(Complex z) {
  Complex retval;
  double exp_zx = exp(z.x);
  retval.x = exp_zx * cos(z.y);
  retval.y = exp_zx * sin(z.y);
  return retval;
}

double z_abs(Complex z) {
  return sqrt(z.x * z.x + z.y * z.y);
}

int z_eq(Complex z1, Complex z2) {
  return(z1.x == z2.x && z1.y == z2.y);
}
