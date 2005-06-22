/* $Id: complex.c,v 1.1 2005-06-22 07:11:20 acs Exp $
   Written by Adam Siepel, Summer 2005
   Copyright 2005, Adam Siepel, University of California
*/

/** \file complex_matrix.c
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
