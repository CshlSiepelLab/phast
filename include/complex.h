/* $Id: complex.h,v 1.1 2005-06-22 07:11:20 acs Exp $
   Written by Adam Siepel, Summer 2005
   Copyright 2005, Adam Siepel, University of California
*/

/** \file complex_matrix.h
    Complex numbers
    \ingroup base
*/

#ifndef COMPLEX_H
#define COMPLEX_H

/** Structure representing complex number */
typedef struct {
  double x;			/* real component */
  double y;			/* imaginary component */
} Complex;

Complex z_set(double x, double y);
Complex z_add(Complex z1, Complex z2);
Complex z_sub(Complex z1, Complex z2);
Complex z_mul(Complex z1, Complex z2);
Complex z_div(Complex z1, Complex z2);
Complex z_mul_real(Complex z, double x);
Complex z_exp(Complex z);
double z_abs(Complex z);
int z_eq(Complex z1, Complex z2);

#endif
