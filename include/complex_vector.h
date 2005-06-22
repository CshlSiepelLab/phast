/* $Id: complex_vector.h,v 1.1 2005-06-22 07:11:20 acs Exp $
   Written by Adam Siepel, Summer 2005
   Copyright 2005, Adam Siepel, University of California
*/

/** \file complex_vector.h
    Vectors of complex numbers
    \ingroup base
*/

#ifndef ZVEC_H
#define ZVEC_H

#include <complex.h>
#include <stdio.h>

/** Structure for vector of complex numbers -- array of Complex
    objects and its length */
typedef struct { 
  Complex *data;		/**< array of Complex objects */ 
  int size;			/**< length of array */ 
} Zvector; 

Zvector *zvec_new(int size);
void zvec_free(Zvector *v);
Complex zvec_get(Zvector *v, int i);
void zvec_set(Zvector *v, int i, Complex val);
void zvec_set_all(Zvector *v, Complex val);
void zvec_copy(Zvector *dest, Zvector *src);
Zvector* zvec_create_copy(Zvector *src);
void zvec_print(Zvector *v, FILE *F);
void zvec_read(Zvector *v, FILE *F);
Zvector *zvec_new_from_file(FILE *F, int size);
void zvec_zero(Zvector *v);
void zvec_plus_eq(Zvector *thisv, Zvector *addv);
void zvec_minus_eq(Zvector *thisv, Zvector *subv);
void zvec_scale(Zvector *v, double scale_factor);

#endif
