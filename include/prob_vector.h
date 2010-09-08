/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#ifndef PROB_VECTOR
#define PROB_VECTOR

#include <vector.h>

typedef enum {LOWER, UPPER, TWOTAIL} p_val_type;

void pv_stats(Vector *p, double *mean, double *var);
void pv_confidence_interval(Vector *p, double size, int *interval_min, 
                            int *interval_max);
int* pv_quantiles(Vector *p);
double pv_p_value(Vector *distrib, double x_0, p_val_type side);
void pv_p_values(Vector *distrib, double *x_0, int n, double *pvals,
                 p_val_type side);
void pv_normalize(Vector *p);
Vector *pv_convolve(Vector *p, int n, double epsilon);
Vector **pv_convolve_save(Vector *p, int n, double epsilon);
Vector *pv_convolve_many(Vector **p, int *counts, int n, double epsilon);
Vector *pv_poisson(double lambda, double epsilon);
Vector *pv_convolve_fast(Vector *p, int n, double epsilon);
Vector *pv_cdf(Vector *pdf, p_val_type side);
int pv_draw_idx_arr(double *arr, int n);
int pv_draw_idx(Vector *pdf);

#endif
