/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#ifndef PROB_MATRIX
#define PROB_MATRIX

#include <vector.h>
#include <matrix.h>

#define PM_EPS 1e-10

void pm_mean(Matrix *p, double *mean_x, double *mean_y);
Vector *pm_marg_x(Matrix *p);
Vector *pm_marg_y(Matrix *p);
Vector *pm_marg_tot(Matrix *p);
Vector *pm_x_given_tot(Matrix *p, int tot);
Vector *pm_y_given_tot(Matrix *p, int tot);
void pm_normalize(Matrix *p);
Matrix *pm_convolve(Matrix *p, int n, double epsilon);
Matrix **pm_convolve_save(Matrix *p, int n, double epsilon);
Matrix *pm_convolve_many(Matrix **p, int *counts, int n, double epsilon);
Matrix *pm_convolve_many_fast(Matrix **p, int n, int max_nrows, int max_ncols);
Matrix *pm_convolve_fast(Matrix *p, int n, double epsilon);
void pm_stats(Matrix *p, double *mean_x, double *mean_y, double *var_x, 
              double *var_y, double *covar);
Vector *pm_x_given_tot_bvn(int tot, double mu_x, double mu_y, 
                           double sigma_x, double sigma_y, double rho);
Vector *pm_y_given_tot_bvn(int tot, double mu_x, double mu_y, 
                           double sigma_x, double sigma_y, double rho);
Vector *pm_x_given_tot_indep(int tot, Vector *marg_x, Vector *marg_y);
Vector *pm_y_given_tot_indep(int tot, Vector *marg_x, Vector *marg_y);

#endif
