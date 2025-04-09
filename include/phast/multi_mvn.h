/* PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#ifndef MMVN_H
#define MMVN_H

#include <stdio.h>
#include <phast/matrix.h>
#include <phast/mvn.h>

/* bundle of MVNs sharing the same covariance matrix */
typedef struct {
  int n; /* number of points */
  int d; /* number of dimensions, equal to number of component MVNs */
  Vector **mu; /* array of mean vectors */
  MVN *mvn; /* MVN object with shared covariance and associated
               preprocessing */
  enum mvn_type type;
} multi_MVN;

multi_MVN *mmvn_new(int n, int d, enum mvn_type type);

void mmvn_sample(multi_MVN *mmvn, Vector *retval);

void mmvn_combine_means(multi_MVN *mmvn, Vector *mu);

double mmvn_log_dens(multi_MVN *mmvn, Vector *x);

double mmvn_log_det(multi_MVN *mmvn);

void mmvn_project_down(multi_MVN *mmvn, Vector *x_full,
                               Vector *x_d, int d);

void mmvn_project_up(multi_MVN *mmvn, Vector *x_d,
                             Vector *x_full, int d);

void mmvn_rederive_std(multi_MVN *mmvn, Vector *points,
                               Vector *points_std);

double mmvn_trace(multi_MVN *mmvn);

double mmvn_mu2(multi_MVN *mmvn);

void mmvn_save_mu(multi_MVN *mmvn, Vector *mu_saved);

void mmvn_restore_mu(multi_MVN *mmvn, Vector *mu_saved);

void mmvn_print(multi_MVN *mmvn, FILE *F, unsigned int in_line,
                        unsigned int do_covariance);

double mmvn_get_mu_el(multi_MVN *mmvn, int i);

void mmvn_set_mu_el(multi_MVN *mmvn, int i, double val);

#endif
