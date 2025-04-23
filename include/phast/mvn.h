/***************************************************************************
 * PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#ifndef MVN_H
#define MVN_H

#include <stdio.h>
#include <phast/matrix.h>

enum mvn_type {MVN_STD, MVN_IDENTITY, MVN_DIAG, MVN_LOWR, MVN_GEN};

typedef struct MVN MVN;

struct MVN {
  int dim; /* dimensionality */
  Vector *mu; /* mean vector */
  Matrix *sigma;   /* covariance matrix*/
  Matrix *cholL;  /* lower triangular matrix from Cholesky
                     decomposition of covariance matrix */
  Vector *evals;  /* eigendecomposition of covariance matrix */
  Matrix *evecs; /* for cases in which Cholesky cannot be obtained */
  Matrix *lowR; /* for low-rank parameterization */
  MVN *lowRmvn; /* embedded low-rank MVN */
  enum mvn_type type;
};

MVN *mvn_new(int dim, Vector *mu, Matrix *sigma);

MVN *mvn_new_LOWR(int dim, Vector *mu, Matrix *R);

void mvn_reset_LOWR(MVN *mvn);

void mvn_free(MVN *mvn);

void mvn_update_type(MVN *mvn);

void mvn_sample_std(Vector *retval);

void mvn_map_std(MVN *mvn, Vector *rv, Vector *lowrv);

void mvn_sample(MVN *mvn, Vector *retval);

void mvn_sample_anti(MVN *mvn, Vector *retval1, Vector *retval2);

void mvn_sample_anti_keep(MVN *mvn, Vector *retval1, Vector *retval2,
                          Vector *origstd);

void mvn_preprocess(MVN *mvn, unsigned int force_eigen);

double mvn_log_dens(MVN *mvn, Vector *x);

double mvn_log_det(MVN *mvn);

void mvn_print(MVN *mvn, FILE *F);

void mvn_rederive_std(MVN *mvn, Vector *x, Vector *x_std);

double mvn_trace(MVN *mvn);

double mvn_mu2(MVN *mvn);

void mvn_project_LOWR(MVN *mvn, Vector *z, Vector *a);

#endif
