/* PHylogenetic Analysis with Space/Time models
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

enum mvn_type {MVN_STD, MVN_IDENTITY, MVN_DIAG, MVN_GEN};

typedef struct {
  int dim; /* dimensionality */
  Vector *mu; /* mean vector */
  Matrix *sigma;   /* covariance matrix*/
  Matrix *cholL;  /* lower triangular matrix from Cholesky
                     decomposition of covariance matrix */
  Vector *evals;  /* eigendecomposition of covariance matrix */
  Matrix *evecs; /* for cases in which Cholesky cannot be obtained */
  enum mvn_type type;
} MVN;

MVN *mvn_new(int dim, Vector *mu, Matrix *sigma);

void mvn_free(MVN *mvn);

void mvn_update_type(MVN *mvn);

void mvn_sample_std(Vector *retval);

void mvn_sample(MVN *mvn, Vector *retval);

void mvn_preprocess(MVN *mvn, unsigned int force_eigen);

double mvn_log_dens(MVN *mvn, Vector *x);

double mvn_log_det(MVN *mvn);

void mvn_print(MVN *mvn, FILE *F);

#endif
