/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* Multivariate normal distributions */


#include <phast/mvn.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <phast/misc.h>

/* Create a new MVN object of specified dimension.  If mean is NULL,
   initialize to vector of zeroes.  If covariance matrix sigma is
   NULL, initialize to identity matrix. */
MVN *mvn_new(int dim, Vector *mu, Matrix *sigma, enum mvn_type type) {
  MVN *mvn = smalloc(sizeof(MVN));
  mvn->dim = dim;

  if (mu == NULL) { /* standard normal case */
    if (sigma != NULL)
      die("ERROR in mvn_new: if mu is NULL then sigma must be NULL also.\n");

    mvn->mu = vec_new(dim);
    vec_zero(mvn->mu);
    mvn->type = MVN_STD;
  }
  else {   /* mu is specified; sigma may or may not be specified */
    if (mu->size != dim)
      die("ERROR in mvn_new: bad dimension in mean vector.\n");
    mvn->mu = mu;
  }
  
  if (sigma == NULL) {  /* use identity covariance */
      mvn->sigma = mat_new(dim, dim);
      mat_set_identity(mvn->sigma);
      if (mu != NULL)  /* if mu is NULL, already set to MVN_STD */
        mvn->type = MVN_IDENTITY;
  }    
  else {  /* sigma is specified; mu may or may not be */
    if (sigma->nrows != dim || sigma->ncols != dim)
      die("ERROR in mvn_new: bad dimension in covariance matrix.\n");
    mvn->sigma = sigma;
    mvn->type = type;  /* note user may still specify something other than MVN_GEN */
  }

  mvn->cholL = NULL;   /* will be updated as needed */

  return mvn;
}
             
/* Sample a vector from a standard multivariate normal distribution,
   with zero mean and identity covariance.  */
void mvn_sample_std(Vector *retval) {
  int i;
  static int seeded = 0;
  double u1, u2, z1, z2;
  
  if (!seeded) {
    srandom((unsigned int)time(NULL));
    seeded = 1;
  }

  /* draw indep samples from standard normal using Box-Muller transform */
  for (i = 0; i < retval->size; i += 2) {
    u1 = unif_rand();
    u2 = unif_rand();
    z1 = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
    z2 = sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
    vec_set(retval, i, z1);
    if (i+1 < retval->size)
      vec_set(retval, i+1, z2);
  }  
}

/* Sample a vector from a multivariate normal distribution with mean
   mu and covariance sigma. */
void mvn_sample(MVN *mvn, Vector *retval) {
  int i, j;

  if (mvn->dim != retval->size)
    die("ERROR in mvn_sample: bad dimensions.\n");
  
  mvn_sample_std(retval);

  /* do nothing if mvn->type == MVN_STD */

  if (mvn->type == MVN_IDENTITY)
    vec_plus_eq(retval, mvn->mu);  
  else if (mvn->type == MVN_DIAG) {
    for (i = 0; i < mvn->dim; i++) 
      vec_set(retval, i, vec_get(mvn->mu, i) + sqrt(mat_get(mvn->sigma, i, i)) *
              vec_get(retval, i));
  }
  else {   /* general covariance matrix.  Assume Cholesky already
              updated (responsibility of calling code) */
    for (i = 0; i < mvn->dim; i++) {
      double covarsum = 0;
      for (j = 0; j < mvn->dim; j++) {
        covarsum += mat_get(mvn->cholL, i, j) * vec_get(retval, j);    /* CHECK!!! */
      }
      vec_set(retval, i, vec_get(mvn->mu, i) + covarsum);
    }
  }       
}

void mvn_update_cholesky(MVN *mvn) {
  int retval = mat_cholesky(mvn->cholL, mvn->sigma);
  if (retval != 0)
    die("ERROR in nj_laplacian_pinv. Cannot compute Cholesky decomposition of Laplacian pseudoinverse.\n");
}

/* return log density function of MVN for a given vector x. */
double mvn_log_dens(MVN *mvn, Vector *x) {
  double retval = -x->size/2 * log(2 * M_PI) - 0.5 * mvn_log_det(mvn);
  int i;

  if (x->size != mvn->dim)
    die("ERROR in mvn_dens: bad dimension\n");

  if (mvn->type == MVN_STD || mvn->type == MVN_IDENTITY || mvn->type == MVN_DIAG) {
    for (i = 0; i < x->size; i++) {
      retval -= 0.5 * pow(vec_get(x, i) - vec_get(mvn->mu, i), 2) *
        mat_get(mvn->sigma, i, i);
    }
  }
  else {
    /* compute the quadratic form using the Cholesky lower triangular matrix */
    Vector *z = vec_create_copy(x), *y = vec_new(x->size);
    vec_minus_eq(z, mvn->mu);
    mat_forward_subst(mvn->cholL, z, y);
    retval -= 0.5 * log(vec_norm(y));
    vec_free(z); vec_free(y);
  }

  return retval;
}

/* return log determinant of covariance matrix */
double mvn_log_det(MVN *mvn) {
  int i;
  double retval = 0;  /* this will be the answer if type MVN_STD or MVN_IDENTITY */

  if (mvn->type == MVN_DIAG) {
    for (i = 0; i < mvn->dim; i++)
      retval += log(mat_get(mvn->sigma, i, i));
  }
  else if (mvn->type == MVN_GEN) { /* in this case, use the Cholesky factorization */
    if (mvn->cholL == NULL)
      die("ERROR in mvn_log_det: Cholesky factorization not avaialable.\n");

    for (i = 0; i < mvn->dim; i++)
      retval += 2*log(mat_get(mvn->cholL, i, i));   /* determinant is sum of squares on main diag */
  }

  return retval;
}
