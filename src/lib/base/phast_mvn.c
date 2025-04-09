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
#include <phast/eigen.h>

/* Create a new MVN object of specified dimension.  If mean is NULL,
   initialize to vector of zeroes.  If covariance matrix sigma is
   NULL, initialize to identity matrix. */
MVN *mvn_new(int dim, Vector *mu, Matrix *sigma) {
  MVN *mvn = smalloc(sizeof(MVN));
  mvn->dim = dim;

  assert(mvn->dim > 0);
  
  if (mu == NULL) { /* standard normal case */
    mvn->mu = vec_new(dim);
    vec_zero(mvn->mu);
  }
  else {  
    if (mu->size != dim)
      die("ERROR in mvn_new: bad dimension in mean vector.\n");
    mvn->mu = mu;
  }
  
  if (sigma == NULL) {  /* use identity covariance */
      mvn->sigma = mat_new(dim, dim);
      mat_set_identity(mvn->sigma);
  }    
  else {  
    if (sigma->nrows != dim || sigma->ncols != dim)
      die("ERROR in mvn_new: bad dimension in covariance matrix.\n");
    mvn->sigma = sigma;
  }

  mvn_update_type(mvn);
  mvn->cholL = NULL;   /* will be updated as needed */
  mvn->evals = NULL;
  mvn->evecs = NULL;
  
  return mvn;
}

MVN *mvn_copy(MVN *mvn) {
  return mvn_new(mvn->dim, vec_create_copy(mvn->mu), mat_create_copy(mvn->sigma));
}

void mvn_free(MVN *mvn) {
  vec_free(mvn->mu);
  mat_free(mvn->sigma);
  if (mvn->cholL != NULL)
    mat_free(mvn->cholL);
  if (mvn->evals != NULL)
    vec_free(mvn->evals);
  if (mvn->evecs != NULL)
    mat_free(mvn->evecs);
  free(mvn);
}

void mvn_update_type(MVN *mvn) {
  int mean_zero = TRUE, diag_covar = TRUE, ident_covar = TRUE;
  int i, j;
  
  for (i = 0; mean_zero == TRUE && i < mvn->dim; i++)
    if (vec_get(mvn->mu, i) != 0)
      mean_zero = FALSE;

  for (i = 0; diag_covar == TRUE && i < mvn->dim; i++) 
    for (j = 0; diag_covar == TRUE && j < i; j++) 
      if (mat_get(mvn->sigma, i, j) != 0 || mat_get(mvn->sigma, j, i) != 0)
        diag_covar = FALSE;

  if (diag_covar == FALSE)
    ident_covar = FALSE;
  else  /* diagonal; now check whether identity */
    for (i = 0; ident_covar && i < mvn->dim; i++)
      if (mat_get(mvn->sigma, i, i) != 1)
        ident_covar = FALSE;
      
  if (mean_zero && ident_covar)
    mvn->type = MVN_STD;
  else if (ident_covar)
    mvn->type = MVN_IDENTITY;
  else if (diag_covar)
    mvn->type = MVN_DIAG;
  else
    mvn->type = MVN_GEN;

  /* notice that if dim == 1, MVN_GEN is not possible but the other three are possible */
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

  printf("orig sample: ");
  vec_print(retval, stdout);
  
  /* do nothing if MVN_STD */
  
  if (mvn->type == MVN_IDENTITY)
    vec_plus_eq(retval, mvn->mu);  
  else if (mvn->type == MVN_DIAG) {
    for (i = 0; i < mvn->dim; i++) 
      vec_set(retval, i, vec_get(mvn->mu, i) + sqrt(mat_get(mvn->sigma, i, i)) *
              vec_get(retval, i));
  }
  else if (mvn->type == MVN_GEN) {
    /* general covariance matrix. Assume Cholesky or eigencomposition
       already updated (responsibility of calling code) */
    Vector *tmp = vec_create_copy(retval);

    if (mvn->cholL != NULL) { /* use Cholesky if possible */
      for (i = 0; i < mvn->dim; i++) {
        double covarsum = 0;
        for (j = 0; j <= i; j++) 
          covarsum += mat_get(mvn->cholL, i, j) * vec_get(tmp, j);
        vec_set(retval, i, vec_get(mvn->mu, i) + covarsum);
      }
    }
    else if (mvn->evals != NULL) {
      for (i = 0; i < mvn->dim; i++) {
        double covarsum = 0;
        for (j = 0; j < mvn->dim; j++) 
          covarsum += mat_get(mvn->evecs, i, j) * sqrt(vec_get(mvn->evals, j)) * vec_get(tmp, j); 
        vec_set(retval, i, vec_get(mvn->mu, i) + covarsum);
      }
    }
    else
      die("ERROR in mvn_sample: must have either Cholesky or eigendecomposition.  Call mvn_preprocess first.\n");

    vec_free(tmp);
  }       
}

/* pre-calculate Cholesky decomposition or eigendecomposition of
   covariance matrix for faster sampling and density calculation.  Use
   force_eigen = TRUE to force the eigendecomposition */
void mvn_preprocess(MVN *mvn, unsigned int force_eigen) {
  int retval = 1;

  if (force_eigen == FALSE) {
    if (mvn->cholL == NULL)
      mvn->cholL = mat_new(mvn->dim, mvn->dim);

    retval = mat_cholesky(mvn->cholL, mvn->sigma);
  
    if (retval == 0) {
      /* Cholesky successful.  Don't bother with eigendecomposition but
         make sure values are NULL so an old one doesn't get used */
      if (mvn->evals != NULL) {
        vec_free(mvn->evals);
        mvn->evals = NULL;
      }
      if (mvn->evecs != NULL) {
        mat_free(mvn->evecs);
        mvn->evecs = NULL;
      }
    }
  }

  if (force_eigen == TRUE || retval != 0) {
    /* do eigendecomposition instead */
    int i;

    if (mvn->cholL != NULL) {
      mat_free(mvn->cholL);
      mvn->cholL = NULL;
    }

    if (mvn->evals == NULL)
      mvn->evals = vec_new(mvn->dim);
    if (mvn->evecs == NULL)
      mvn->evecs = mat_new(mvn->dim, mvn->dim);

    mat_diagonalize_sym(mvn->sigma, mvn->evals, mvn->evecs);

    /* prohibit zero eigenvalues */
    for (i = 0; i < mvn->dim; i++) {
      if (fabs(vec_get(mvn->evals, i)) < 1e-6)
        vec_set(mvn->evals, i, 1e-6);
      if (vec_get(mvn->evals, i) < 0)
        die("ERROR in mvn_preprocess: covariance matrix is not positive definite.\n");
    }
  }
}

/* return log density function of MVN for a given vector x. */
double mvn_log_dens(MVN *mvn, Vector *x) {
  double retval = -x->size/2.0 * log(2 * M_PI) - 0.5 * mvn_log_det(mvn);
  int i, j;

  if (x->size != mvn->dim)
    die("ERROR in mvn_dens: bad dimension\n");

  if (mvn->type == MVN_STD || mvn->type == MVN_IDENTITY || mvn->type == MVN_DIAG) {
    for (i = 0; i < x->size; i++) {
      retval -= 0.5 * pow(vec_get(x, i) - vec_get(mvn->mu, i), 2) /
        mat_get(mvn->sigma, i, i);
    }
  }
  else {
    /* more complicated for MVN_GEN */
    Vector *z = vec_create_copy(x), *y = vec_new(x->size);
    double dotprod = 0;
    vec_minus_eq(z, mvn->mu);

    if (mvn->cholL != NULL) {
      /* compute the quadratic form by forward substitution from the Cholesky lower triangular matrix */     
      mat_forward_subst(mvn->cholL, z, y);
      for (i = 0; i < mvn->dim; i++)
        dotprod += pow(vec_get(y, i), 2);
    }
    else if (mvn->evals != NULL) {
      /* compute the quadratic form from the eigendecomposition */
      vec_zero(y);
      for (i = 0; i < mvn->dim; i++) {
        for (j = 0; j < mvn->dim; j++)
          vec_set(y, i, vec_get(y, i) + mat_get(mvn->evecs, j, i) * vec_get(z, j));  /* note must use transpose of evecs here */
        dotprod += pow(vec_get(y, i), 2) / vec_get(mvn->evals, i); /* weight by inverse eigenvalue */
      }
    }
    else
      die("ERROR in mvn_log_dens: must have either Cholesky or eigendecomposition.  Call mvn_preprocess first.\n");
    retval -= 0.5 * dotprod;
    vec_free(z); vec_free(y);
  }

  return retval;
}

/* return log determinant of covariance matrix */
double mvn_log_det(MVN *mvn) {
  int i;
  double retval = 0;  /* this will be the answer if type MVN_STD or MVN_IDENTITY */

  if (mvn->type == MVN_DIAG) { /* easy if diagonal */
    for (i = 0; i < mvn->dim; i++)
      retval += log(mat_get(mvn->sigma, i, i)); 
  }
  else if (mvn->type == MVN_GEN) { /* harder in this case */
    if (mvn->cholL != NULL) { /* use Cholesky if available */
      for (i = 0; i < mvn->dim; i++)
        retval += 2*log(mat_get(mvn->cholL, i, i));   /* log determinant is sum of log squares on main diag */
    }
    else if (mvn->evals != NULL) { /* otherwise use eigenvalues */
      for (i = 0; i < mvn->dim; i++)
        retval += log(vec_get(mvn->evals, i));   /* log determinant is sum of log eigenvalues */
    }
    else
      die("ERROR in mvn_log_det: must have either Cholesky or eigendecomposition.  Call mvn_preprocess first.\n");
  }

  return retval;
}

void mvn_print(MVN *mvn, FILE *F){
  fprintf(F, "MVN (dim %d, type %d)\n", mvn->dim, mvn->type);
  fprintf(F, "mu: ");
  vec_print(mvn->mu, F);
  fprintf(F, "sigma:\n");
  mat_print(mvn->sigma, F);
  if (mvn->cholL != NULL) {
    fprintf(F, "cholL:\n");
    mat_print(mvn->cholL, F);
  }
  if (mvn->evals != NULL) {
    fprintf(F, "evals: ");
    vec_print(mvn->evals, F);
  }
  if (mvn->evecs != NULL) {
    fprintf(F, "evecs:\n");
    mat_print(mvn->evecs, F);
  }
}

/* obtain the underlying standard mvn used to generate a
   non-standard mvn.  Useful for calculations based on the
   reparameterization trick */ 
void mvn_rederive_std(MVN *mvn, Vector *x, Vector *x_std) {
  int i, j;
  if (mvn->dim != x->size || x->size != x_std->size)
    die("ERROR in mvn_rederive_std: bad dimension.\n");
  if (mvn->type == MVN_GEN) {
    Vector *y = vec_create_copy(x);
    vec_minus_eq(y, mvn->mu);
    if (mvn->cholL != NULL)  /* can do by forward substitution */
      mat_forward_subst(mvn->cholL, y, x_std);
    else if (mvn->evals != NULL) { /* here need full matrix mult */
      vec_zero(x_std);
      for (i = 0; i < mvn->dim; i++) {
        for (j = 0; j < mvn->dim; j++) {
          double val = mat_get(mvn->evecs, j, i) * vec_get(y, j);
          /* note implicit transpose of eigenvector matrix */
          vec_set(x_std, i, vec_get(x_std, i) + val);
        }
        vec_set(x_std, i, vec_get(x_std, i)/sqrt(vec_get(mvn->evals, i)));
      }
    }
    else
      die("ERROR in mvn_derive_std: must have either Cholesky or eigendecomposition.  Call mvn_preprocess first.\n");
    vec_free(y);
  }
  else { /* in diagonal cases can be done element by element */
    for (i = 0; i < mvn->dim; i++)
      vec_set(x_std, i, (vec_get(x, i) - vec_get(mvn->mu, i)) / sqrt(mat_get(mvn->sigma, i, i)));
  }
}

/* return trace of covariance matrix */
double mvn_trace(MVN *mvn) {
  double retval = 0;
  int i;
  for (i = 0; i < mvn->dim; i++)
    retval += mat_get(mvn->sigma, i, i);
  return retval;
}

/* return mu \dot mu */
double mvn_mu2(MVN *mvn) {
  double retval = 0;
  int i;
  for (i = 0; i < mvn->dim; i++)
    retval += pow(vec_get(mvn->mu, i), 2);
  return retval;
}
