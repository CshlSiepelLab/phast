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
#include <assert.h>
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
  mvn->lowR = NULL;
  mvn->lowRmvn = NULL;
  mvn->lowR_invRtR = NULL;
  
  return mvn;
}

MVN *mvn_new_LOWR(int dim, Vector *mu, Matrix *R) {
  MVN *retval = mvn_new(dim, mu, NULL);
  if (dim != mu->size || mu->size != R->nrows || R->ncols > R->nrows)
    die("ERROR in mvn_new_LOWR: bad dimension.\n");
  retval->type = MVN_LOWR;
  mvn_reset_LOWR(retval, R);
  return retval;
}

/* reset covariance matrix and other auxiliary data based on low-rank
   representation */
void mvn_reset_LOWR(MVN *mvn, Matrix *R) {
  assert(mvn->type == MVN_LOWR && mvn->sigma->nrows == R->nrows &&
         R->ncols <= R->nrows);

  if (mvn->lowR == NULL)
    mvn->lowR = mat_create_copy(R);
  else
    mat_copy(mvn->lowR, R);
  
  if (mvn->lowRmvn == NULL) 
    mvn->lowRmvn = mvn_new(mvn->lowR->ncols, NULL, NULL);
  else
    assert(mvn->lowRmvn->dim == mvn->lowR->ncols);

  if (mvn->lowR_invRtR == NULL)
    mvn->lowR_invRtR = mat_new(mvn->lowR->ncols, mvn->lowR->ncols);
  else
    assert(mvn->lowR_invRtR->nrows == mvn->lowR->ncols);

  mvn->lowRmvn->type = MVN_GEN;
  mat_set_gram(mvn->sigma, mvn->lowR);
  mat_set_gram_col(mvn->lowRmvn->sigma, mvn->lowR);
  mvn_preprocess(mvn->lowRmvn, FALSE); /* will do Cholesky */
  mvn_invert_RtR_LOWR(mvn); /* uses Cholesky */
}

MVN *mvn_copy(MVN *mvn) {
  if (mvn->type == MVN_LOWR)
    return mvn_new_LOWR(mvn->dim, vec_create_copy(mvn->mu), mat_create_copy(mvn->lowR));
  else
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
  if (mvn->lowR != NULL)
    mat_free(mvn->lowR);
  if (mvn->lowRmvn != NULL)
    mvn_free(mvn->lowRmvn);
  if (mvn->lowR_invRtR != NULL)
    mat_free(mvn->lowR_invRtR);
  free(mvn);
}

void mvn_update_type(MVN *mvn) {
  int mean_zero = TRUE, diag_covar = TRUE, ident_covar = TRUE;
  int i, j;

  assert(mvn->sigma != NULL);

  if (mvn->mu == NULL)
    mean_zero = FALSE; /* can happen with multi-MVNs */
  else
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
  /* note that MVN_LOWR will never be used by this function; has to be set separately */
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
  if (mvn->dim != retval->size)
    die("ERROR in mvn_sample: bad dimensions.\n");

  Vector *lowrv = NULL;

  if (mvn->type == MVN_LOWR) { /* in this case sample in the lower
                                  dim */
    lowrv = vec_new(mvn->lowR->ncols);
    mvn_sample_std(lowrv);
  }
  else 
    mvn_sample_std(retval);    

  if (mvn->type != MVN_STD)  /* already correct if MVN_STD */
    mvn_map_std(mvn, retval, lowrv);

  if (lowrv != NULL)
    vec_free(lowrv);
}

/* Like mvn_sample, but sample a pair of antithetic random variates to
   reduce variance */
void mvn_sample_anti(MVN *mvn, Vector *retval1, Vector *retval2) {
  if (mvn->dim != retval1->size || mvn->dim != retval2->size)
    die("ERROR in mvn_sample_anti: bad dimensions.\n");

  if (mvn->type != MVN_LOWR) {
    mvn_sample_std(retval1);
    vec_copy(retval2, retval1);
    vec_scale(retval2, -1.0);

    if (mvn->type != MVN_STD) {
      mvn_map_std(mvn, retval1, NULL);
      mvn_map_std(mvn, retval2, NULL);
    }
  }
  else { /* LOWR case -- sample in k dim and project up */
    int k = mvn->lowR->ncols;
    Vector *proj1 = vec_new(k), *proj2 = vec_new(k);
    mvn_sample_std(proj1);
    vec_copy(proj2, proj1);
    vec_scale(proj2, -1.0);
    mvn_map_std(mvn, retval1, proj1);
    mvn_map_std(mvn, retval2, proj2);
    vec_free(proj1);
    vec_free(proj2);
  }
}

/* Like mvn_sample_anti, but keeps track of the original standard
   normal variate for use in downstream calculations */
void mvn_sample_anti_keep(MVN *mvn, Vector *retval1,
                          Vector *retval2, Vector *origstd) {
  if (mvn->dim != retval1->size || mvn->dim != retval2->size)
    die("ERROR in mvn_sample_anti_keep: bad dimensions.\n");

  if (mvn->type == MVN_LOWR && origstd->size != mvn->lowR->ncols)
    die("ERROR in mvn_sample_anti_keep: origstd must have low-rank dimension in MVN_LOWR case.\n");
  else if (mvn->type != MVN_LOWR && origstd->size != mvn->dim)
    die("ERROR in mvn_sample_anti_keep: bad dimension in origstd.\n");
  
  if (mvn->type != MVN_LOWR) {
    mvn_sample_std(retval1);

    vec_copy(origstd, retval1);
    vec_copy(retval2, retval1);
    vec_scale(retval2, -1.0);

    if (mvn->type != MVN_STD) {
      mvn_map_std(mvn, retval1, NULL);
      mvn_map_std(mvn, retval2, NULL);
    }
  }
  else { /* LOWR case -- sample in k dim and project up */
    int k = mvn->lowR->ncols;
    Vector *proj1 = vec_new(k), *proj2 = vec_new(k);
    mvn_sample_std(proj1);
    vec_copy(origstd, proj1);
    vec_copy(proj2, proj1);
    vec_scale(proj2, -1.0);
    mvn_map_std(mvn, retval1, proj1);
    mvn_map_std(mvn, retval2, proj2);
    vec_free(proj1);
    vec_free(proj2);
  }
}

/* Map a standard MVN variable to a general MVN variable.  On input,
   rv should be a draw from mvn_sample_std, on output it will be a
   random variate from the distribution defined by mvn.  In the
   MVN_LOWR case, must supply a separate input rv of the lower
   dimension k (lowrv).  In other cases, pass NULL for this
   variable. */
void mvn_map_std(MVN *mvn, Vector *rv, Vector *lowrv) {
  int i, j;
  
  if (mvn->type == MVN_IDENTITY)
    vec_plus_eq(rv, mvn->mu);  
  else if (mvn->type == MVN_DIAG) {
    for (i = 0; i < mvn->dim; i++) 
      vec_set(rv, i, vec_get(mvn->mu, i) + sqrt(mat_get(mvn->sigma, i, i)) *
              vec_get(rv, i));
  }
  else if (mvn->type == MVN_GEN) {
    /* general covariance matrix. Assume Cholesky or eigencomposition
       already updated (responsibility of calling code) */
    Vector *tmp = vec_create_copy(rv);

    if (mvn->cholL != NULL) { /* use Cholesky if possible */
      for (i = 0; i < mvn->dim; i++) {
        double covarsum = 0;
        for (j = 0; j <= i; j++) 
          covarsum += mat_get(mvn->cholL, i, j) * vec_get(tmp, j);
        vec_set(rv, i, vec_get(mvn->mu, i) + covarsum);
      }
    }
    else if (mvn->evals != NULL) {
      for (i = 0; i < mvn->dim; i++) {
        double covarsum = 0;
        for (j = 0; j < mvn->dim; j++) 
          covarsum += mat_get(mvn->evecs, i, j) * sqrt(vec_get(mvn->evals, j)) * vec_get(tmp, j);
        
        vec_set(rv, i, vec_get(mvn->mu, i) + covarsum);
      }
    }
    else
      die("ERROR in mvn_map_std: must have either Cholesky or eigendecomposition.  Call mvn_preprocess first.\n");
    
    vec_free(tmp);
  }
  else {
    int k = mvn->lowR->ncols;
    assert(mvn->type == MVN_LOWR);
    if (lowrv == NULL || lowrv->size != k) 
      die("ERROR in mvn_map_std: must supply lowrv of correct dimension in MVN_LOWR case.\n");
    mat_vec_mult(rv, mvn->lowR, lowrv); /* project up via lowR */
    vec_plus_eq(rv, mvn->mu); /* add mu */
  }
}

/* pre-calculate Cholesky decomposition or eigendecomposition of
   covariance matrix for faster sampling and density calculation.  Use
   force_eigen = TRUE to force the eigendecomposition */
void mvn_preprocess(MVN *mvn, unsigned int force_eigen) {
  int retval = 1;

  if (mvn->type != MVN_GEN)
    return; /* NOTE: includes MVN_LOWR; in this case we don't want to
               preprocess the outer MVN, only the embedded one */
  
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

    if (mat_diagonalize_sym(mvn->sigma, mvn->evals, mvn->evecs) != 0)
      die("ERROR in mvn_preprocess: matrix diagonalization failed.\n");

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
    die("ERROR in mvn_log_dens: bad dimension\n");

  if (mvn->type == MVN_STD || mvn->type == MVN_IDENTITY || mvn->type == MVN_DIAG) {
    for (i = 0; i < x->size; i++) {
      retval -= 0.5 * pow(vec_get(x, i) - vec_get(mvn->mu, i), 2) /
        mat_get(mvn->sigma, i, i);
    }
  }
  else if (mvn->type == MVN_LOWR) {
    /* in this case, we need to project x down to k dimensions and
       compute its density under the embedded low-rank MVN */
    Vector *z = vec_create_copy(x), *a = vec_new(mvn->lowR->ncols);
    vec_minus_eq(z, mvn->mu);
    mvn_project_LOWR(mvn, z, a);
    retval = mvn_log_dens(mvn->lowRmvn, a);
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
  else if (mvn->type == MVN_LOWR) /* just delegate to low-rank MVN */
    retval = mvn_log_det(mvn->lowRmvn);
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
  if (mvn->mu == NULL)
    fprintf(F, "<NULL>\n");
  else
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
  else if (mvn->type != MVN_LOWR) { /* in diagonal cases can be done element by element */
    for (i = 0; i < mvn->dim; i++)
      vec_set(x_std, i, (vec_get(x, i) - vec_get(mvn->mu, i)) / sqrt(mat_get(mvn->sigma, i, i)));
  }
  else if (mvn->type == MVN_LOWR)
    /* we'll punt in this case for now, don't need */
    die("ERROR in mvn_rederive_std: MVN_LOWR case not supported.\n");
}

/* return trace of covariance matrix */
double mvn_trace(MVN *mvn) {
  double retval = 0;
  int i;

  if (mvn->type == MVN_LOWR) /* in this case we can just use the trace of
                                the embedded low-rank MVN */
    return mvn_trace(mvn->lowRmvn);
  
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

/* project a vector in n dimensions representing a zero-mean MVN with
   covariance R x R^T down to the embedded k-dimensional MVN that lies
   in the column space of R.  The new vector a can be assumed to be an
   MVN with zero mean and covariance R^T x R.  */
void mvn_project_LOWR(MVN *mvn, Vector *z, Vector *a) {
  int n = mvn->lowR->nrows, k = mvn->lowR->ncols;
  Vector *b = vec_new(k);
  
  assert(z->size == n && a->size == k && mvn->lowRmvn != NULL &&
         mvn->lowR_invRtR != NULL);
    
  mat_vec_mult_transp(b, mvn->lowR, z); /* start by computing b = R^T z */

  /* now we can just pre-multiply by the precomputed inverse of R^T x R */
  mat_vec_mult(a, mvn->lowR_invRtR, b);

  vec_free(b);
}

/* compute the inverse of R^T x R for use in projecting n-dimensional
   vectors down to the k-dimensional subspace spanned by R.  Can be
   precomputed and stored for efficiency.  Uses Cholesky decomposition
   of R^T x R. */
void mvn_invert_RtR_LOWR(MVN *mvn) {
  int k = mvn->lowR->ncols;
  Vector *b = vec_new(k), *y = vec_new(k), *x = vec_new(k);
  assert(mvn->lowRmvn->cholL != NULL);
  
  for (int j = 0; j < k; j++) { /* columns */
    vec_zero(b);
    vec_set(b, j, 1.0);
    mat_forward_subst(mvn->lowRmvn->cholL, b, y);  /* solve L y = b for y */
    mat_backward_subst(mvn->lowRmvn->cholL, y, x); /* solve L^T x = y for x */

    /* fill out appropriate column of inv matrix */
    for (int i = 0; i < k; i++)   /* rows */
      mat_set(mvn->lowR_invRtR, i, j, vec_get(x, i));
  }
  vec_free(b);
  vec_free(y);
  vec_free(x);
}
