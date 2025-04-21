/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "phast/misc.h"
#include "phast/stringsplus.h"
#include "phast/mvn.h"
#include "phast/multi_mvn.h"

/* This is a wrapper for the MVN class that gains efficiency by
   representing a large MVN as a collection of smaller MVNs that share
   the same covariance structure but may have different means.  It is
   useful for representing embeddings of points in a d-dimensional
   space such that the d coordinate axes are independent, but along
   each axis the points share the same covariance structure.
   Efficiency is gained by keeping the covariance matrix small and
   reusing precomputed Cholesky or eigendecompositions.  The routines
   are written in a general enough manner to handle diagonal
   covariance matrices as well, in which case the underlying MVN is
   left intact (because there are no real gains in efficiency in this
   case) */

/* create a new multi-MVN of desired type from scratch, for n points
   and dimension d. If type != MVN_GEN, will just wrap a single MVN */
multi_MVN *mmvn_new(int n, int d, enum mvn_type type) {
  int i;
  multi_MVN *retval = smalloc(sizeof(multi_MVN));

  if (type != MVN_GEN) {
    retval->n = n*d;
    retval->d = 1;
  }
  else {
    retval->n = n;
    retval->d = d;
  }

  retval->mvn = mvn_new(retval->n, NULL, NULL);
  retval->mvn->type = type;
  retval->type = type;
  
  if (type == MVN_GEN) {
    /* clear the orig mu vector; will handle it just by pointer swapping */
    vec_free(retval->mvn->mu);
    retval->mvn->mu = NULL;
    retval->mu = smalloc(d * sizeof(Vector*));
    for (i = 0; i < d; i++)
      retval->mu[i] = vec_new(retval->n);
  }
  else 
    retval->mu = NULL;
  
  return retval;
}

/* set the mean parameters */
void mmvn_set_mu(multi_MVN *mmvn, Vector *mu) {
  assert(mu->size == mmvn->d * mmvn->n);
  if (mmvn->type != MVN_GEN)
    vec_copy(mmvn->mvn->mu, mu);
  else {
    for (int d = 0; d < mmvn->d; d++)
      mmvn_project_down(mmvn, mu, mmvn->mu[d], d);
  }
}

/* set shared covariance */
void mmvn_set_sigma(multi_MVN *mmvn, Matrix *sigma) {
  assert(sigma->nrows == mmvn->n &&
         sigma->nrows == sigma->ncols);
  mat_copy(mmvn->mvn->sigma, sigma);
  mmvn_preprocess(mmvn, TRUE);
}

void mmvn_preprocess(multi_MVN *mmvn, unsigned int force_eigen) {
  mvn_preprocess(mmvn->mvn, force_eigen);
}

void mmvn_sample(multi_MVN *mmvn, Vector *retval) {
  if (mmvn->type != MVN_GEN)  /* in this case it's just as efficient
                                 to use the full MVN */
    mvn_sample(mmvn->mvn, retval);
  else {
    Vector *xcomp = vec_new(mmvn->n);
    int d;
    assert(retval->size == mmvn->d * mmvn->n);
    for (d = 0; d < mmvn->d; d++) {
      mmvn->mvn->mu = mmvn->mu[d]; /* swap in the appropriate mean */
      mvn_sample(mmvn->mvn, xcomp);
      mmvn_project_up(mmvn, xcomp, retval, d);
    }
    vec_free(xcomp);
  }
}

void mmvn_sample_anti(multi_MVN *mmvn, Vector *retval1, Vector *retval2) {
  if (mmvn->type != MVN_GEN)  /* in this case it's just as efficient
                                 to use the full MVN */
    mvn_sample_anti(mmvn->mvn, retval1, retval2);
  else {
    Vector *xcomp1 = vec_new(mmvn->n), *xcomp2 = vec_new(mmvn->n);
    int d;
    assert(retval1->size == mmvn->d * mmvn->n &&
           retval2->size == mmvn->d * mmvn->n);
    for (d = 0; d < mmvn->d; d++) {
      mmvn->mvn->mu = mmvn->mu[d]; /* swap in the appropriate mean */
      mvn_sample_anti(mmvn->mvn, xcomp1, xcomp2);
      mmvn_project_up(mmvn, xcomp1, retval1, d);
      mmvn_project_up(mmvn, xcomp2, retval2, d);
    }
    vec_free(xcomp1);
    vec_free(xcomp2);
  }
}

void mmvn_map_std(multi_MVN *mmvn, Vector *retval) {
  if (mmvn->type != MVN_GEN) 
    mvn_map_std(mmvn->mvn, retval);
  else {
    Vector *xcomp = vec_new(mmvn->n);
    int d;
    assert(retval->size == mmvn->d * mmvn->n);
    for (d = 0; d < mmvn->d; d++) {
      mmvn_project_down(mmvn, retval, xcomp, d);
      mmvn->mvn->mu = mmvn->mu[d]; /* swap in the appropriate mean */
      mvn_map_std(mmvn->mvn, xcomp);
      mmvn_project_up(mmvn, xcomp, retval, d);
    }
  }
}

double mmvn_log_dens(multi_MVN *mmvn, Vector *x) {
  double retval = 0;
  if (mmvn->type != MVN_GEN)  /* in this case it's just as efficient
                                 to use the full MVN */
    retval = mvn_log_dens(mmvn->mvn, x);
  else {
    Vector *x_d = vec_new(mmvn->n);
    int d;
    assert(x->size == mmvn->d * mmvn->n);
    for (d = 0; d < mmvn->d; d++) {
      mmvn_project_down(mmvn, x, x_d, d);
      mmvn->mvn->mu = mmvn->mu[d]; /* swap in the appropriate mean */
      retval += mvn_log_dens(mmvn->mvn, x_d); /* total log density is
                                                 sum of component log
                                                 densities */
    }
    vec_free(x_d);
  }  
  return retval;
}

double mmvn_log_det(multi_MVN *mmvn) {
  if (mmvn->type != MVN_GEN)
    return mvn_log_det(mmvn->mvn);
  else 
    return mmvn->d * mvn_log_det(mmvn->mvn);
}

/* take a vector in the full n * d dimensional space and project it down to the
   dth dimension.  Essentially pulls out the coordinate axis of each
   of n points corresponding to dimension d */
void mmvn_project_down(multi_MVN *mmvn, Vector *x_full, Vector *x_d, int d) {
  int i;
  assert(x_d->size == mmvn->n &&
         x_full->size == mmvn->n * mmvn->d && mmvn->type == MVN_GEN);
  for (i = 0; i < mmvn->n; i++) 
    vec_set(x_d, i, vec_get(x_full, i*mmvn->d + d));
}

/* take a vector representing coordinates along dimension d and
   project them up to the full n*d dimensional space */
void mmvn_project_up(multi_MVN *mmvn, Vector *x_d, Vector *x_full, int d) {
  int i;
  assert(x_d->size == mmvn->n &&
         x_full->size == mmvn->n * mmvn->d && mmvn->type == MVN_GEN);
  for (i = 0; i < mmvn->n; i++) 
    vec_set(x_full, i*mmvn->d + d, vec_get(x_d, i));
}

/* obtain underlying standard mvn from a multi MVN.  Note that this is
   not guaranteed to return the same standard mvn as would
   mvn_rederive_std on the full-dimensional MVN, because of
   differences in eigenvector directionality.  Both values, however,
   will lie on the same contour of the density function. */
void mmvn_rederive_std(multi_MVN *mmvn, Vector *points, Vector *points_std) {
  if (mmvn->type != MVN_GEN) 
    mvn_rederive_std(mmvn->mvn, points, points_std);
  else {
    /* rederive the component standard mvns and project each one up */
    int d;
    Vector *points_d = vec_new(mmvn->n), *points_d_std = vec_new(mmvn->n);
    for (d = 0; d < mmvn->d; d++) {
      mmvn_project_down(mmvn, points, points_d, d);
      mmvn->mvn->mu = mmvn->mu[d];
      mvn_rederive_std(mmvn->mvn, points_d, points_d_std);
      mmvn_project_up(mmvn, points_d_std, points_std, d);
    }
    vec_free(points_d); vec_free(points_d_std);
  }
}

/* return trace of covariance for a multi MVN */
double mmvn_trace(multi_MVN *mmvn) {
  if (mmvn->type != MVN_GEN)
    return mvn_trace(mmvn->mvn);
  else
    return mmvn->d * mvn_trace(mmvn->mvn);
}

/* return trace of covariance for a multi MVN */
double mmvn_mu2(multi_MVN *mmvn) {
  if (mmvn->type != MVN_GEN)
    return mvn_mu2(mmvn->mvn);
  else {
    double retval = 0;
    int d;
    for (d = 0; d < mmvn->d; d++) {
      mmvn->mvn->mu = mmvn->mu[d];
      retval += mvn_mu2(mmvn->mvn);
    }
    return retval;
  }
}

/* save the mean parameters from a multi MVN in a vector for future
   retrieval */
void mmvn_save_mu(multi_MVN *mmvn, Vector *mu_saved) {
  if (mmvn->type != MVN_GEN)
    vec_copy(mu_saved, mmvn->mvn->mu);
  else {
    int d;
    assert(mu_saved->size == mmvn->d * mmvn->n);
    for (d = 0; d < mmvn->d; d++)
      mmvn_project_up(mmvn, mmvn->mu[d], mu_saved, d);
  }
}

void mmvn_print(multi_MVN *mmvn, FILE *F, unsigned int in_line,
                unsigned int verbose) {
  int j;
  Vector *mu_full = vec_new(mmvn->d * mmvn->n);
  mmvn_save_mu(mmvn, mu_full);
  if (!in_line)
    fprintf(F, "mu: ");
  for (j = 0; j < mu_full->size; j++)
    fprintf(F, "%f\t", vec_get(mu_full, j));
  if (!in_line)
    fprintf(F, "\n");
  if (verbose) {
    fprintf(F, "Component MVN:\n");
    mvn_print(mmvn->mvn, F);
  }
  vec_free(mu_full);
}

double mmvn_get_mu_el(multi_MVN *mmvn, int i) {
  if (mmvn->type != MVN_GEN)
    return vec_get(mmvn->mvn->mu, i);
  else {
    int d = i % mmvn->d, taxon = i / mmvn->d;
    return vec_get(mmvn->mu[d], taxon);   
  }
}

void mmvn_set_mu_el(multi_MVN *mmvn, int i, double val) {
  if (mmvn->type != MVN_GEN)
    vec_set(mmvn->mvn->mu, i, val);
  else {
    int d = i % mmvn->d, taxon = i / mmvn->d;
    vec_set(mmvn->mu[d], taxon, val);   
  }
}
