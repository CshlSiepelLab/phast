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

/* wrapper functions to take advantage of MVNs that represent points
   in d dimensions, such that the coordinate axes are independent, but
   along each axis the points share the same covariance structure.
   This representation is used by the neighbor-joining/variational
   code. */

/* create a new multi-MVN of desired type from scratch */
multi_MVN *mmvn_new(int n, int d, enum covar_type type) {
  MVN *mvn = NULL;
  if (type == DIAG) 
    mvn = mvn_new(n*d, NULL, NULL);
  else
    mvn = mvn_new(n, NULL, NULL);
  return mmvn_new_from_mvn(mvn, d, type);
}

/* create a multi-MVN from an MVN */
multi_MVN *mmvn_new_from_mvn(MVN *mvn, int d, enum covar_type type) {
  multi_MVN *retval = smalloc(sizeof(multi_MVN));
  retval->d = d;
  retval->ntips = (type == DIST ? mvn->dim : mvn->dim / d);
  retval->mvn = mvn;
  retval->type = type;

  if (type == DIST) {
    int i;
    /* clear the orig mu vector; will handle it just by pointer swapping */
    vec_free(retval->mvn->mu);
    retval->mu = smalloc(d * sizeof(Vector*));
    for (i = 0; i < d; i++)
      retval->mu[i] = vec_new(retval->ntips);
  }
  else
    retval->mu = NULL;
  
  return retval;
}

void mmvn_sample(multi_MVN *mmvn, Vector *retval) {
  if (mmvn->type == DIAG)  /* in this case it's just as efficient to use the full MVN */
    mvn_sample(mmvn->mvn, retval);
  else {
    Vector *xcomp = vec_new(mmvn->ntips);
    int d;
    assert(retval->size == mmvn->d * mmvn->ntips);
    for (d = 0; d < mmvn->d; d++) {
      mmvn->mvn->mu = mmvn->mu[d]; /* swap in the appropriate mean */
      mvn_sample(mmvn->mvn, xcomp);
      mmvn_project_up(mmvn, xcomp, retval, d);
    }
  }
}

double mmvn_log_dens(multi_MVN *mmvn, Vector *x) {
  double retval = 0;
  if (mmvn->type == DIAG)  /* in this case it's just as efficient to use the full MVN */
    retval = mvn_log_dens(mmvn->mvn, x);
  else {
    Vector *x_d = vec_new(mmvn->ntips);
    int d;
    assert(x->size == mmvn->d * mmvn->ntips);
    for (d = 0; d < mmvn->d; d++) {
      mmvn_project_down(mmvn, x, x_d, d);
      mmvn->mvn->mu = mmvn->mu[d]; /* swap in the appropriate mean */
      retval += mvn_log_dens(mmvn->mvn, x_d); /* total log density is sum of component log densities */
    }
    vec_free(x_d);
  }  
  return retval;
}

double mmvn_log_det(multi_MVN *mmvn) {
  if (mmvn->type == DIAG)
    return mvn_log_det(mmvn->mvn);
  else 
    return mmvn->d * mvn_log_det(mmvn->mvn);
}

/* take a vector in the full n * d dimensional space and project it down to the
   dth dimension.  Essentially pulls out the coordinate axis of each
   of n points corresponding to dimension d */
void mmvn_project_down(multi_MVN *mmvn, Vector *x_full, Vector *x_d, int d) {
  int i;
  assert(x_d->size == mmvn->ntips && x_full->size == mmvn->ntips * mmvn->d && mmvn->type == DIST);
  for (i = 0; i < mmvn->ntips; i++) 
    vec_set(x_d, i, vec_get(x_full, i*mmvn->d + d));
}

/* take a vector in representing coordinates along dimension d and
   project them up to the full n*d dimensional space */
void mmvn_project_up(multi_MVN *mmvn, Vector *x_d, Vector *x_full, int d) {
  int i;
  assert(x_d->size == mmvn->ntips && x_full->size == mmvn->ntips * mmvn->d && mmvn->type == DIST);
  for (i = 0; i < mmvn->ntips; i++) 
    vec_set(x_full, i*mmvn->d + d, vec_get(x_d, i));
}

/* obtain underlying standard mvn from a multi MVN */
void mmvn_rederive_std(multi_MVN *mmvn, Vector *points, Vector *points_std) {
  if (mmvn->type == DIAG) 
    mvn_rederive_std(mmvn->mvn, points, points_std);
  else {
    /* rederive the component standard mvns and project each one up */
    int d;
    Vector *points_d = vec_new(mmvn->ntips), *points_d_std = vec_new(mmvn->ntips);
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
  if (mmvn->type == DIAG)
    return mvn_trace(mmvn->mvn);
  else
    return mmvn->d * mvn_trace(mmvn->mvn);
}

/* return trace of covariance for a multi MVN */
double mmvn_mu2(multi_MVN *mmvn) {
  if (mmvn->type == DIAG)
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
  if (mmvn->type == DIAG)
    vec_copy(mu_saved, mmvn->mvn->mu);
  else {
    int d;
    assert(mu_saved->size == mmvn->d * mmvn->ntips);
    for (d = 0; d < mmvn->d; d++)
      mmvn_project_up(mmvn, mmvn->mu[d], mu_saved, d);
  }
}

/* restore the mean parameters from a saved version */
void mmvn_restore_mu(multi_MVN *mmvn, Vector *mu_saved) {
  if (mmvn->type == DIAG)
    vec_copy(mmvn->mvn->mu, mu_saved);
  else {
    int d;
    assert(mu_saved->size == mmvn->d * mmvn->ntips);
    for (d = 0; d < mmvn->d; d++)
      mmvn_project_down(mmvn, mu_saved, mmvn->mu[d], d);
  }
}

void mmvn_print(multi_MVN *mmvn, FILE *F, unsigned int in_line,
                unsigned int verbose) {
  int j;
  Vector *mu_full = vec_new(mmvn->d * mmvn->ntips);
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
  if (mmvn->type == DIAG)
    return vec_get(mmvn->mvn->mu, i);
  else {
    int d = i % mmvn->d, taxon = i / mmvn->d;
    return vec_get(mmvn->mu[d], taxon);   /* CHECK */
  }
}

void mmvn_set_mu_el(multi_MVN *mmvn, int i, double val) {
  if (mmvn->type == DIAG)
    vec_set(mmvn->mvn->mu, i, val);
  else {
    int d = i % mmvn->d, taxon = i / mmvn->d;
    vec_set(mmvn->mu[d], taxon, val);   /* CHECK */
  }
}
