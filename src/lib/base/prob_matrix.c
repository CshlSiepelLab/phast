/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: prob_matrix.c,v 1.15 2008-11-12 02:07:59 acs Exp $ */

/* Matrices representing discrete bivariate probability distributions
   defined for non-negative integers.  General idea is element (x,y)
   of matrix M (x,y >= 0) represents p(x,y).  With long-tailed
   distributions, matrix dimensions are truncated at size x_max, y_max
   such that p(x,y) < epsilon for x >= x_max, y >= x_max, where
   epsilon is an input parameter. */

#include <prob_matrix.h>
#include <prob_vector.h>
#include <misc.h>

void pm_mean(Matrix *p, double *mean_x, double *mean_y) {
  int x, y;
  *mean_x = *mean_y = 0;
  for (x = 0; x < p->nrows; x++) {
    for (y = 0; y < p->ncols; y++) {
      *mean_x += x * p->data[x][y];
      *mean_y += y * p->data[x][y];
    }
  }
}

/* return marginal distribution for x, sum_y p(x, y) */
Vector *pm_marg_x(Matrix *p) {
  int x, y;
  Vector *marg = vec_new(p->nrows);
  vec_zero(marg);
  for (x = 0; x < p->nrows; x++)
    for (y = 0; y < p->ncols; y++)
      marg->data[x] += p->data[x][y];
  return marg;
}

/* return marginal distribution for y, sum_x p(x, y) */
Vector *pm_marg_y(Matrix *p) {
  int x, y;
  Vector *marg = vec_new(p->ncols);
  vec_zero(marg);
  for (x = 0; x < p->nrows; x++)
    for (y = 0; y < p->ncols; y++)
      marg->data[y] += p->data[x][y];
  return marg;
}

/* return marginal distribution for x + y, p(n) = sum_{x,y:x+y=n} p(x,y) */
Vector *pm_marg_tot(Matrix *p) {
  int x, y;
  Vector *marg = vec_new(p->ncols + p->nrows - 1);
  vec_zero(marg);
  for (x = 0; x < p->nrows; x++)
    for (y = 0; y < p->ncols; y++)
      marg->data[x+y] += p->data[x][y];
  return marg;
}

/* return conditional probability distribution of x given x+y */
Vector *pm_x_given_tot(Matrix *p, int tot) {
  int x;
  Vector *cond = vec_new(min(p->nrows, tot+1));
  vec_zero(cond);

  if (!(tot <= p->nrows + p->ncols - 2))
    die("ERROR pm_x_given_tot: tot=%i, p->nrows=%i, p->ncols=%i\n",
	tot, p->nrows, p->ncols);
  /* joint distribution has to allow for a total this large */

  for (x = max(0, tot - p->ncols + 1); x < p->nrows; x++) {
    if (x > tot) break;
    cond->data[x] = p->data[x][tot-x];
  }
  pv_normalize(cond);
  return cond;
}

/* return conditional probability distribution of y given x+y */
Vector *pm_y_given_tot(Matrix *p, int tot) {
  int y;
  Vector *cond = vec_new(min(p->ncols, tot+1));
  vec_zero(cond);

  if (!(tot <= p->nrows + p->ncols - 2))
    die("ERROR: pm_y_given_tot: tot=%i, p->nrows=%i, p->ncols=%i\n",
	tot, p->nrows, p->ncols);
  /* joint distribution has to allow for a total this large */

  for (y = max(0, tot - p->nrows + 1); y < p->ncols; y++) {
    if (y > tot) break;
    cond->data[y] = p->data[tot-y][y];
  }
  pv_normalize(cond);
  return cond;
}

/* normalize such that all elements sum to one */
void pm_normalize(Matrix *p) {
  int x, y;
  double sum = 0;
  for (x = 0; x < p->nrows; x++) 
    for (y = 0; y < p->ncols; y++) 
      sum += p->data[x][y];
  mat_scale(p, 1/sum);
}

/* convolve distribution n times */
Matrix *pm_convolve(Matrix *p, int n, double epsilon) {
  int i, j, k, x, y;
  Matrix *q_i, *q_i_1;
  double mean, var, max_nsd;
  int max_nrows = p->nrows * n, max_ncols = p->ncols * n;

  if (n <= 0)
    die("ERROR pm_convlve: n=%i\n", n);

  if (n == 1) 
    return mat_create_copy(p);

  if (n > 50) {
    /* use central limit theorem to limit size of matrix to
       keep track of.  Here work with marginal distributions */
    Vector *marg_x = pm_marg_x(p);
    Vector *marg_y = pm_marg_y(p);
    pv_stats(marg_x, &mean, &var);
    max_nsd = -inv_cum_norm(epsilon) + 1; 
    max_nrows = ceil(n * mean + max_nsd * sqrt(n * var)) + 1;
    pv_stats(marg_y, &mean, &var);
    max_ncols = ceil(n * mean + max_nsd * sqrt(n * var)) + 1;
    vec_free(marg_x);
    vec_free(marg_y);
  }

  q_i = mat_new(max_nrows, max_ncols);
  q_i_1 = mat_new(max_nrows, max_ncols);

  /* compute convolution recursively */
  mat_zero(q_i_1);
  for (x = 0; x < p->nrows; x++)
    for (y = 0; y < p->ncols; y++)
      q_i_1->data[x][y] = p->data[x][y];

  for (i = 1; i < n; i++) {
    mat_zero(q_i);
    for (x = 0; x < q_i->nrows; x++) {
      for (y = 0; y < q_i->ncols; y++) 
        for (j = max(0, x - p->nrows + 1); j <= x; j++) 
          for (k = max(0, y - p->ncols + 1); k <= y; k++) 
            q_i->data[x][y] += q_i_1->data[j][k] * p->data[x - j][y - k];
    }
    mat_copy(q_i_1, q_i);
  }

  mat_free(q_i_1);

  /* trim dimension before returning */
  max_nrows = max_ncols = -1;
  for (x = q_i->nrows - 1; max_nrows == -1 && x >= 0; x--) 
    for (y = 0; max_nrows == -1 && y < q_i->ncols; y++) 
      if (q_i->data[x][y] > epsilon) 
        max_nrows = x+1;      
  for (y = q_i->ncols - 1; max_ncols == -1 && y >= 0; y--) 
    for (x = 0; max_ncols == -1 && x < q_i->nrows; x++) 
      if (q_i->data[x][y] > epsilon) 
        max_ncols = y+1;
  mat_resize(q_i, max_nrows, max_ncols);

  pm_normalize(q_i);
  return q_i;
}

/* convolve distribution n times and keep all intermediate
   distributions.  Return value is an array q such that q[i] (1 <= i
   <= n) is the ith convolution of p (q[0] will be NULL) */
Matrix **pm_convolve_save(Matrix *p, int n, double epsilon) {
  int i, j, k, x, y;
  double mean, var, max_nsd;
  int max_nrows = p->nrows * n, max_ncols = p->ncols * n;
  Matrix **q = smalloc((n+1) * sizeof(void*));

  if (n <=0)
    die("ERROR pm_convolve_save: n=%i\n", n);

  q[0] = NULL;			/* placeholder */

  if (n == 1) {
    q[1] = mat_create_copy(p);
    return q;
  }

  if (n > 50) {
    /* use central limit theorem to limit size of matrix to keep track
       of.  Here work with marginal distributions. */
    Vector *marg_x = pm_marg_x(p);
    Vector *marg_y = pm_marg_y(p);
    pv_stats(marg_x, &mean, &var);
    max_nsd = -inv_cum_norm(epsilon) + 1; 
    max_nrows = max(1, ceil(n * mean + max_nsd * sqrt(n * var)));
    pv_stats(marg_y, &mean, &var);
    max_ncols = max(1, ceil(n * mean + max_nsd * sqrt(n * var)));
    vec_free(marg_x);
    vec_free(marg_y);
  }

  /* compute convolution recursively */
  q[1] = mat_new(max_nrows, max_ncols);
  mat_zero(q[1]);
  for (x = 0; x < p->nrows; x++)
    for (y = 0; y < p->ncols; y++)
      q[1]->data[x][y] = p->data[x][y];

  for (i = 2; i <= n; i++) {
    q[i] = mat_new(max_nrows, max_ncols);
    mat_zero(q[i]);
    for (x = 0; x < q[i]->nrows; x++) 
      for (y = 0; y < q[i]->ncols; y++) 
        for (j = max(0, x - p->nrows + 1); j <= x; j++) 
          for (k = max(0, y - p->ncols + 1); k <= y; k++) 
            q[i]->data[x][y] += q[i-1]->data[j][k] * p->data[x - j][y - k];
  }

  /* trim dimension before returning */
  for (i = 1; i <= n; i++) {
    max_nrows = max_ncols = -1;
    for (x = q[i]->nrows - 1; max_nrows == -1 && x >= 0; x--) 
      for (y = 0; max_nrows == -1 && y < q[i]->ncols; y++) 
        if (q[i]->data[x][y] > epsilon) 
          max_nrows = x+1;      
    for (y = q[i]->ncols - 1; max_ncols == -1 && y >= 0; y--) 
      for (x = 0; max_ncols == -1 && x < q[i]->nrows; x++) 
        if (q[i]->data[x][y] > epsilon) 
          max_ncols = y+1;
    mat_resize(q[i], max_nrows, max_ncols);
    pm_normalize(q[i]);
  }

  return q;
}

/* take convolution of a set of probability matrices.  If counts is
   NULL, then each distrib is assumed to have multiplicity 1 */
Matrix *pm_convolve_many(Matrix **p, int *counts, int n, double epsilon) {
  int i, j, k, l, x, y, max_nrows, max_ncols, count, tot_count = 0,
    this_max_nrows, this_max_ncols;
  Matrix *q_i, *q_i_1;
  double max_nsd;

  max_nrows = max_ncols = 0; 
  for (i = 0; i < n; i++) {
    count = (counts == NULL ? 1 : counts[i]);
    tot_count += count;
    max_nrows += p[i]->nrows;
    max_ncols += p[i]->ncols;
  }

  if (n == 1 && (counts == NULL || counts[0] == 1))
    /* no convolution necessary */
    return mat_create_copy(p[0]);

  if (tot_count > 50) {
    /* as above, use central limit theorem to limit size of matrix to
       keep track of */
    double tot_mean_x = 0, tot_var_x = 0, tot_mean_y = 0, tot_var_y = 0,
      mean, var;
    for (i = 0; i < n; i++) {
      Vector *marg_x = pm_marg_x(p[i]);
      Vector *marg_y = pm_marg_y(p[i]);
      pv_stats(marg_x, &mean, &var);
      count = (counts == NULL ? 1 : counts[i]);
      tot_mean_x += mean * count;
      tot_var_x += var * count;
      pv_stats(marg_y, &mean, &var);
      tot_mean_y += mean * count;
      tot_var_y += var * count;
      vec_free(marg_x); vec_free(marg_y);
    }

    max_nsd = -inv_cum_norm(epsilon) + 1; 
    max_nrows = ceil(tot_mean_x + max_nsd * sqrt(tot_var_x)) + 1;
    max_ncols = ceil(tot_mean_y + max_nsd * sqrt(tot_var_y)) + 1;
  }

  q_i = mat_new(max_nrows, max_ncols);
  q_i_1 = mat_new(max_nrows, max_ncols);

  /* compute convolution recursively */
  mat_zero(q_i_1);
  this_max_nrows = min(p[0]->nrows, max_nrows);
  this_max_ncols = min(p[0]->ncols, max_ncols);
  for (x = 0; x < this_max_nrows; x++)
    for (y = 0; y < this_max_ncols; y++)
      q_i_1->data[x][y] = p[0]->data[x][y];
 
  this_max_nrows = p[0]->nrows;
  this_max_ncols = p[0]->ncols;
  for (i = 0; i < n; i++) {
    count = (counts == NULL ? 1 : counts[i]);
    if (i == 0) count--; /* initialization takes care of first one */
    this_max_nrows = min(max_nrows, this_max_nrows + p[i]->nrows);
    this_max_ncols = min(max_ncols, this_max_ncols + p[i]->ncols);
    for (l = 0; l < count; l++) {
      mat_zero(q_i);
      for (x = 0; x < this_max_nrows; x++) {
        for (y = 0; y < this_max_ncols; y++) 
          for (j = max(0, x - p[i]->nrows + 1); j <= x; j++) 
            for (k = max(0, y - p[i]->ncols + 1); k <= y; k++) 
              q_i->data[x][y] += q_i_1->data[j][k] * p[i]->data[x - j][y - k];
      }
      mat_copy(q_i_1, q_i);
    }
  }

  mat_free(q_i_1);

  /* trim dimension before returning */
  max_nrows = max_ncols = -1;
  for (x = q_i->nrows - 1; max_nrows == -1 && x >= 0; x--) 
    for (y = 0; max_nrows == -1 && y < q_i->ncols; y++) 
      if (q_i->data[x][y] > epsilon) 
        max_nrows = x+1;      
  if (max_nrows == -1) max_nrows = q_i->nrows;
  for (y = q_i->ncols - 1; max_ncols == -1 && y >= 0; y--) 
    for (x = 0; max_ncols == -1 && x < q_i->nrows; x++) 
      if (q_i->data[x][y] > epsilon) 
        max_ncols = y+1;
  if (max_ncols == -1) max_ncols = q_i->ncols;
  mat_resize(q_i, max_nrows, max_ncols);

  pm_normalize(q_i);
  return q_i;
}

/* take convolution of a set of probability matrices, avoiding some
   overhead of function above; does not take counts, does not
   normalize, does not trim dimension, allows max size to be
   specified */
Matrix *pm_convolve_many_fast(Matrix **p, int n, int max_nrows, int max_ncols) {
  int i, j, k, x, y, this_max_nrows, this_max_ncols;
  Matrix *q_i, *q_i_1;

  if (n == 1)
    /* no convolution necessary */
    return mat_create_copy(p[0]);

  q_i = mat_new(max_nrows, max_ncols);
  q_i_1 = mat_new(max_nrows, max_ncols);

  /* compute convolution recursively */
  mat_zero(q_i_1);
  this_max_nrows = min(p[0]->nrows, max_nrows);
  this_max_ncols = min(p[0]->ncols, max_ncols);
  for (x = 0; x < this_max_nrows; x++)
    for (y = 0; y < this_max_ncols; y++)
      q_i_1->data[x][y] = p[0]->data[x][y];
 
  this_max_nrows = p[0]->nrows;
  this_max_ncols = p[0]->ncols;
  for (i = 1; i < n; i++) {
    this_max_nrows = min(max_nrows, this_max_nrows + p[i]->nrows);
    this_max_ncols = min(max_ncols, this_max_ncols + p[i]->ncols);
    mat_zero(q_i);
    for (x = 0; x < this_max_nrows; x++) {
      for (y = 0; y < this_max_ncols; y++) 
        for (j = max(0, x - p[i]->nrows + 1); j <= x; j++) 
          for (k = max(0, y - p[i]->ncols + 1); k <= y; k++) 
            q_i->data[x][y] += q_i_1->data[j][k] * p[i]->data[x - j][y - k];
    }
    mat_copy(q_i_1, q_i);
  }

  mat_free(q_i_1);
  return q_i;
}

/* convolve distribution n times, using a faster algorithm than the
   ones above; time is proportional to log(n) rather than n */
Matrix *pm_convolve_fast(Matrix *p, int n, double epsilon) {
  int i, j, checksum;
  int logn = log2_int(n);
  Matrix *pow_p[64], *pows[64];
  Matrix *retval;
  double mean, var, max_nsd;
  int max_nrows = p->nrows * n, max_ncols = p->ncols * n;

  if (n == 1)
    return mat_create_copy(p);

  if (n > 50) {
    /* use central limit theorem to limit size of matrix to
       keep track of.  Here work with marginal distributions */
    Vector *marg_x = pm_marg_x(p);
    Vector *marg_y = pm_marg_y(p);
    pv_stats(marg_x, &mean, &var);
    max_nsd = -inv_cum_norm(epsilon) + 1; 
    max_nrows = ceil(n * mean + max_nsd * sqrt(n * var)) + 1;
    pv_stats(marg_y, &mean, &var);
    max_ncols = ceil(n * mean + max_nsd * sqrt(n * var)) + 1;
    vec_free(marg_x);
    vec_free(marg_y);
  }

  /* Let p^i be the ith convolution of p (i.e., p o p o ... o p, i
     times, where 'o' is the convolution operator).  We use the fact
     that p^i o p^j = p^(i+j) to compose p^n from p^1, p^2, p^4, p^8,
     ... in log2(k) steps */

  /* compute "powers" of p */
  pow_p[0] = p;
  for (i = 1; i <= logn; i++) 
    pow_p[i] = pm_convolve(pow_p[i-1], 2, epsilon);

  /* now combine powers to get desired convolution */
  j = checksum = 0;
  for (i = 0; i <= logn; i++) {
    unsigned bit_i = (n >> i) & 1;
    if (bit_i) {
      pows[j++] = pow_p[i];
      checksum += int_pow(2, i);
    }
  }
  if (n != checksum)
    die("ERROR pm_convolve_fast: n (%i) != checksum (%i)\n", n, checksum);

  retval = pm_convolve_many_fast(pows, j, max_nrows, max_ncols);

  for (i = 1; i <= logn; i++) 
    mat_free(pow_p[i]);

  return retval;
}

/* compute means, variances, and covariance */
void pm_stats(Matrix *p, double *mean_x, double *mean_y, double *var_x, 
              double *var_y, double *covar) {
  int x, y;
  *mean_x = *mean_y = *var_x = *var_y = *covar = 0;
  for (x = 0; x < p->nrows; x++) {
    for (y = 0; y < p->ncols; y++) {
      *mean_x += x * p->data[x][y];
      *mean_y += y * p->data[x][y];
      *var_x += x * x * p->data[x][y]; 
      *var_y += y * y * p->data[x][y]; 
      *covar += x * y * p->data[x][y];
    }
  }
  *var_x -= (*mean_x * *mean_x);
  *var_y -= (*mean_y * *mean_y);
  *covar -= (*mean_x * *mean_y);
}

/* version of pm_x_given_tot that assumes bivariate normal with given
   means, standard deviations, and correlation coefficient.  Only the
   region of the bivariate normal with x,y >= 0 is considered.  For
   use in central limit theorem approximations. */
Vector *pm_x_given_tot_bvn(int tot, double mu_x, double mu_y, 
                           double sigma_x, double sigma_y, double rho) {
  int x;
  Vector *cond = vec_new(tot+1); 
  for (x = 0; x <= tot; x++) 
    cond->data[x] = bvn_p(x, tot-x, mu_x, mu_y, sigma_x, sigma_y, rho);
  pv_normalize(cond);
  return cond;  
}

/* version of pm_y_given_tot that assumes bivariate normal with given
   means, standard deviations, and correlation coefficient.  Only the
   region of the bivariate normal with x,y >= 0 is considered.  For
   use in central limit theorem approximations. */
Vector *pm_y_given_tot_bvn(int tot, double mu_x, double mu_y, 
                           double sigma_x, double sigma_y, double rho) {
  int y;
  Vector *cond = vec_new(tot+1); 
  for (y = 0; y <= tot; y++) 
    cond->data[y] = bvn_p(tot-y, y, mu_x, mu_y, sigma_x, sigma_y, rho);
  pv_normalize(cond);
  return cond;  
}

/* version of pm_x_given_tot that assumes independence of x and y, and
   computes the desired conditional distribution from their
   marginals */
Vector *pm_x_given_tot_indep(int tot, Vector *marg_x, Vector *marg_y) {
  int x;
  Vector *cond = vec_new(min(marg_x->size, tot+1));
  vec_zero(cond);
  for (x = max(0, tot - marg_y->size + 1); x < marg_x->size; x++) {
    if (x > tot) break;
    cond->data[x] = marg_x->data[x] * marg_y->data[tot-x];
  }
  pv_normalize(cond);
  return cond;  
}

/* version of pm_y_given_tot that assumes independence of x and y, and
   computes the desired conditional distribution from their
   marginals */
Vector *pm_y_given_tot_indep(int tot, Vector *marg_x, Vector *marg_y) {
  int y;
  Vector *cond = vec_new(min(marg_y->size, tot+1));
  vec_zero(cond);
  for (y = max(0, tot - marg_x->size + 1); y < marg_y->size; y++) {
    if (y > tot) break;
    cond->data[y] = marg_x->data[tot-y] * marg_y->data[y];
  }
  pv_normalize(cond);
  return cond;  
}
