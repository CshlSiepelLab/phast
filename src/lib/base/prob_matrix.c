/* $Id: prob_matrix.c,v 1.3 2005-08-12 22:47:43 acs Exp $ 
   Written by Adam Siepel, 2005
   Copyright 2005, Adam Siepel, University of California 
*/

/* Matrices representing discrete bivariate probability distributions
   defined for non-negative integers.  General idea is element (x,y)
   of matrix M (x,y >= 0) represents p(x,y).  With long-tailed
   distributions, matrix dimensions are truncated at size x_max, y_max
   such that p(x,y) < PM_EPS for x >= x_max, y >= x_max. */

#include <prob_matrix.h>
#include <prob_vector.h>
#include <misc.h>
#include <assert.h>

void pm_mean(Matrix *p, double *mean_x, double *mean_y) {
  int x, y;
  *mean_x = *mean_y = 0;
  for (x = 0; x < p->nrows; x++) {
    for (y = 0; y < p->ncols; y++) {
      *mean_x *= x * p->data[x][y];
      *mean_y *= y * p->data[x][y];
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
  Vector *cond = vec_new(p->nrows);
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
  Vector *cond = vec_new(p->ncols);
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
Matrix *pm_convolve(Matrix *p, int n) {
  int i, j, k, x, y;
  Matrix *q_i, *q_i_1;
  double mean, var;
  int max_nrows = p->nrows * n, max_ncols = p->ncols * n;

  assert(n > 0);

  if (n == 1) 
    return mat_create_copy(p);

  if (n > 50) {
    /* use central limit theorem to limit size of matrix to
       keep track of.  Here work with marginal distributions */
    Vector *marg_x = pm_marg_x(p);
    Vector *marg_y = pm_marg_y(p);
    pv_stats(marg_x, &mean, &var);
    max_nrows = ceil(n * mean + 6 * sqrt(n * var));
    pv_stats(marg_y, &mean, &var);
    max_ncols = ceil(n * mean + 6 * sqrt(n * var));
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
      if (q_i->data[x][y] > PM_EPS) 
        max_nrows = x+1;      
  for (y = q_i->ncols - 1; max_ncols == -1 && y >= 0; y--) 
    for (x = 0; max_ncols == -1 && x < q_i->nrows; x++) 
      if (q_i->data[x][y] > PM_EPS) 
        max_ncols = y+1;
  mat_resize(q_i, max_nrows, max_ncols);

  pm_normalize(q_i);
  return q_i;
}

/* convolve distribution n times and keep all intermediate
   distributions.  Return value is an array q such that q[i] (1 <= i
   <= n) is the ith convolution of p (q[0] will be NULL) */
Matrix **pm_convolve_save(Matrix *p, int n) {
  int i, j, k, x, y;
  double mean, var;
  int max_nrows = p->nrows * n, max_ncols = p->ncols * n;
  Matrix **q = smalloc((n+1) * sizeof(void*));

  assert(n > 0);

  q[0] = NULL;			/* placeholder */

  if (n == 1) {
    q[1] = mat_create_copy(p);
    return q;
  }

  if (n > 50) {
    /* use central limit theorem to limit size of matrix to
       keep track of.  Here work with marginal distributions */
    Vector *marg_x = pm_marg_x(p);
    Vector *marg_y = pm_marg_y(p);
    pv_stats(marg_x, &mean, &var);
    max_nrows = ceil(n * mean + 6 * sqrt(n * var));
    pv_stats(marg_y, &mean, &var);
    max_ncols = ceil(n * mean + 6 * sqrt(n * var));
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
        if (q[i]->data[x][y] > PM_EPS) 
          max_nrows = x+1;      
    for (y = q[i]->ncols - 1; max_ncols == -1 && y >= 0; y--) 
      for (x = 0; max_ncols == -1 && x < q[i]->nrows; x++) 
        if (q[i]->data[x][y] > PM_EPS) 
          max_ncols = y+1;
    mat_resize(q[i], max_nrows, max_ncols);
    pm_normalize(q[i]);
  }

  return q;
}

/* take convolution of a set of probability matrices.  If counts is
   NULL, then each distrib is assumed to have multiplicity 1 */
Matrix *pm_convolve_many(Matrix **p, int *counts, int n) {
  int i, j, k, l, x, y, max_nrows, max_ncols, count, tot_count = 0;
  Matrix *q_i, *q_i_1;
  double tot_mean_x = 0, tot_var_x = 0, tot_mean_y = 0, tot_var_y = 0,
    mean, var;
  Vector *marg_x, *marg_y;

  max_nrows = max_ncols = 0; 
  for (i = 0; i < n; i++) {
    count = (counts == NULL ? 1 : counts[i]);
    tot_count += count;
    max_nrows += p[i]->nrows;
    max_ncols += p[i]->ncols;
    marg_x = pm_marg_x(p[i]);
    marg_y = pm_marg_y(p[i]);
    pv_stats(marg_x, &mean, &var);
    tot_mean_x += mean * count;
    tot_var_x += var * count;
    pv_stats(marg_y, &mean, &var);
    tot_mean_y += mean * count;
    tot_var_y += var * count;
    vec_free(marg_x); vec_free(marg_y);
  }

  if (n == 1 && (counts == NULL || counts[0] == 1))
    /* no convolution necessary */
    return mat_create_copy(p[0]);

  if (n > 50) {
    /* as above, use central limit theorem to limit size of matrix to
       keep track of */
    max_nrows = ceil(tot_mean_x + 6 * sqrt(tot_var_x));
    max_ncols = ceil(tot_mean_y + 6 * sqrt(tot_var_y));
  }

  q_i = mat_new(max_nrows, max_ncols);
  q_i_1 = mat_new(max_nrows, max_ncols);

  /* compute convolution recursively */
  mat_zero(q_i_1);
  for (x = 0; x < p[0]->nrows; x++)
    for (y = 0; y < p[0]->ncols; y++)
      q_i_1->data[x][y] = p[0]->data[x][y];
 
  for (i = 0; i < n; i++) {
    count = (counts == NULL ? 1 : counts[i]);
    if (i == 0) count--; /* initialization takes care of first one */
    for (l = 0; l < count; l++) {
      mat_zero(q_i);
      for (x = 0; x < q_i->nrows; x++) {
        for (y = 0; y < q_i->ncols; y++) 
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
      if (q_i->data[x][y] > 1e-10) 
        max_nrows = x+1;      
  for (y = q_i->ncols - 1; max_ncols == -1 && y >= 0; y--) 
    for (x = 0; max_ncols == -1 && x < q_i->nrows; x++) 
      if (q_i->data[x][y] > 1e-10) 
        max_ncols = y+1;
  mat_resize(q_i, max_nrows, max_ncols);

  pm_normalize(q_i);
  return q_i;
}
