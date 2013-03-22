/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: prob_vector.c,v 1.14 2008-11-12 02:07:59 acs Exp $ */

/* Vectors representing discrete probability distributions over
   non-negative integers.  General idea is element x of vector v (x >=
   0) represents p(x).  With long-tailed distributions (e.g.,
   Poisson), vectors are truncated at size x_max such that p(y) <
   epsilon for y >= x_max, where epsilon is an input parameter. */

#include <prob_vector.h>
#include <misc.h>

/* compute mean and variance */
void pv_stats(Vector *p, double *mean, double *var) {  
  int x;
  *mean = 0;
  *var = 0;
  for (x = 0; x < p->size; x++) {
    *mean += x * p->data[x];
    *var += x * x * p->data[x];
  }
  *var -= (*mean * *mean);
}

/* compute min and max of specified confidence interval.  Size of
   confidence interval (given by 'size') should be between 0 and 1.  A
   central interval is assumed (equal probability in tails).  */
void pv_confidence_interval(Vector *p, double size, int *interval_min, 
                            int *interval_max) {
  int x;
  double tail_size = (1-size)/2;
  double cum = 0;
  if (p->size <= 0)
    die("ERROR pv_confidence_interval got p->size=%i\n", p->size);
  if (!(size > 0 && size < 1))
    die("ERROR pv_confidence_interval got size=%f\n", size);
  
  for (x = 0; x < p->size; x++) {
    cum += p->data[x];
    if (cum >= tail_size) {     /* total prob up to and including x */
      if (cum == tail_size) *interval_min = x + 1;
      else *interval_min = x;
      /* round towards tail -- ensures total probability of interval
         is at least size */
      break;
    }
  }
  cum = 0;
  for (x = p->size-1; x >= 0; x--) {
    cum += p->data[x];
    if (cum >= tail_size) {
      if (cum == tail_size) *interval_max = x - 1;
      else *interval_max = x;
      break;
    }
  }

  /* correct for possible overflow */
  if (*interval_min >= p->size) *interval_min = p->size - 1;
  if (*interval_max < 0) *interval_max = 0;
}

/* Compute quantiles 0.00, 0.01, 0.02, ..., 1.00.  Returns an array of
   length 101, such that element x contains the x/100th quantile.  */
int* pv_quantiles(Vector *p) {
  int *retval = smalloc(101 * sizeof(int));
  int qidx = 0, x;
  double cum = 0;
  for (x = 0; x < p->size; x++) {
    cum += p->data[x];
    while (cum >= 1.0*qidx/100)
      retval[qidx++] = x;
  }
  if (qidx < 101) retval[100] = p->size-1; /* rounding error */
  return retval;
}

/* return one-sided p-value: p(x <= x_0) if side == LOWER, p(x >= x_0)
   if side == UPPER.  If side == TWOTAIL, heuristically returns 2 *
   min(p(x <= x_0, p(x >= x_0)) (for drawbacks of this approach and
   discussion see Dunne et al., The Statistician, 1996) */  
double pv_p_value(Vector *distrib, double x_0, p_val_type side) {
  double lretval = 0, hretval = 0;
  int x;
  if (side != UPPER)
    for (x = 0; x < distrib->size && x <= x_0; x++) 
      lretval += distrib->data[x];
  if (side != LOWER)
    for (x = distrib->size-1; x >= 0 && x >= x_0; x--) 
      hretval += distrib->data[x];

  if (side == LOWER) return lretval;
  else if (side == UPPER) return hretval;
  return 2 * (min(lretval, hretval));
}

/* compute one-sided p-values for array of values.  Like pv_p_value, but
   saves time by computing CDF and using for all pvals */
void pv_p_values(Vector *distrib, double *x_0, int n, double *pvals,
                 p_val_type side) {
  Vector *lcdf = NULL, *hcdf = NULL;
  int i;

  if (side != UPPER)
    lcdf = pv_cdf(distrib, LOWER); /* CDF for lower tail */
  if (side != LOWER)
    hcdf = pv_cdf(distrib, UPPER); /* CDF for upper tail */
  
  /* look up tail probabilities from CDF */
  for (i = 0; i < n; i++) {
    if (side == LOWER)
      pvals[i] = lcdf->data[(int)floor(x_0[i])];
    else if (side == UPPER)
      pvals[i] = hcdf->data[(int)ceil(x_0[i])];
    else                        /* side == TWOTAIL */
      pvals[i] = 2*(min(lcdf->data[(int)floor(x_0[i])], 
                        hcdf->data[(int)ceil(x_0[i])]));
  }

  if (lcdf != NULL) vec_free(lcdf);
  if (hcdf != NULL) vec_free(hcdf);
}

/* normalize distribution */
void pv_normalize(Vector *p) {
  int x;
  double sum = 0;
  for (x = 0; x < p->size; x++) {
    if (p->data[x] < 0) die("ERROR pv_normalize got pv->data[%i]=%f\n",
			    x, p->data[x]);
    sum += p->data[x];
  }
  vec_scale(p, 1/sum);
}

/* convolve distribution n times */
Vector *pv_convolve(Vector *p, int n, double epsilon) {
  int i, j, x;
  Vector *q_i, *q_i_1;
  double mean, var, max_nsd;
  int max_x = p->size * n;

  if (n == 1)
    return vec_create_copy(p);

  if (n > 50) {
    /* use central limit theorem to limit size of vector to keep track
       of; convolution should be approx normal with mean n * mean of p
       and variance n * variance of p.  Note that the CLT is better
       near the mean than on the tails, but this seems okay for
       purposes of bounding. */
    pv_stats(p, &mean, &var);
    max_nsd = -inv_cum_norm(epsilon) + 1; 
    max_x = max((int)ceil(n * mean + max_nsd * sqrt(n * var)), p->size);
  }

  q_i = vec_new(max_x);
  q_i_1 = vec_new(max_x);

  /* compute convolution recursively */
  vec_zero(q_i_1);
  for (x = 0; x < p->size; x++)
    q_i_1->data[x] = p->data[x];

  for (i = 1; i < n; i++) {
    vec_zero(q_i);
    for (x = 0; x < q_i->size; x++) {
      for (j = max(0, x - p->size + 1); j <= x; j++) 
        q_i->data[x] += q_i_1->data[j] * p->data[x - j];
    }
    if (i < n - 1) vec_copy(q_i_1, q_i);
  }

  vec_free(q_i_1);

  /* trim very small values off tail before returning */
  for (x = q_i->size - 1; x >= 0; x--) {
    if (q_i->data[x] > epsilon) {
      q_i->size = x+1;
      break;
    }
  }

  pv_normalize(q_i);
  return q_i;
}

/* convolve distribution n times and keep all intermediate
   distributions.  Return value is an array q such that q[i] (1 <= i <=
   n) is the ith convolution of p (q[0] will be NULL) */
Vector **pv_convolve_save(Vector *p, int n, double epsilon) {
  int i, j, x;
  double mean, var, max_nsd;
  int max_x = p->size * n, newsize;
  Vector **q = smalloc((n+1) * sizeof(void*));

  q[0] = NULL;			/* placeholder */

  if (n == 1) {
    q[1] = vec_create_copy(p);
    return q;
  }

  if (n > 50) {
    /* use central limit theorem to limit size of vector to keep track
       of; convolution should be approx normal with mean n * mean
       of p and variance n * variance of p.  We'll go to 6
       standard deviations beyond the mean, which should ensure any
       omitted value has prob on the order of 1e-10 or less.  Note
       that the CLT is better near the mean than on the tails, but
       this seems okay for purposes of bounding. */
    pv_stats(p, &mean, &var);
    max_nsd = -inv_cum_norm(epsilon) + 1; 
    max_x = (int)ceil(n * mean + max_nsd * sqrt(n * var));
  }

  /* compute convolution recursively */
  q[1] = vec_new(max_x);
  vec_zero(q[1]);
  for (x = 0; x < p->size; x++)
    q[1]->data[x] = p->data[x];

  for (i = 2; i <= n; i++) {
    q[i] = vec_new(max_x);
    vec_zero(q[i]);
    for (x = 0; x < q[i]->size; x++) 
      for (j = max(0, x - p->size + 1); j <= x; j++) 
        q[i]->data[x] += q[i-1]->data[j] * p->data[x - j];
  }

  /* trim very small values off tail before returning */
  for (i = 1; i <= n; i++) {
    newsize = -1;
    for (x = q[i]->size - 1; newsize == -1 && x >= 0; x--) 
      if (q[i]->data[x] > epsilon) 
        newsize = x+1;
    q[i]->size = newsize;       /* maybe should realloc? */
    pv_normalize(q[i]);
  }

  return q;
}

/* take convolution of a set of probability vectors.  If counts is
   NULL, then each distrib is assumed to have multiplicity 1 */
Vector *pv_convolve_many(Vector **p, int *counts, int n, double epsilon) {
  int i, j, k, x, max_x = 0, tot_count = 0, count, thismax;
  Vector *q_i, *q_i_1;
  double mean, var, max_nsd;

  for (i = 0; i < n; i++) {
    count = (counts == NULL ? 1 : counts[i]);
    tot_count += count;
    max_x += p[i]->size;
  }

  if (n == 1 && (counts == NULL || counts[0] == 1)) 
    /* no convolution necessary */
    return vec_create_copy(p[0]);

  if (tot_count > 50) {
    /* use (Lyapunov's or Lindeberg's) central limit theorem to reduce
       size of vector to keep track of; convolution should be approx
       normal with mean tot_mean and variance tot_var.  Note that the
       CLT is better near the mean than on the tails, but this seems
       okay here.  We're also implicitly assuming some regularity
       conditions (see http://en.wikipedia.org/wiki/Central_limit_theorem). */
    double tot_mean = 0, tot_var = 0;
    for (i = 0; i < n; i++) {
      pv_stats(p[i], &mean, &var);
      count = (counts == NULL ? 1 : counts[i]);
      tot_mean += mean * count;
      tot_var += var * count;
    }
    max_nsd = -inv_cum_norm(epsilon) + 1; 
    max_x = (int)ceil(tot_mean + max_nsd * sqrt(tot_var));
  }

  q_i = vec_new(max_x);
  q_i_1 = vec_new(max_x);

  /* compute convolution recursively */
  vec_zero(q_i_1);
  thismax = min(p[0]->size, max_x);
  for (x = 0; x < thismax; x++)
    q_i_1->data[x] = p[0]->data[x];

  for (i = 0; i < n; i++) {
    count = (counts == NULL ? 1 : counts[i]);
    if (i == 0) count--; /* initialization takes care of first one */
    thismax = min(max_x, thismax + p[i]->size);
    for (k = 0; k < count; k++) {
      vec_zero(q_i);
      for (x = 0; x < thismax; x++) {
        for (j = max(0, x - p[i]->size + 1); j <= x; j++) 
          q_i->data[x] += q_i_1->data[j] * p[i]->data[x - j];
      }
      vec_copy(q_i_1, q_i);
    }
  }

  vec_free(q_i_1);

  /* trim very small values off tail before returning */
  for (x = q_i->size - 1; x >= 0; x--) {
    if (q_i->data[x] > epsilon) {
      q_i->size = x+1;
      break;
    }
  }

  pv_normalize(q_i);
  return q_i;
}

/* compute and return a probability vector giving Pois(x | lambda) up to
   point where < epsilon */
Vector *pv_poisson(double lambda, double epsilon) {
  int j;
  Vector *pois = vec_new(max((int)(10.0 * lambda), 50));
  vec_zero(pois);
  pois->data[0] = exp(-lambda);
  for (j = 1; j < pois->size; j++) {
    pois->data[j] = pois->data[j-1] * lambda / j;
    if (pois->data[j] < epsilon) {
      pois->size = j+1;
      break;
    }
  }
  return pois;
}

/* convolve distribution n times, using a faster algorithm than the
   ones above; time is proportional to log(n) rather than n */
Vector *pv_convolve_fast(Vector *p, int n, double epsilon) {
  int i, j, checksum;
  int logn = log2_int(n);
  Vector *pow_p[64], *pows[64];
  Vector *retval;

  if (n == 1)
    return vec_create_copy(p);

  /* Let p^i be the ith convolution of p (i.e., p o p o ... o p, i
     times, where 'o' is the convolution operator).  We use the fact
     that p^i o p^j = p^(i+j) to compose p^n from p^1, p^2, p^4, p^8,
     ... in log2(k) steps */

  /* compute "powers" of p */
  pow_p[0] = p;
  for (i = 1; i <= logn; i++) 
    pow_p[i] = pv_convolve(pow_p[i-1], 2, epsilon);

  /* now combine powers to get desired convolution */
  j = checksum = 0;
  for (i = 0; i <= logn; i++) {
    unsigned bit_i = (n >> i) & 1;
    if (bit_i) {
      pows[j++] = pow_p[i];
      checksum += int_pow(2, i);
    }
  }
  if (checksum != n)
    die("ERROR pv_convolve_fast checksum (%i) != n (%i)\n", checksum, n);
  
  retval = pv_convolve_many(pows, NULL, j, epsilon);

  for (i = 1; i <= logn; i++) 
    vec_free(pow_p[i]);

  return retval;
}

/* compute CDF based on probability vector.  If side == UPPER,
   computes cumulative probabilities for right tail rather than
   left */
Vector *pv_cdf(Vector *pdf, p_val_type side) {
  Vector *cdf = vec_new(pdf->size);
  int x;
  if (side == LOWER) {
    cdf->data[0] = pdf->data[0];
    for (x = 1; x < pdf->size; x++) 
      cdf->data[x] = cdf->data[x-1] + pdf->data[x];
  }
  else {
    cdf->data[pdf->size-1] = pdf->data[pdf->size-1];
    for (x = pdf->size-2; x >= 0; x--) 
      cdf->data[x] = cdf->data[x+1] + pdf->data[x];
  }
  return cdf;
}
/** Given a probability vector, draw an index.  Call srandom externally */

int pv_draw_idx_arr(double *arr, int n) {
  int i;
  double sum, r;
  sum = 0;
  r = unif_rand();
  for (i = 0; i < n; i++) {
    sum += arr[i];
    if (r < sum) break;
  }
  if (i == n) i = n-1;    /* to be safe */
  return i;
}
  

/** Given a probability vector, draw an index.  Call srandom externally */
int pv_draw_idx(Vector *pv) {
  return pv_draw_idx_arr(pv->data, pv->size);
}


