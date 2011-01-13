/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file prob_vector.h
    Vectors representing discrete probability distributions over
    non-negative integers.  General idea is element x of vector v (x >=
    0) represents p(x).  With long-tailed distributions (e.g.,
    Poisson), vectors are truncated at size x_max such that p(y) <
    epsilon for y >= x_max, where epsilon is an input parameter. 
    @ingroup base
*/

#ifndef PROB_VECTOR
#define PROB_VECTOR

#include <vector.h>

/** Type of p-value calculated */
typedef enum {LOWER, /**< Lower tail p-value */
 UPPER, /**< Upper tail p-value */ 
 TWOTAIL /**< Two tail p-value */
} p_val_type;

/** Compute Mean and Variance of probability vector.
  @param p Probability vector
  @param mean Mean value of p
  @param var Variance of p
 */
void pv_stats(Vector *p, double *mean, double *var);

/** Compute min and max of specified confidence interval.  
  @param[in] p Probability vector
  @param[in] size Size of confidence interval (between 0 and 1)
  @param[out] interval_min Computed lower cap on confidence interval
  @param[out] interval_max Computed high cap on of confidence interval
 */
void pv_confidence_interval(Vector *p, double size, int *interval_min, 
                            int *interval_max);
/** Compute quantiles 0.00, 0.01, 0.02, ..., 1.00. based on probability vector 
    @param p Probability vector
    @result Array of length 101, such that element x contains the x/100th quantile. 
*/
int* pv_quantiles(Vector *p);


/** Return one-sided p-value: p(x <= x_0) if side == LOWER, p(x >= x_0)
   if side == UPPER.  If side == TWOTAIL, heuristically returns 2 *
   min(p(x <= x_0, p(x >= x_0)) 
   @param[in] distrib Probability vector
   @param[in] x_0 Index of distrib (between 0 and distrib->size)
   @param[in] side Type of p-value i.e. lower, upper, two tail
   @result p-value
   @note For drawbacks of this approach and discussion see Dunne et al.,
   The Statistician, 1996) 
*/
double pv_p_value(Vector *distrib, double x_0, p_val_type side);

/** Compute one-sided p-values for array of values.  
   @param[in] distrib Probability vector
   @param[in] x_0 Array of indices of distrib (between 0 and distrib->size)
   @param[in] n Number of elements in x_0
   @param[out] pvals Array of p-values corresponding to indices of x_0
   @Like pv_p_value, but saves time by computing CDF and using for all pvals
*/
void pv_p_values(Vector *distrib, double *x_0, int n, double *pvals,
                 p_val_type side);

/** Normalize distribution 
    @param p Distribution to normalize
*/
void pv_normalize(Vector *p);

/** \name Convolve vector functions 
\{ */
/** Convolve distribution 'n' times (slower)
  @param p Distribution to convolve
  @param n Number of times to convolve
  @param epsilon Trim values less than epsilon off tail
  @result Convolved vector
*/
Vector *pv_convolve(Vector *p, int n, double epsilon);

/** Convolve distribution 'n' times and keep all intermediate distributions
  @param p Distribution to convolve
  @param n Number of times to convolve
  @param epsilon Trim values less than epsilon off tail
  @result Array (q) of convolved vectors s.t. q[i] ( 1 <= i <= n) is the ith convolution of p (q[0] is null)
*/
Vector **pv_convolve_save(Vector *p, int n, double epsilon);

/** Take convolution of a set of probability vectors.  
  @param p Array of probability vectors
  @param counts (Optional) Array of multiplicities, one for each distribution in p; Defaults to 1 per dist.
  @param epsilon Trim values less than epsilon off tail 
  @result Convolved vector
 */
Vector *pv_convolve_many(Vector **p, int *counts, int n, double epsilon);

/** Convolve distribution 'n' times (faster)
  @param p Distribution to convolve
  @param n Number of times to convolve
  @param epsilon Trim values less than epsilon off tail
  @result Convolved vector 
*/
Vector *pv_convolve_fast(Vector *p, int n, double epsilon);

/** \} */

/** Compute and return a probability vector giving Pois(x | lambda) up to
   point where < epsilon 
  @param lambda Description of distribution
  @param epsilon Lowest allowed value in probability vector
  @result Probability vector
*/
Vector *pv_poisson(double lambda, double epsilon);

/** Compute CDF based on probability vector.  
  @param pdf Probability vector
  @param side Type of p-value; If side == UPPER, computes 
  cumulative probabilities for right tail rather than left
  @result CDF
*/
Vector *pv_cdf(Vector *pdf, p_val_type side);

/** Given a probability array, draw an index
   @pre Call srandom externally
   @param arr Probability to draw from
   @param n Number of elements in arr
   @result Draw from probability vector
 */
int pv_draw_idx_arr(double *arr, int n);

/** Given a probability vector, draw an index. 
   @pre Call srandom externally
   @param pdf Probability to draw from
   @result Draw from probability vector
*/
int pv_draw_idx(Vector *pdf);

#endif
