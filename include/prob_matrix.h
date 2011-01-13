/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file prob_matrix.h
   Libraries and functions to operate on matrices representing discrete 
   bivariate probability distributions defined for non-negative integers.  

   General idea is element (x,y)
   of matrix M (x,y >= 0) represents p(x,y).  With long-tailed
   distributions, matrix dimensions are truncated at size x_max, y_max
   such that p(x,y) < epsilon for x >= x_max, y >= x_max, where
   epsilon is an input parameter. 
   
   A stochastic matrix describes a Markov chain over a finite state space.
   If the probability of moving from x to y in one time step is Pr(y | x) = Px,y, the stochastic matrix P is given by using Pi,j as the ith row and jth column element

   @ingroup base
*/

#ifndef PROB_MATRIX
#define PROB_MATRIX

#include <vector.h>
#include <matrix.h>

/** Probability matrix epsilon value */
#define PM_EPS 1e-10

/** Calculate mean x & mean y for probability matrix 
    @param[in] p Probability matrix
    @param[out] mean_x Mean of X probabilities
    @param[out] mean_y Mean of Y probabilities
*/
void pm_mean(Matrix *p, double *mean_x, double *mean_y);

/** \name Marginal distribution functions 
\{ */
/** Return marginal distribution for x, sum_y p(x, y) as vector
  @param p Probability Matrix 
  @result Sum of p[x][y] grouped by x with each element in the vector being a different x
*/
Vector *pm_marg_x(Matrix *p);

/** Return marginal distribution for y, sum_x p(x, y) as vector
  @param p Probability Matrix 
  @result Sum of p[x][y] grouped by y with each element in the vector being a different y
*/
Vector *pm_marg_y(Matrix *p);

/** Return marginal distribution for x + y, p(n) = sum_{x,y:x+y=n} p(x,y) as vector
  @param p Probability Matrix
  @result Sum of p[x][y] grouped by (x+y)
*/
Vector *pm_marg_tot(Matrix *p);

/** \} */

/** Normalize matrix so that all elements sum to one.
  @param[in,out] p Probability Matrix 
*/
void pm_normalize(Matrix *p);


/** Compute means, variances, and covariance of probability matrix.
   @param[in] p Matrix containing data to analyze
   @param[out] mean_x Mean of X probabilities
   @param[out] mean_y Mean of Y probabilities
   @param[out] var_x Variance of X probabilities
   @param[out] var_y Variance of Y probabilities
   @param[out] covar Covariance of probability matrix
*/  
void pm_stats(Matrix *p, double *mean_x, double *mean_y, double *var_x, 
              double *var_y, double *covar);

/** \name Conditional Probability Distribution functions
\{ */


/** Return conditional probability distribution of x given x+y  
  @param p Probability Matrix
  @param tot Given total, x+y
*/
Vector *pm_x_given_tot(Matrix *p, int tot);

/** Return conditional probability distribution of y given x+y
   @param p Probability Matrix
   @param tot Given total, x+y
*/
Vector *pm_y_given_tot(Matrix *p, int tot);

/** Return conditional probability distribution of x given x+y assuming bivariate normal
  @param p Probability Matrix
  @param tot Given total, x+y
  @param mu_x Mean X
  @param mu_y Mean Y
  @param sigma_x Standard deviation X
  @param sigma_y Standard deviation Y
  @param rho Correlation coefficient
  @note Only the region of the bivariate normal with x,y >= 0 is considered.  
  @note For use in central limit theorem approximations. */
Vector *pm_x_given_tot_bvn(int tot, double mu_x, double mu_y, 
                           double sigma_x, double sigma_y, double rho);

/** Return conditional probability distribution of y given x+y assuming bivariate normal
  @param p Probability Matrix
  @param tot Given total, x+y
  @param mu_x Mean X
  @param mu_y Mean Y
  @param sigma_x Standard deviation X
  @param sigma_y Standard deviation Y
  @param rho Correlation coefficient
  @note Only the region of the bivariate normal with x,y >= 0 is considered.  
  @note For use in central limit theorem approximations. */
Vector *pm_y_given_tot_bvn(int tot, double mu_x, double mu_y, 
                           double sigma_x, double sigma_y, double rho);

/** Return conditional probability distribution of x given x+y assuming independence of x and y
  @param p Probability Matrix
  @param tot Given total, x+y
  @param mu_x Mean X
  @param mu_y Mean Y
  @param sigma_x Standard deviation X
  @param sigma_y Standard deviation Y
  @param rho Correlation coefficient
  @note Computes desired conditional distribution from their marginals
*/
Vector *pm_x_given_tot_indep(int tot, Vector *marg_x, Vector *marg_y);

/** Return conditional probability distribution of y given x+y assuming independence of x and y
  @param p Probability Matrix
  @param tot Given total, x+y
  @param mu_x Mean X
  @param mu_y Mean Y
  @param sigma_x Standard deviation X
  @param sigma_y Standard deviation Y
  @param rho Correlation coefficient
  @note Computes desired conditional distribution from their marginals
*/
Vector *pm_y_given_tot_indep(int tot, Vector *marg_x, Vector *marg_y);
/** \} */


/** \name Convolution functions
\{ */

/** Convolve distribution 'n' times. (Slower)
  @param p Probability Matrix
  @param n Amount of times to convolve distribution
  @param epsilon Trim values less than epsilon off tail
  @result Convolved matrix
*/
Matrix *pm_convolve(Matrix *p, int n, double epsilon);

/** Convolve distribution n times and keep all intermediate
   distributions.
  @param p Probability Matrix
  @param n Amount of times to convolve distribution
  @param epsilon Trim values less than epsilon off tail 
  @result Array q such that q[i] (1 <= i <= n) is the ith convolution of p (q[0] will be NULL)
 */
Matrix **pm_convolve_save(Matrix *p, int n, double epsilon);

/** Take convolution of a set of probability matrices (slower)
   @param p Array of Probability Matrices
   @param counts (Optional) Array of multiplicities, one of each distribution in p; Defaults to 1 per dist.
   @param n Amount of times to convolve distribution
   @param epsilon Trim values less than epsilon off tail
   @result Convolved matrix
 */
Matrix *pm_convolve_many(Matrix **p, int *counts, int n, double epsilon);

/** Take convolution of a set of probability matrices (Faster)
    @param p Array of probability Matrices
    @param n Number of times to convolve distribution
    @param max_nrows Maximum number of rows of result matrix
    @param max_ncols maximum number of columns of result matrix
    @result Convolved matrix
    @note Dos not take counts, normalize, or trim dimension
*/
Matrix *pm_convolve_many_fast(Matrix **p, int n, int max_nrows, int max_ncols);

/** Convolve distribution 'n' times. (Faster)
  @param p Probability Matrix
  @param n Amount of times to convolve distribution
  @param epsilon Trim values less than epsilon off tail
  @result Convolved matrix
*/
Matrix *pm_convolve_fast(Matrix *p, int n, double epsilon);

/** \} */

#endif
