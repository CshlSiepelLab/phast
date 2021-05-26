/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/**
   @file dgamma.h
   Gamma distribution related computations.
   @ingroup phylo
*/

#ifndef DGAMMA_H
#define DGAMMA_H

#include <stdlib.h>

/** Inverse CDF Gamma (quantile function).
    @note See Yang 1994, p. 308 (discrete Gamma paper)
*/
#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))

/** Inverse CDF chi sq (quantile function).
   @param prob Probability
   @param v Degree of freedom???
   @result -1 if error otherwise z so that Prob{x<z}=prob where x is Chi2 distributed with df=v; 0.000002<prob<0.999998
   */
double PointChi2 (double prob, double v);

/** Calculate percentage of points on a normal distribution with a given probability.
   @param prob Probability to get point estimation for
   @result -999 if error otherwise z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   @note Used by PointChi2 */
double PointNormal (double prob);

/** Calculate ln(gamma(x)) for x>0.
   @note Accurate to 10 decimal places.
   @note Used by PointChi2
   @result ln(gamma(x)) */
double LnGamma (double x);

/** Incomplete gamma ratio, otherwise known as complementary normalized incomplete gamma function.
   @param x Upper limit of integration
   @param alpha Shape parameter
   @param ln_gamma_alpha ln(gamma(alpha))
   @result -1 if error, otherwise Incomplete gamma ratio I(x,alpha)
 */
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);

/** Discretization of gamma distribution with equal proportions in each category.
   @param freqK Array of frequencies, one for each category
   @param rK Rate for each category
   @param alpha Shapes distribution
   @param beta Shapes distribution
   @param K Number of equally probable rate categories
   @param median Median mutation rate for substitution model used
   @note See Yang, 1994, p. 308
 */
int DiscreteGamma (double freqK[], double rK[], 
                   double alpha, double beta, int K, int median);

#endif
