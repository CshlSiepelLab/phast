/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/**
   @file fit_feature.h
   Functions to compute likelihoods, estimate scale factors, perform
   LRTs, score tests, etc. for multi-column features.  Generalization
   of fit_columm.c to features.  Makes use of several of the
   single-column functions. 
   @ingroup phylo
*/

#ifndef FIT_FEAT_H
#define FIT_FEAT_H

#include <fit_column.h>

/** metadata for fitting scale factors to individual alignment columns. */
typedef struct {
  GFF_Feature *feat; /**< Feature */
  ColFitData *cdata; /**< Column Fit */
} FeatFitData;

/** \name Feature Fit Data allocation function
 \{ */

/** Create Feature Fit Data object with metadata and scratch memory for fitting scale factors.
  @param[in] mod Tree Model to initialize Fit Fit Data in
  @param[in] msa Multiple Sequence Alignment sequence data
  @param[in] stype How much of the tree to use i.e. ALL, SUBTREE
  @param[in] mode which type of scoring to use i.e. CON, ACC, NNEUT, CONACC   
  @param[in] second_derivs initial second derivative for Feature Fit Data
*/
FeatFitData *ff_init_fit_data(TreeModel *mod,  MSA *msa, scale_type stype, 
                              mode_type mode, int second_derivs);

/** \name Feature Fit Data cleanup function
 \{ */

/** Free metadata and memory for fitting scale factors */
void ff_free_fit_data(FeatFitData *d);

/** \name Feature Fit Data likelihood calculation functions 
 \{ */

/** Estimate parameters
  @param params scales used with the model
  @param data Feature Fit data to estimate parameters with
  @result data must be able to cast to type FeatFitData
*/
double ff_likelihood_wrapper(Vector *params, void *data);

/** Compute and return the log likelihood of a tree model with respect
   to a given feature.  
   @param mod Tree Model to perform log likelihood on
   @param msa Multiple Sequence Alignment sequence data
   @param feat feature to compute log likelihood with respect of
   @param scratch pre-allocated memory to calculate log likelihood */ 
double ff_compute_log_likelihood(TreeModel *mod, MSA *msa, GFF_Feature *feat,
                                 double **scratch);

/** \name Feature Fit Data likelihood ratio test functions 
 \{ */

/**  Perform a likelihood ratio test for multiple features.  
  Compares the given null model with an alternative model
  that has free scaling parameter for all branches.  
  @param[in] mod Tree Model to perform likelihood ratio test on
  @param[in] msa Multiple Sequence Alignment sequence data.
  @param[in] feats Features to perform likelihood ratio tests for
  @param[in] mode Which type of scoring to use
  @param[out] feat_pvals (Optional) Computed p-values for each feature using chi-squared distribution
  @param[out] feat_scales (Optional) Computed scale factors for each feature
  @param[out] feat_llrs (Optional) raw likelihood ratios
  @param logf Location to save output
  @note Assumes a 0th order model, leaf-to-sequence mapping 
    already available, probability matrices computed, sufficient statistics available.
*/
void ff_lrts(TreeModel *mod, MSA *msa, GFF_Set *feats, mode_type mode, 
             double *feat_pvals, double *feat_scales, double *feat_llrs, 
             FILE *logf);

/**  Perform a likelihood ratio test for multiple features on a subtree.  
  Compares the given null model with an alternative model
  that has free scaling parameter for all branches.  
  @param[in] mod Tree Model to perform likelihood ratio test on
  @param[in] msa Multiple Sequence Alignment sequence data.
  @param[in] gff Features to perform likelihood ratio tests for
  @param[in] mode Which type of scoring to use
  @param[out] feat_pvals (Optional) Computed p-values for each feature using chi-squared distribution
  @param[out] feat_null_scales (Optional) Scales for null hypothesis
  @param[out] feat_scales (Optional) Computed scale factors for each feature
  @param[out] feat_sub_scales (Optional) Scales for sub optimal alternative hypothesis
  @param[out] feat_llrs (Optional) raw likelihood ratios
  @param logf Location to save output
  @note Assumes a 0th order model, leaf-to-sequence mapping 
    already available, probability matrices computed, sufficient statistics available.

*/
void ff_lrts_sub(TreeModel *mod, MSA *msa, GFF_Set *gff, mode_type mode, 
                 double *feat_pvals, double *feat_null_scales, 
                 double *feat_scales, double *feat_sub_scales, 
                 double *feat_llrs, FILE *logf);

/** \name Feature Fit Data derivative calculation functions 
 \{ */

/**  Compute the first and (optionally) second derivatives with respect
   to the scale parameter for the single-feature log likelihood.
  @param[in] d Feature Data to analyze 
  @param[out] first_deriv first derivative of column data
  @param[out] second_deriv (Optional) second derivative of feature data
  @param scratch pre allocated memory for performing derivative calculations
  @result log likelihood
  @note This version assumes a single scale parameter; see below for the subtree version   
 */
double ff_scale_derivs(FeatFitData *d, double *first_deriv,
                       double *second_deriv, double ***scratch);

/**  Compute the first and (optionally) second derivatives with respect to the scale parameters for the single-feature log likelihood.
  @param[in] d Feature Data to analyze
  @param[out] gradient gradient from first partial derivative
  @param[out] hessian (Optional) second order partial derivative
  @param[in] scratch pre-allocated memory for performing derivative calculations
  @result log likelihood
  @note This version assumes scale parameters for the whole tree and for the subtree
*/
double ff_scale_derivs_subtree(FeatFitData *d, Vector *gradient, 
                                 Matrix *hessian, double ***scratch);

/** \name Feature Fit Data gradient calculation functions 
 \{ */

/** Estimate gradient based on feature data.
    Used in parameter estimation
    @parm grad[out] gradient
    @param params NOT USED
    @param[in] data feature data used to determine gradient
    @param lb NOT USED
    @param ub NOT USED
    @warning data must be able to cast to type FeatFitData
*/
void ff_grad_wrapper(Vector *grad, Vector *params, void *data, 
                     Vector *lb, Vector *ub);

/** \name Feature Fit Data scoring functions 
 \{ */


/** Calculate scores of full tree using feature fit data
  @param[in] mod Tree Model to perform likelihood test on  
  @param[in] msa Multiple Sequence Alignment sequence data.
  @param[out] feat_pvals (Optional) Computed p-values 
  @param[out] feat_derivs (Optional) Computed first derivatives by feature
  @param[out] feat_teststats (Optional) Statistics for each test (first_derivative^2/fim)
  @see col_score_tests_sub
*/
void ff_score_tests(TreeModel *mod, MSA *msa, GFF_Set *gff, mode_type mode, 
                    double *feat_pvals, double *feat_derivs, 
                    double *feat_teststats);

/** Calculate scores of subtree using feature fit data.
  @param mod[in Tree model to perform likelihood test on
  @param[in] msa Multiple Sequence Alignment sequence data.
  @param[in] mode What type of score to produce i.e. CON, ACC, NNEUT, CONACC
  @param[out] feat_pvals (Optional) computed p-values using chi-squared distribution
  @param[out] feat_null_scales (Optional) computed null hypothesis scales
  @param[out] feat_derivs (Optional) first derivatives by tuple column
  @param[out] feat_sub_derivs (Optional) derivatives for sub optimal
  @param[out] feat_teststats (Optional) statistics or each test (first_derivative^2/fit)
*/
void ff_score_tests_sub(TreeModel *mod, MSA *msa, GFF_Set *gff, mode_type mode,
                        double *feat_pvals, double *feat_null_scales, 
                        double *feat_derivs, double *feat_sub_derivs, 
                        double *feat_teststats, FILE *logf);


/** Perform a GERP-like computation to compute conservation scores for each feature.
   @param[in] mod Tree Model to analyze
   @param[in] msa Multiple Sequence Alignment sequence data
   @param[in] mode What type of score to produce i.e. CON, ACC, NNEUT, CONACC
   @param[out] feat_nneut (Optional) exact number of substitutions under neutrality
   @param[out] feat_nobs (Optional) expected number of substitutions after re-scaling
   @param[out] feat_nrejected (Optional) expected number of rejected substitutions
   @param[out] feat_nspec (Optional) number of species with data
   @note Gaps and missing data are handled by working with the induced subtree.
 */
void ff_gerp(TreeModel *mod, MSA *msa, GFF_Set *gff, mode_type mode, 
             double *feat_nneut, double *feat_nobs, double *feat_nrejected, 
             double *feat_nspec, FILE *logf);

/** \name Feature Fit Data check sufficient data to perform analysis functions
 \{ */

/** Identify branches wrt which a given column feature are not informative,
   in the sense that all leaves beneath these branches having only missing
   data.  Will set (preallocated) array has_data[i] = I(branch above
   node i is informative).  Will also set *nspec equal to number of
   @param mod Tree Model to use
   @param msa Multiple Sequence Alignment sequence data
   @param feat feature we are interested in
   @param has_data which nodes have data (indexed by node id)
   @param nspec number of leaves that have data
 */
void ff_find_missing_branches(TreeModel *mod, MSA *msa, GFF_Feature *feat, 
                              int *has_data, int *nspec);


/** Check if a feature has enough data to perform an analysis.
    @param mod Tree Model 
    @param msa Multiple Sequence Alignment sequence data
    @param feature Feature to check for enough data
    @result TRUE if at least one column in the feature has two or more 
    actual bases (not gaps or missing data), otherwise FALSE   
 */
int ff_has_data(TreeModel *mod, MSA *msa, GFF_Feature *f);

/**  Check if a feature has enough data in the subtree of interest to perform an analysis.
    @param mod Tree Model
    @param msa Multiple Sequence Alignment sequence data
    @param f Feature to check for enough data
    @param inside list of nodes inside the subtree
    @param outside list of nodes outside the subtree
    @result TRUE if at least one column in feature has at least one base
    in the subtree of interest, at least one in the supertree of interest, 
    and at least three bases total (the minimum required for a meaningful 
    subtree test), otherwise returns FALSE.  
    @warning If inside and outside are both NULL, then returns TRUE if there
   are at least three bases total.
*/
int ff_has_data_sub(TreeModel *mod, MSA *msa, GFF_Feature *f, List *inside,
		                    List *outside);
/** \} */
#endif
