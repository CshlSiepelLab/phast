/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/


/**
  @file indel_mod.h
  Simple model of insertions and deletions, assumes given indel history 
  @ingroup phylo
*/

#ifndef INDMOD_H
#define INDMOD_H

#include <indel_history.h>

/** Number of indel states i.e. insert, delete, base */
#define NINDEL_STATES 3
//need to make sure ERROR isn't already defined from R libraries
#undef ERROR
/** Possible indel results when comparing two sites. */
typedef enum {MATCH,  /**< Parent and child both 'base' at same site. */ 
	      CHILDINS, /**< Parent has 'insertion' while child has 'base' at same site. */
	      CHILDDEL, /**< Parent has 'base' while child has 'insertion' at same site. */
	      SKIP,   /**< Parent and child both have ('insertion' at same site) or ('deletion' at same site). */
	      ERROR  /**< There was an error in comparing sites between parent and child. */
	     } col_type;

/** Indel model for a single branch */
typedef struct {
  double alpha; /**< Rate of insertion per expected substitution per site for this branch */
  double beta; /**< Rate of deletion per expected substitution per site for this branch */
  double tau; /**< Roughly equal to the inverse of the expected indel length (modulo
    adjustments required to make probabilities sum to one) for this branch */
  double t; /** ??? */
  MarkovMatrix *probs; /**< Indel substitution pattern for this branch*/
  Matrix *log_probs;  /**< Log of substitution pattern probabilities for this branch*/
  Vector *beg_probs;  /**< Beginning probabilities for this branch */
  Vector *beg_log_probs; /**< Log of beg_probs for this branch */
} BranchIndelModel;

/** Indel model for a dataset */
typedef struct {
  double alpha; /**< Overall rate of insertion per expected substitution per site */
  double beta; /**< Overall rate of deletion per expected substitution per site */
  double tau;  /**< Roughly equal to the inverse of the expected indel length (modulo 
    adjustments required to make probabilities sum to one) for this branch */
  double training_lnl; /**< Log likelihood from performing training (if any) */
  TreeNode *tree; /**< Tree representing structure of dataset */
  BranchIndelModel **branch_mods; /**< Indel models for each branch */
} IndelModel;

/** Model sufficient statistics for per branch */
typedef struct {
  Matrix *trans_counts; /**< Transition counts */
  Vector *beg_counts;   /**< Beginning counts */
} BranchIndelSuffStats;

/** Model sufficient statistics per dataset */
typedef struct {
  TreeNode *tree; /**< Tree structure */
  BranchIndelSuffStats **branch_counts; /**< List of model statistics per branch */
} IndelSuffStats;

/** \name Indel (Branch) Model allocation functions 
 \{ */

/** Create new Indel model for a branch 
    @param alpha Rate of insertion per expected substitution per site
    @param beta Rate of deletion per expected substitution per site 
    @param tau Roughly equal to the inverse of the expected indel length 
               (modulo adjustments required to make probabilities sum to one) for this branch 
    @param t t
    @result New branch indel model
*/
BranchIndelModel *im_new_branch(double alpha, double beta, double tau,
                                double t);

/** Create new Indel model containing all branches in a tree, each initialized to the exact same alpha, beta, tau.
    @param alpha Rate of insertion per expected substitution per site
    @param beta Rate of deletion per expected substitution per site
    @param tau Roughly equal to the inverse of the expected indel length 
               (modulo adjustments required to make probabilities sum to one) for this branch 
    @param tree Tree structure
    @result New branch indel model
*/
IndelModel *im_new_all(double alpha, double beta, double tau,  TreeNode *tree);

/** Create new Indel model containing all branches in a tree, each with its own alpha, beta, tau.
    @param alpha Rates of insertion per expected substitution per site
    @param beta Rates of deletion per expected substitution per site 
    @param tau Roughly equal to the inverse of the expected indel length 
               (modulo adjustments required to make probabilities sum to one) for this branch 
    @param tree Tree structure
    @result New branch indel model
*/
IndelModel *im_new(double *alpha, double *beta, double *tau, 
                   TreeNode *tree);

/** \name Indel (Branch) Model cleanup functions 
 \{ */

/**  Free a branch indel model. 
     @param bim Branch Model to free
*/
void im_free_branch(BranchIndelModel *bim);

/** Free a indel model.
    @param im Indel Model to free
*/
void im_free(IndelModel *im);

/** Free a statistics object. 
    @param iss Statistics object to free
*/
void im_free_suff_stats(IndelSuffStats *iss);


/** \name Indel (Branch) Model modification functions 
 \{ */

/**  Set values of alpha, beta, tau, t for a branch indel model
     @param bim Branch Indel Model to modify
     @param alpha rate of insertion per expected substitution per site
     @param beta rate of deletion per expected substitution per site 
     @param tau Roughly equal to the inverse of the expected indel length 
               (modulo adjustments required to make probabilities sum to one) for this branch 
*/
void im_set_branch(BranchIndelModel *bim, double alpha, double beta, 
                   double tau, double t);

/**  Set values of alpha, beta, tau, t for EVERY branch in indel model
    @param im Indel Model to modify
    @param alpha Rate of insertion per expected substitution per site
    @param beta Rate of deletion per expected substitution per site 
    @param tau Roughly equal to the inverse of the expected indel length 
               (modulo adjustments required to make probabilities sum to one) for this branch 
    @param tree Tree structure
*/
void im_set_all(IndelModel *im, double alpha, double beta, double tau, 
                TreeNode *tree);

/**  Set per branch values of alpha, beta, tau, tree for EVERY branch in indel model
    @param im Indel Model to modify
    @param alpha Rates of insertion per expected substitution per site
    @param beta Rates of deletion per expected substitution per site 
    @param tau Roughly equal to the inverse of the expected indel length 
               (modulo adjustments required to make probabilities sum to one) for this branch 
    @param tree Tree structure
*/
void im_set(IndelModel *im, double *alpha, double *beta, double *tau, 
            TreeNode *tree);

/** \name Indel (Branch) Model Log Likelihood functions 
 \{ */

/** Calculate log likelihood for a branch.
    @param[in] ih Indel History
    @param[in] bim Branch to compute log likelihood for
    @param[in] child ID of child within branch
    @param[out] col_logl Computed log likelihoods per column in branch
    @result Log likelihood for entire branch
 */
double im_branch_column_logl(IndelHistory *ih, BranchIndelModel *bim, 
                             int child_id, double *col_logl);

/** Calculate total log likelihood.
    @param[in] ih Indel History
    @param[in] im Indel Model to get log likelihood for
    @param[out] col_logl Computed log likelihoods for each branch
    @result Column log likelihood
*/
double im_column_logl(IndelHistory *ih, IndelModel *im, double *col_logl);

/** Calculate total log likelihood.
    @param[in] im Indel Model to get log likelihood for
    @param[in] iss Indel Sufficient Statistics
    @result Tree log likelihood
*/
double im_likelihood(IndelModel *im, IndelSuffStats *iss);

/** \name Indel (Branch) Model sufficient statistics functions
 \{ */

/** Create Branch Indel Sufficient Statistics.
    @param ih Indel History to calculate statistics on
    @param child_id ID of first node of branch to get statistics for in tree
    @result Branch indel sufficient statistics
*/
BranchIndelSuffStats *im_suff_stats_branch(IndelHistory *ih, int child_id);

/** Create Indel Sufficient Statistics.
    @param ih Indel History to create statistics for
    @result Indel sufficient statistics
*/
IndelSuffStats *im_suff_stats(IndelHistory *ih);


/* Not implemented 
double im_simulate_history(IndelModel *tim, int ncols);
*/


/** Collect sufficient statistics for a branch, considering only sites in the specified category.
    @param ih Indel History to calculate Sufficient Statistics for
    @param child_id First node of branch to obtain sufficient statistics for
    @param categories Integers specifying which category each column is in
    @param do_cat Integer specifying which category is accepted 
    @result Collected sufficient statistics for the branch specified
*/
BranchIndelSuffStats *im_suff_stats_branch_cat(IndelHistory *ih, int child_id,
                                               int *categories, int do_cat);

/** Collect sufficient statistics for an entire tree, considering only sites in the specified category.
    @param ih Indel History to calculate Sufficient Statistics for
    @param categories Integers specifying which category each column is in
    @param do_cat Integer specifying which category is accepted 
    @result Collected sufficient statistics for the tree specified
*/
IndelSuffStats *im_suff_stats_cat(IndelHistory *ih, int *categories, 
                                  int do_cat);
/** \} */

/** Estimate alpha, beta, and tau from an indel history by maximum likelihood.
    @param im Indel Model 
    @param ih Indel History
    @param ss Indel Sufficient Statistics
    @param logf Log file to write to
*/
void im_estimate(IndelModel *im, IndelHistory *ih, IndelSuffStats *ss, 
                 FILE *logf);

#endif
