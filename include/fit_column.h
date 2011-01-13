/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file fit_column.h
   Functions to compute likelihoods for individual alignment columns,
   estimate column-by-column scale factors by maximum likelihood,
   perform single-base LRTs, score tests, etc
   @ingroup phylo
*/ 

#ifndef FIT_COL_H
#define FIT_COL_H

#include <misc.h>
#include <tree_model.h>
#include <msa.h>
#include <vector.h>
#include <matrix.h>
#include <complex_matrix.h>

/** Portions of tree that can be used. */
typedef enum {ALL, /**< Use entire tree. */
	      SUBTREE /**< Use only part of the tree. */
	     } scale_type;

/** Types of scores that can be produced. */
typedef enum {CON, /**< Conservation scores */
	      ACC, /**< Acceleration scores */
	      NNEUT, /**< Non-Neutrality scores */
	      CONACC /**< Summarize conservation and acceleration scores */
	     } mode_type;

/** Metadata for fitting scale factors to individual alignment columns. */
typedef struct {
  TreeModel *mod;             /**< Pointer to Tree Model this column belongs to */
  MSA *msa;                   /**< Pointer to MSA containing this column */
  int tupleidx;				/**< Specifies which column this is in the MSA */
  scale_type stype;             /**< Whether doing all-branches or
                                   supertree/subtree estimation.  */
  mode_type mode;               /**< Type of parameter bounding.  */
  int second_derivs;            /**< Whether or not second derivatives
                                   need to be computed. */
  Vector *params;		/** Parameters */
  Vector *lb;			/**< Lower bound of parameters. */
  Vector *ub;			/**< Upper bound  of parameters. */
  double init_scale;		/**< Initial scale for parameters. */
  double init_scale_sub;	/**< Initial amount to subtract for scale. */
  Matrix ***PP;			
  Matrix ***PPP;
  Matrix ***QQ;
  Matrix ***QQQ;
  Matrix ***RRR;
  Zvector *expdiag_z;		/**< Complex number exponential diagonalized matrix. */
  Vector *expdiag_r;		/**< Real number exponential diagonalized matrix. */

  double ***fels_scratch;       /**< Scratch memory for Felsenstein's
                                   Algorithm (likelihoods and derivatives). */
  int nfels_scratch;            /**< Number of scratch arrays (depends on mode). */
  Zmatrix *mat_scratch_z;         /**< Scratch memory for complex number matrix manipulation. */
  Zvector *vec_scratch1_z, *vec_scratch2_z; /**< Scratch memory for complex number vector manipulation. */
  Vector *vec_scratch1_r, *vec_scratch2_r; /**< Scratch memory for real number vector manipulation. */
  double deriv2;                /**< second derivative for 1d case. */
} ColFitData;

/* data for grid of pre-computed Fisher Information Matrices */
#define GRIDSIZE1 0.02
#define GRIDSIZE2 0.05
#define GRIDMAXLOG 3

/** Grid of Fisher Information Matrices. */
typedef struct {
  double *scales;               /**< Scale factors in grid */
  int ngrid1;                   /**< Number of grid points between 0 and
                                   1 (linear) */
  int ngrid2;                   /**< Number of grid points above 1 (log
                                   linear) */
  int ngrid;                    /**< Total number of grid points */
  Matrix **fim;                 /**< Pre-computed FIMs */
} FimGrid;

/** \name Column Fit Data allocation function
 \{ */

/** Create Column Fit Data object with metadata and scratch memory for fitting scale factors.
  @param[in] mod Tree Model to initialize Column Fit Data in
  @param[in] msa Multiple Sequence Alignment sequence data
  @param[in] stype Use whole tree or subtree
  @param[in] mode which type of scoring to use i.e. CON, ACC, NNEUT, CONACC   
  @param[in] second_derivs initial second derivative for Column Fit Data
*/
ColFitData *col_init_fit_data(TreeModel *mod, MSA *msa, scale_type stype,
                              mode_type mode, int second_derivs);

/** \name Column Fit Data cleanup function
 \{ */

/** Free metadata and memory for fitting scale factors
   @param d Column Fit Data to free */
void col_free_fit_data(ColFitData *d);

/** \name Column Fit Data likelihood calculation functions
 \{ */

/** Estimate parameters
  @param Parameter scales to use with model
  @param data Column Fit Data to estimate parameters with
  @result Estimated log likelihood
  @warning data must be able to cast to type ColFitData
*/
double col_likelihood_wrapper(Vector *params, void *data);

/** Wrapper for likelihood function for use in parameter estimation;
   version for use with opt_newton_1d 
  @param x scale to use with model 
  @param data Column Fit Data
  @result Estimated log likelihood
  @warning data must be able to cast to type ColFitData
*/
double col_likelihood_wrapper_1d(double x, void *data);


/** Compute and return the log likelihood of a tree model with respect
   to a single column tuple in an alignment.  

   This is a pared-down
   version of tl_compute_log_likelihood for use in estimation of
   base-by-base scale factors.  It assumes a 0th order model,
   leaf-to-sequence mapping already available, probability matrices computed,
   sufficient stats already available.  This function does allow for
   rate variation.
  @param mod Substitution model, rates and its metadata
  @param msa Sequence data and its metadata
  @param tupleidx Which column to compute log likelihood for
  @param scratch Pre-allocated memory used as scratch space when computing likelihood 
  @result Estimated log likelihood
  @note This function is simply a wrapper for col_compute_likelihood and computes log() on the result
  @note Uses log rather than log2
*/

double col_compute_log_likelihood(TreeModel *mod, MSA *msa, int tupleidx,
                                  double **scratch);

/** Compute and return the likelihood of a tree model with respect
   to a single column tuple in an alignment.  

   This is a paired-down
   version of tl_compute_likelihood for use in estimation of
   base-by-base scale factors.  It assumes a 0th order model,
   leaf-to-sequence mapping already available, probability matrices computed,
   sufficient statistics already available.  This function does allow for
   rate variation.
  @param mod Substitution model, rates and its metadata
  @param msa Sequence data and its metadata
  @param tupleidx Which column to compute likelihood for
  @param scratch Pre-allocated memory used as scratch space when computing likelihood 
  @result Estimated Log likelihood
*/
double col_compute_likelihood(TreeModel *mod, MSA *msa, int tupleidx,
		              double **scratch);

/** \name Column Fit Data likelihood ratio test functions
 \{ */

/** Perform a likelihood ratio test for each column tuple in an
   alignment.

   Compare the given null model with an alternative model
   that has a free scaling parameter for all branches.  Assumes a 0th
   order model, leaf-to-sequence mapping already available, probability
   matrices computed, sufficient statistics available.
   @param[in] mod Tree Model to perform likelihood ratio test on
   @param[in] msa Multiple Sequence Alignment sequence data
   @param[out] tuple_pvals (Optional) Computed p-values using chi-squared distribution
   @param[out] tuple_scales (Optional) Computed individual scale factors
   @param[out] tuple_llrs (Optional) raw likelihood ratios
   @param logf Location to save output
   @note Must define mode as CON (for 0 <= scale <= 1), ACC
   (for 1 <= scale), NNEUT (0 <= scale), or CONACC (0 <= scale) */ 

void col_lrts(TreeModel *mod, MSA *msa, mode_type mode, double *tuple_pvals, 
              double *tuple_scales, double *tuple_llrs, FILE *logf);

/** Perform a likelihood ratio test for each column tuple in a subtree of an alignment.
    @param[in] mod Tree Model to perform likelihood test on
    @param[in] msa Multiple Sequence Alignment sequence data
    @param[out] tuple_pvals (Optional) computed p-values using chi-squared distribution
    @param[out] tuple_null_scales (Optional) Scales for Null hypothesis
    @param[out] tuple_scales (Optional) Scales for optimal alternative hypothesis
    @param[out] tuple_sub_scales (Optional) Scales for sub optimal alternative hypothesis
    @param[out] tuple_llrs (Optional) Log Likelihood RS ratio
    @param[in] logf output file to write to
    @see col_grad_wrapper
*/ 
void col_lrts_sub(TreeModel *mod, MSA *msa, mode_type mode, 
                  double *tuple_pvals, double *tuple_null_scales, 
                  double *tuple_scales, double *tuple_sub_scales, 
                  double *tuple_llrs, FILE *logf);


/** \name Column Fit Data derivative calculation functions
 \{ */


void col_scale_derivs_num(ColFitData *d, double *first_deriv, 
                          double *second_deriv);

void col_scale_derivs_subtree_num(ColFitData *d, Vector *gradient, 
                                  Matrix *hessian);

/**
   Compute the first and (optionally) second derivatives with respect
   to the scale parameter for the single-column log likelihood.

   @param[in] d Column Data to analyze
   @param[out] first_deriv first derivative of column data
   @param[out] second_deriv (optional) second derivative of column data if == NULL will not be computed
   @param[in] scratch pre allocated memory for performing derivative calculations
   @result log likelihood
   @note This version assumes a single scale parameter; see below for the subtree version.
*/
double col_scale_derivs(ColFitData *d, double *first_deriv, 
                        double *second_deriv, double ***scratch);

/** Compute the first and (optionally) second derivatives with respect
   to the scale parameters for the single-column log likelihood
   function (col_compute_log_likelihood).
   @param[in] d Column Data to analyze
   @param[out] gradient gradient from first partial derivative
   @param[out] hessian (optional) second order partial derivative
   @param[in] scratch pre-allocated memory for performing derivative calculations
   @result Estimated log likelihood
   @note  This version assumes scale parameters for the whole tree and for the subtree.
  */
double col_scale_derivs_subtree(ColFitData *d, Vector *gradient, 
                                Matrix *hessian, double ***scratch);

/** \name Column Fit Data gradient calculation functions
 \{ */

/** Estimate gradient based on column data.
  Used in parameter estimation.
  @param[out] grad gradient
  @param params NOT USED
  @param[in] data column data used to determine gradient 
  @param lb NOT USED
  @param UB NOT USED
  @result Estimated log likelihood
  @warning data must be able to cast to type ColFitData 
 */
void col_grad_wrapper(Vector *grad, Vector *params, void *data, 
                      Vector *lb, Vector *ub);


/** \name Column Fit Data scoring functions
 \{ */

/** Calculate scores of full tree using column fit data
  @param[in] mod Tree Model to perform likelihood test on  
  @param[in] msa Multiple Sequence Alignment sequence data.
  @param[out] tuple_pvals (Optional) Computed p-values 
  @param[out] tuple_derivs (Optional) Computed first derivatives by tuple column
  @param[out] tuple_teststats (Optional) Statistics for each test  (first_derivative^2 / fim)
  @see col_score_tests_sub
*/
void col_score_tests(TreeModel *mod, MSA *msa, mode_type mode, 
                     double *tuple_pvals, double *tuple_derivs, 
                     double *tuple_teststats);


/** Calculate scores of subtree using column fit data.
  @param mod[in Tree model to perform likelihood test on
  @param[in] msa Multiple Sequence Alignment sequence data.
  @param[in] mode What type of score to produce i.e. CON, ACC, NNEUT, CONACC
  @param[out] tuple_pvals (Optional) computed p-values using chi-squared distribution
  @param[out] tuple_null_scales (Optional) computed null hypothesis scales
  @param[out] tuple_derivs (Optional) first derivatives by tuple column
  @param[out] tuple_sub_derivs (Optional) derivatives for sub optimal
  @param[out] tuple_teststats (Optional) statistics or each test (first_derivative^2 / fim)
*/
void col_score_tests_sub(TreeModel *mod, MSA *msa, mode_type mode,
                         double *tuple_pvals, double *tuple_null_scales, 
                         double *tuple_derivs, double *tuple_sub_derivs, 
                         double *tuple_teststats, FILE *logf);



/** Perform a GERP-like computation to compute conservation scores
   for each tuple.
   @param[in] mod Tree Model to analyze
   @param[in] msa Multiple Sequence Alignment sequence data
   @param[in] mode What type of score to produce i.e. CON, ACC, NNEUT, CONACC
   @param[out] tuple_nneut (Optional) exact number of substitutions under neutrality
   @param[out] tuple_nobs (Optional) expected number of substitutions after re-scaling
   @param[out] tuple_nrejected (Optional) expected number of rejected substitutions
   @param[out] tuple_nspecies (Optional) number of species with data
   @note Gaps and missing data are handled by working with the induced subtree.
 */
void col_gerp(TreeModel *mod, MSA *msa, mode_type mode, double *tuple_nneut, 
              double *tuple_nobs, double *tuple_nrejected, 
              double *tuple_nspecies, FILE *logf);


/** \} \name Column Fit Data fisher information matrix functions
 \{ */

/** Estimate 2x2 Fisher Information Matrix (expected value of the
   negative Hessian) for the subtree case, based on a particular value
   of the scale parameter (set in calling code).  Estimation is done
   by sampling. */
Matrix *col_estimate_fim_sub(TreeModel *mod);


/** Pre-compute estimates of Fisher Information Matrix for a grid of
   possible scale parameters (subtree case) */
FimGrid *col_fim_grid_sub(TreeModel *mod);

/** Free FimGrid object  */
void col_free_fim_grid(FimGrid *g);

/** Estimate Fisher Information Matrix for the non-subtree case.
   This version does not depend on any free parameters, so no grid is
   required.  Estimation is done by sampling */
double col_estimate_fim(TreeModel *mod);

/** Retrieve estimated Fisher Information Matrix for given scale 
   function; uses linear interpolation from pre-computed grid
   @param g Fisher Information Matrix(FIM) Grid holding FIM to return
   @param scale Only return matrix with this scale
   @result Fisher Information Matrix with a scale equal to that in parameters
*/
Matrix *col_get_fim_sub(FimGrid *g, double scale);

/** \name Column Fit Data check sufficient data to perform analysis functions
 \{ */

/** Check if a column has enough data to perform an analysis.
    @param mod Tree Model 
    @param msa Multiple Sequence Alignment sequence data
    @param tupleidx Column to check for enough data
    @result TRUE if column has two or more actual bases (not gaps or 
    missing data), otherwise FALSE   

 */
int col_has_data(TreeModel *mod, MSA *msa, int tupleidx);

/**  Check if a column has enough data in the subtree of interest to perform an analysis.
    @param mod Tree Model
    @param msa Multiple Sequence Alignment sequence data
    @param tupleidx Column to check for enough data
    @param inside list of nodes inside the subtree
    @param outside list of nodes outside the subtree
    @result TRUE if column has at least one base in the subtree of
   interest, at least one in the supertree of interest, and at least
   three bases total (the minimum required for a meaningful subtree
   test), otherwise returns FALSE.  
    @warning If inside and outside are both NULL, then returns TRUE if there
   are at least three bases total.
*/
int col_has_data_sub(TreeModel *mod, MSA *msa, int tupleidx, List *inside,
		                     List *outside);

/** Identify branches wrt which a given column tuple are not informative,
   in the sense that all leaves beneath these branches having missing
   data.
   @param mod Tree Model to use
   @param msa Multiple Sequence Alignment sequence data
   @param tupleidx column tuple we are interested in
   @param[out] has_data which nodes have data (indexed by node id) array has_data[i] = I(branch above
   node i is informative).
   @param nspec number of leaves that have data
 */
void col_find_missing_branches(TreeModel *mod, MSA *msa, int tupleidx, 
                               int *has_data, int *nspec);
/** \} */
#endif
