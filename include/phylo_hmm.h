/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#ifndef PHMM_H
#define PHMM_H

#include <category_map.h>
#include <hmm.h>
#include <tree_model.h>
#include <lists.h>
#include <gap_patterns.h>

/** @file phylo_hmm.h
    Phylogenetic Hidden Markov Model library used for: Model fitting, Probability estimation, Indel estimation and more.
  Allows for automatic expansion of the state
   space to accommodate features on the reverse strand, and for the
     indel model described in Siepel & Haussler, RECOMB '04.  Also
   allows for cross-product constructions involving functional states
   and rate categories (Siepel & Haussler, RECOMB '03).  
   @ingroup phylo_hmm
*/

/** Defines how to model indels in a given HMM. */
typedef enum {MISSING_DATA, /**< Indels treated as missing data */
 PARAMETERIC, /**< Indels treated as parametric data */
 NONPARAMETERIC /**< Indels treated as non-parametric data */
} indel_mode_type;

/** Package of data used during EM parameter estimation */
typedef struct {
  MSA *msa;                     /**< Data to which model is fit */
  int fix_functional,		/**< Can be used to force functional transitional probabilities to be fixed */
  fix_indel;			/**< Can be used to force indel transition probabilities to be fixed */
  double rho;			/**< Scaling constant for two-state
                                   rate-variation phylo-HMM */
  double gamma;			/**< Target coverage for two-state
                                   rate-variation phylo-HMM */
  Matrix *H;                    /**< Inverse Hessian for BFGS  */
} EmData;

/** Phylo HMM object */
typedef struct {
  CategoryMap *cm;              /**< Category map */
  HMM *hmm;                     /**< Final HMM, after reflection and
                                   allowance for rate categories  */
  HMM *functional_hmm,		/**< One of two original HMMs used to create cross
                                   product (NULL if no rate
                                   categories) */
  *autocorr_hmm;                /**< One of two original HMMs used to create cross
                                   product (NULL if no rate
                                   categories) */
  TreeModel **mods;             /**< Array of tree models, after
                                   allowance for rate categories  */
  int nmods;                    /**< Number of tree models (length of
                                   array mods) */
  GapPatternMap* gpm;           /**< Gap pattern mapping data; NULL if
                                   no indel modeling */
  int *state_to_mod;            /**< Mapping of HMM state number to tree
                                   model number */
  int *state_to_cat;            /**< Mapping of HMM state number to (spooled)
                                   category number */
  int *reverse_compl;           /**< Array of length hmm->nstates
                                   with value 1 for each state
                                   that corresponds to the reverse
                                   strand and value 0 otherwise */
  List **cat_to_states;         /**< One to many mapping */
  int *state_to_pattern;        /**< Gap pattern associated with each
                                   state, when modeling indels (-1 for
                                   no gap pattern) */
  int nratecats;                /**< Number of rate categories (for
                                   cross-product constructions) */
  int reflected;                /**< Whether "reflected" for reverse complement */
  double **emissions;           /**< Values computed by
                                   phmm_compute_emissions */
  double **forward;             /**< Forward scores */
  int alloc_len;                /**< Length for which emissions and/or
                                   forward are (or are to be) allocated */
  int *state_pos, 		/**< Contain positive tracking data for emissions */
  *state_neg;   		/**< Contain negative tracking data for emissions */
  indel_mode_type indel_mode;   /**< Indel mode in use */
  TreeNode *topology;           /**< Representative tree from tree
                                   models, used to define topology
                                   with indel model (in this case, all
                                   topologies must be the same) */
  double *alpha,		/**< Category-specific indel
                                   parameter for parametric indel
                                   model */
  *beta,			/**< Category-specific indel
                                   parameter for parametric indel
                                   model */
  *tau;   			/**< Category-specific indel
                                   parameter for parametric indel
                                   model */
  double **T,			/**< Branch-length factor used in
				   Parametric indel model */
  **t;        		        /**< Branch-length factor used in
                                   Parametric indel model */
  EmData *em_data;              /**< Used in parameter estimation by EM  */
} PhyloHmm;

/** Package of data used in estimation of indel parameters */
typedef struct {
  double *u_alpha,  
  *u_beta, 
  *u_tau;
  double **u_self, 
  **T, 				/**< Branch-length factor used in Parametric indel model */
  **fcounts;			/**< Counts of functional states, first index by state */
  int nfunctional_states, 	/**< Number of functional states */
  current_dest_cat;         	/**< Destination Category */
  GapPatternMap *gpm;		/**< Gap pattern mapping data; NULL if no indel modeling   */
  PhyloHmm *phmm;		/**< Phylo-HMM data */
} IndelEstimData;

/** \name Phylo-HMM Allocate/Cleanup functions 
\{ */

/** Create a new PhyloHmm object.
    @param hmm If indel_mode == MISSING_DATA or indel_mode ==
     PARAMETERIC, then an HMM for functional categories (with
     nstates == cm->ncats + 1) should be passed in; otherwise (indel_mode ==
     NONPARAMETERIC), an hmm for the full gapped (but not 'reflected')
     phylo-HMM should be passed in. NULL may also be given, in which
     case a trivial, single-state HMM will be created.  Any HMM that is
     passed in will be left unchanged (a copy will be created) 
    @param tree_models Array of TreeModel objects.
     Number of elements is assumed to equal number of categories
     (cm->ncats+1)
   @param cm CategoryMap.  A copy is created,
    Original isn't touched.  If NULL, a trivial map is created, with a
    separate category for every HMM state.  Must be non-NULL if
    indel_mode != MISSING_DATA. 
   @param pivot_cats Categories (by name or number)
    about which to "reflect" the HMM. Allows for prediction on both
    strands (see hmm_reverse_compl). Pass NULL for no reflection
   @param indel_mode How to model indels.  Allowable
    values are MISSING_DATA, PARAMETERIC, and NONPARAMETERIC.
 Optionally expands original HMM to
    allow for features on both the positive and negative strands. */
PhyloHmm *phmm_new(HMM *hmm, TreeModel **tree_models, CategoryMap *cm, 
                   List *pivot_cats, indel_mode_type indel_mode);
/** Free Phylo-HMM object
   @param phmm Phylo-HMM object to free
*/
void phmm_free(PhyloHmm *phmm);

/** \} */

/** Reflect PhyloHmm about pivot states corresponding to specified list of
   category names
   @param phmm Phylo-HMM to reflect
   @param pivot_cats Category names to pivot on
   @note Updates all mappings accordingly
 */
void phmm_reflect_hmm(PhyloHmm *phmm, List *pivot_cats);

/** Create an HMM representing the autocorrelation model of Felsenstein
   and Churchill
   @param[in,out] hmm HMM object already allocated with desired number of states
   @param[in] lambda Probability of entering binding state from background state  
*/
void phmm_create_autocorr_hmm(HMM *hmm, double lambda);

/** Scale tree models and replace HMM by its cross-product with an autocorr HMM defined by 'lambda'. 
   @param phmm Phylo-HMM to modify
   @param nratecats Number of rate categories
   @param lambda Probability of entering binding state from background state
   @param expand_cats Whether to expand categories
   @note Update all mappings appropriately.  
   @note Weight-matrix tree models are not scaled. 
   @see phmm_update_cross_prod
*/
void phmm_rates_cross(PhyloHmm *phmm, int nratecats, double lambda, 
                      int expand_cats);

/** Update cross product HMM using new value of lambda 
   @param phmm Cross product HMM
   @param lambda Updated value of lambda
   @see phmm_rates_cross
*/
void phmm_update_cross_prod(PhyloHmm *phmm, double lambda);



/** Compute emissions for given PhyloHmm and MSA. 
    @param phmm Initialized PhyloHMM
    @param msa Source Alignment
    @param quiet If == 1 don't report progress to stderr
    @note Used as a preprocessor for phmm_viterbi_features,
    phmm_posterior_probs, and phmm_lnl (often only needs to be run once). 
*/
void phmm_compute_emissions(PhyloHmm *phmm, MSA *msa, int quiet);

/** Calculate Log Likelihood for given Phylo-HMM and Lambda.
    @param phmm Phylo-HMM to get LogL for
    @param lambda Lambda probability
 */
double phmm_fit_lambda(PhyloHmm *phmm, double *lambda, FILE *logf);

/** Calculate predictions based on the Viterbi algorithm.
    @pre  Emissions must have already been computed
    @param phmm Initialized Phylo-HMM  
    @param seqname Sequence name for feature set (e.g., "chr1")
    @param grouptag (Optional) tag to use for groups (e.g., "exon_id", "transcript_id")
    @param idpref (Optional) prefix for assigned ids (e.g., "chr1.15")
    @param frame (Optional) Names of features for which to obtain frame 
    @see phmm_compute_emissions
*/
GFF_Set* phmm_predict_viterbi(PhyloHmm *phmm, char *seqname, char *grouptag,
                              char *idpref, List *frame);

/** Run the Viterbi algorithm and return a set of features indicating
    regions in which the Viterbi path remains in states corresponding
    to any of the specified categories.  
    @pre Emissions must have already been computed
    @param phmm Phylo-HMM object
    @param cats List of categories
    @param seqname Sequence name for feature set (e.g., "chr1")
    @param grouptag tag to use for groups (e.g., "exon_id", "transaction_id")
    @param idpref Prefix for assigned ids (e.g., "chr1.15")
    @param frame Names of features for which to obj=tain fname
    @param new_type Replace type of each retained feature with this string if
     non-NULL (old types may no longer make sense, because of merging)
    @see phmm_compute_emissions
*/
GFF_Set* phmm_predict_viterbi_cats(PhyloHmm *phmm, List *cats, char *seqname,
                                   char *grouptag, char *idpref, List *frame, 
                                   char *new_type);
/** Compute Log Likelihood with forward algorithm.
    @pre Emissions must have already been computed 
    @param phmm Phylo-HMM to compute LogL of
    @see phmm_compute_emissions)
 */
double phmm_lnl(PhyloHmm *phmm);

/** Computes posterior probabilities for a PhyloHmm. 
    @pre Emissions must have already been computed 
    @param[in] phmm PhyloHMM object
    @param[out] post_probs Calculated post probabilities
    @result Log likelihod. 
    @see phmm_compute_emissions
    @see phmm_new_postprobs
*/
double phmm_postprobs(PhyloHmm *phmm, double **post_probs);

/** Returns new object with posterior probabilities for a PhyloHmm. 
    @pre Emissions must have already been computed 
    @param[in] phmm PhyloHMM object
    @result Log likelihood. 
    @see phmm_compute_emissions
    @see phmm_postprobs
*/
double **phmm_new_postprobs(PhyloHmm *phmm);

/** Computes and returns an array of length phmm->alloc_len
    representing the marginal posterior prob at every site, summed
    over states corresponding to categories in the specified list.
    @pre Emissions must have already been computed 
    @param[in] phmm Phylo-HMM data for computation
    @param[in] cats Categories of interest by name or number
    @param[out] lnl (Optional) Log Likelihood calculated as side-effect
    @result Array representing marginal posterior probabilities
     at every site, summed over states corresponding to specified categories
    @see phmm_compute_emissions
 */
double* phmm_postprobs_cats(PhyloHmm *phmm, List *cats, double *lnl);

/** Score a set of predicted features using log odds scoring.
    @param phmm Phylo-HMM object 
    @param preds Predicted features to score
    @param score_cats Categories to score (specify by name or number)
    @param helper_cats (Optional) Secondary categories to be 
     included in scoring if adjacent to main cats.  Pass NULL if none.
    @param null_cats (Optional) Categories in null model.  Pass
     NULL to use all categories not in score_cats or helper_cats 
    @param score_only_score_cats Restrict scoring to features
     matching score_cats; if FALSE, all features are scored 
 */
void phmm_score_predictions(PhyloHmm *phmm, GFF_Set *preds, List *score_cats, 
                            List *helper_cats, List *null_cats,
                            int score_only_score_cats);

/** Add specified "bias" to log transition probabilities from
   designated background categories to non-background categories, then
   renormalize. 
   @param phmm Phylo-HMM to add bias to
   @param backgd_cat_names Categories in the background model
   @param bias Bias may be negative (increases specificity) or positive (increases sensitivity)
   @note Provides a simple "knob" for controlling the sensitivity/specificity tradeoff.  
*/
void phmm_add_bias(PhyloHmm *phmm, List *backgd_cat_names, double bias);

/* Special-purpose logging function for phmm_fit_em  
   @param logf File to save log to
   @param logl Log Likelihood that was calculated
   @param hmm HMM to log
   @param data PhyloHMM to log
   @param show_header If == 1 write a header into the log file as well
*/
void phmm_log_em(FILE *logf, double logl, HMM *hmm, void *data, 
                 int show_header);

void phmm_compute_emissions_copy_em(double **emissions, void **models, 
                                    int nmodels, void *data, int sample, 
                                    int length);
/* Wrapper for phmm_compute_emissions for use in EM 
   @param emissions NOT USED
   @param models NOT USED
   @param nmodels NOT USED
   @param data PhyloHMM object with em_data->msa
   @param sample NOT USED
   @param length NOT USED
*/
void phmm_compute_emissions_em(double **emissions, void **models, int nmodels,
                               void *data, int sample, int length);

/** Re-estimate phylogenetic models based on expected counts 
    @param models NOT USED
    @param nmodels NOT USED
    @param data Phylo-HMM object containing phylogenetic information
    @param E Category counts
    @param nobs Number of objects in each category
    @param logf File descriptor to log to
*/
void phmm_estim_mods_em(TreeModel **models, int nmodels, void *data, 
                        double **E, int nobs, FILE *logf);

/** Return observation index associated with given position, here a tuple index.
   @param data Phylo-HMM object
   @param sample NOT USED
   @param position Index of tuple_idx array
   @result Index of tuple
 */
int phmm_get_obs_idx_em(void *data, int sample, int position);

/** General routine to estimate the parameters of a phylo-HMM by EM.
   Can be used with or without the indel model, and for estimation of
   transition params only or transition parameters and tree models.
   @param phmm Phylo-HMM
   @param msa (Optional) Alignment.  NULL means not to estimate tree 
    models (emissions must be pre-computed)
   @param fix_functional Whether to fix transition parameters of functional HMM
   @param fix_indel Whether to fix indel parameters
   @param logf File to log status to
   @result Log likelihood
 */
double phmm_fit_em(PhyloHmm *phmm, MSA *msa, int fix_functional,
                   int fix_indel, FILE *logf);

/** Set up phmm->T and phmm->t, the arrays of branch lengths and branch
   length sums that are used in the parametric indel model.
   @param phmm phylo-HMM to setup branch length factors for
*/
void phmm_set_branch_len_factors(PhyloHmm *phmm);

/** Reset HMM of phylo-HMM, allowing for possible changes in indel parameters.
  @param phmm phylo-HMM to rest   */
void phmm_reset(PhyloHmm *phmm);

/** Maximize all parameters for state transition (M step of EM).
  @param hmm HMM to estimate transition for
  @param data PhyloHmm object 
  @param A Counts of functional states
  @note This function is passed to hmm_train_by_em in phmm_fit_em;
*/
void phmm_estim_trans_em(HMM *hmm, void *data, double **A);


/** \name Indel Estimate Data functions
\{ */
/** Create new indel estimation data object. 
   @param phmm Phylo-HMM object
   @param A Counts of functional states
   @result Indel Estimate Data object
*/
IndelEstimData *phmm_new_ied(PhyloHmm *phmm, double **A);

/** Free Indel Estimation Data object
   @param ied IED to free
 */
void phmm_free_ied(IndelEstimData *ied);

/** Estimate indels in an Indel Estimation Data object
  @param phmm Phylo-HMM
  @param Indel Estimate Data object
*/
void phmm_em_estim_indels(PhyloHmm *phmm, IndelEstimData *ied);

/** \} */

#endif
