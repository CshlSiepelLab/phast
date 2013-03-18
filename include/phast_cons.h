/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file phast_cons.h
  Functions to support the program phastCons.
*/

#ifndef PHAST_CONS_H
#define PHAST_CONS_H


/** Default minimum number of informative sites (see -I) */
#define DEFAULT_NSITES_THRESHOLD 50

/** Default starting alpha for dgamma */
#define DEFAULT_ALPHA 1

#include <tree_model.h>
#include <msa.h>
#include <category_map.h>
#include <stringsplus.h>
#include <lists.h>
#include <gff.h>
#include "phylo_hmm.h"
#include "list_of_lists.h"

/** Default RHO */
#define DEFAULT_RHO 0.3

/** Package holding all phastCons data */
struct phastCons_struct {
  MSA *msa;		/**< Multiple Sequence Alignment */
  int post_probs,	/**< Whether to use posterior probabilities */
    score,		/**< Whether to calculate a score */
    quiet,		/**< Whether to display warnings/errors on stderr */
    gff,		/**< Feature Set */
    FC,			/**< Whether to use Felsenstein/Churchill model */
    estim_lambda,	/**< Whether to estimate lambda via ML or have it provided by user */
    estim_transitions,  /**< Whether to estimate transitions or have it provided by user */
    two_state,		/**< ???*/
    indels,		/**< Whether to expand HMM state space to model indels as described in Siepel & Haussler (2004)*/
    indels_only,	/**< Whether to expand HMM state space to model indels in a single HMM state */
    estim_indels,	/**< Whether to estimate indels */
    estim_trees,	/**< Whether to estimate free parameters of tree models */
    ignore_missing,	/**< Whether to ignore regions of missing data in all sequences but the reference
			   sequence (excluding sequences specified by --not-informative) when 
			   estimating transition probabilities */
    estim_rho,		/**< Whether to estimate the rho parameter */
    set_transitions,	/**< Whether user supplies mu, nu for transition information, otherwise estimated */
    viterbi,		/**< Whether to use Viterbi algorithm to predict discrete elements */
    compute_likelihood; /**< Whether to compute the likelihood */
  int nrates,		/**< Number of rates for first tree model */
    nrates2,		/**< Number of rates for second tree model */
    refidx,		/**< Index of reference sequence */
    max_micro_indel;	/**< Maximum length of an alignment gap, any gap longer is treated as missing data*/
  double lambda,	/**< Lambda parameter value */ 
    mu,			/**< Transitions mu value */
    nu,			/**< Transitions nu value */
    alpha_0,		/**< Rate of insertion events per substitution per site in Conserved state */
    beta_0,		/**< Rate of deletion events per substitution per site in Conserved state */
    tau_0,		/**< Approximately the inverse of the expected indel length in Conserved state */
    alpha_1,		/**< Rate of insertion events per substitution per site in Non-Conserved state */
    beta_1,		/**< Rate of deletion events per substitution per site in Non-Conserved state */
    tau_1,		/**< Approximately the inverse of the expected indel length in Non-Conserved state */
    gc,			/**< G+C Content value */
    gamma,		/**< Gamma parameter value */
    rho,		/**< Rho parameter value */
    omega;		/**< Expected length of a conserved element Omega */
  FILE *viterbi_f,	/**< File descriptor for printing out predictions */
    *lnl_f,		/**< File descriptor to save log likelihood results and parameters used to calculate it */
    *log_f,		/**< File descriptor to save general info */
    *post_probs_f,	/**< File descriptor to save posterior probs */
    *results_f,		/**< File descriptor to save results */
    *progress_f;	/**< File descriptor to save progress */
  List *states,		/**< List of states of interest in the phylo-HMM, specified by number (starts at 0), or if --catmap, by category name */
    *pivot_states,	/**< List of "pivot" states to "reflect" forward strand HMM around specified by number (starts at 0), or if --catmap, by category name*/
    *inform_reqd,	/**< List of states that must have "informative" columns (i.e., columns with more than two non-missing-data characters) use "none" to disable  */
    *not_informative;	/**< List of sequences not to consider when deciding if a column is "informative" or not */ 
  TreeModel **mod;	/**< List of tree models for MSAs */
  int nummod;		/**< Number of tree models for MSAs */
  char *seqname,	/**< Sequence name for reference sequence */
    *idpref,		/**< Prefix for assigned ids */
    *estim_trees_fname_root,	/**< Root part of filename for tree models i.e. %s.cons.mod or %s.noncons.mod */
    *extrapolate_tree_fname,	/**< Filepath to tree file used to extrapolate a larger set of species*/
    *bgc_branch;        /**< If not NULL, assume a two-state HMM with and without bgc on the named branch*/
  HMM *hmm;		       /**< Hidden Markov Model */
  Hashtable *alias_hash;       /**< Sequence name aliases e.g., "hg17=human; mm5=mouse; rn3=rat" */
  TreeNode *extrapolate_tree;	/**< Root of tree used for extrapolation of larget set of species */
  CategoryMap *cm;		/**< Category Map */
  ListOfLists *results;		/**< Holds results of phast_cons analyses */
};

/** \name Main phastCons functions 
\{ */
/** Run phastCons analysis.
   @param p All settings and preferences needed to perform analysis
   @result 0 on success
 */
int phastCons(struct phastCons_struct *p);

/** Create a new object to hold all settings/data for analysis.
  @param rphast If this code is being used in rphast == 1
  @result Newly allocated phastCons settings/data object
 */
struct phastCons_struct* phastCons_struct_new(int rphast);

/* Functions implemented below and used internally */
/** \} \name Supporting functions 
\{ */

/** Set up HMM and category map for two-state case 
  @param[out] hmm HMM to setup
  @param[out] cm Category Map to setup
  @param[in] mu MU transition value
  @param[in] nu NU transition value
*/
void setup_two_state(HMM **hmm, CategoryMap **cm, double mu, double nu);

/** Estimate parameters for the two-state model using an EM algorithm.
   Any or all of the parameters 'mu' and 'nu', the indel parameters, and
   the tree models themselves may be estimated.  
   @param phmm Allocated phylo-HMM
   @param msa Sequence Alignment data
   @param estim_func Whether to estimate functional
   @param estim_indels Whether to estimate indels
   @param estim_trees Whether to estimate trees
   @param estim_rho Whether to estimate rho
   @param estim_trees_partial Whether this is an HMM where we want to estimate parameters which are shared across tree models
   @param mu (Optional) Transitions value
   @param nu (Optional) Transitions value
   @param alpha_0 (Optional) Rate of insertion events per substitution per site in Conserved state 
   @param beta_0 (Optional) Rate of deletion events per substitution per site in Conserved state 
   @param tau_0 (Optional) Approximately the inverse of the expected indel length in Conserved state 
   @param alpha_1 (Optional) Rate of insertion events per substitution per site in Non-Conserved state 
   @param beta_1 (Optional) Rate of deletion events per substitution per site in Non-Conserved state 
   @param tau_1 (Optional) Approximately the inverse of the expected indel length in Non-Conserved state 
   @param rho (Optional) Rho parameter value
   @param gamma Gamma parameter
   @param logf File descriptor of where to save log info 
   @result Log Likelihood. */
double fit_two_state(PhyloHmm *phmm, MSA *msa, int estim_func, int estim_indels,
                     int estim_trees, int estim_rho,
		     double *mu, double *nu, 
                     double *alpha_0, double *beta_0, double *tau_0, 
                     double *alpha_1, double *beta_1, double *tau_1, 
                     double *rho, double gamma, FILE *logf);

/** Re-estimate phylogenetic model based on expected counts (M step of EM) 
   @param models NOT USED
   @param nmodels NOT USED 
   @param PhyloHmm object 
   @param E Category Counts 
   @param nobs Number of objects in E[0]
   @param logf File descriptor of where to save logs
*/
void reestimate_trees(TreeModel **models, int nmodels, void *data, 
                      double **E, int nobs, FILE *logf);

/** re-estimate only scale parameter rho (only the first model for the conserved state needs to be considered)   
   @param models NOT USED
   @param nmodels NOT USED 
   @param PhyloHmm object 
   @param E Category Counts 
   @param nobs Number of objects in E[0]
   @param logf File descriptor of where to save logs
*/
void reestimate_rho(TreeModel **models, int nmodels, void *data, 
		    double **E, int nobs, FILE *logf);

/** Maximize HMM transition parameters subject to constrain implied by
   target coverage (M step of EM).  
   @param hmm NOT USED
   @param data PhyloHmm object
   @param A Counts of functional states
   @note For use with two-state HMM.  
   @note This function is passed to hmm_train_by_em in phmm_fit_em
*/
void phmm_estim_trans_em_coverage(HMM *hmm, void *data, double **A);

void collapse_cats(CategoryMap *cm, List *cats_to_merge);

/** Initialize equilibrium frequencies for tree model
    @param mod Tree Model to initialize equilibrium frequencies for
    @param msa Sequence Alignment data to use for estimation
    @param gc (Optional) If gc == -1 estimate from alignment, otherwise specifies G+C content   
*/
void init_eqfreqs(TreeModel *mod, MSA *msa, double gc);

/** \} */
#endif
