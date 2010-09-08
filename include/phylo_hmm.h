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

/** \file phylo_hmm.h
   Code for phylo-HMMs.  Allows for automatic expansion of the state
   space to accommodate features on the reverse strand, and for the
     indel model described in Siepel & Haussler, RECOMB '04.  Also
   allows for cross-product constructions involving functional states
   and rate categories (Siepel & Haussler, RECOMB '03).  */

typedef enum {MISSING_DATA, PARAMETERIC, NONPARAMETERIC} indel_mode_type;

/* package of data used during EM parameter estimation */
typedef struct {
  MSA *msa;                     /**< Data to which model is fit */
  int fix_functional, fix_indel;
                                /**< Can be used to force either the
                                   functional transition probs or the
                                   indel transition probs to be fixed */
  double rho;			/**< Scaling constant for two-state
                                   rate-variation phylo-HMM */
  double gamma;			/**< Target coverage for two-state
                                   rate-variation phylo-HMM */
  Matrix *H;                /** inverse Hessian for BFGS  */
} EmData;

/** Phylo HMM object */
typedef struct {
  CategoryMap *cm;              /**< category map */
  HMM *hmm;                     /**< final HMM, after reflection and
                                   allowance for rate categories  */
  HMM *functional_hmm, *autocorr_hmm;
                                /**< original HMMs used to create cross
                                   product (NULL if no rate
                                   categories) */
  TreeModel **mods;             /**< array of tree models, after
                                   allowance for rate categories  */
  int nmods;                    /**< number of tree models (length of
                                   array mods) */
  GapPatternMap* gpm;           /**< gap pattern mapping data; NULL if
                                   no indel modeling */
  int *state_to_mod;            /**< mapping of HMM state number to tree
                                   model number */
  int *state_to_cat;            /**< mapping of HMM state number to (spooled)
                                   category number */
  int *reverse_compl;           /**< array of length hmm->nstates
                                   with value 1 for each state
                                   that corresponds to the reverse
                                   strand and value 0 otherwise */
  List **cat_to_states;         /**< one to many mapping */
  int *state_to_pattern;        /**< gap pattern associated with each
                                   state, when modeling indels (-1 for
                                   no gap pattern) */
  int nratecats;                /**< number of rate categories (for
                                   cross-product constructions) */
  int reflected;                /**< whether "reflected" for reverse compl */
  double **emissions;           /**< values computed by
                                   phmm_compute_emissions */
  double **forward;             /**< forward scores */
  int alloc_len;                /**< length for which emissions and/or
                                   forward are (or are to be) allocated */
  int *state_pos, *state_neg;   /**< contain tracking data for emissions */
  indel_mode_type indel_mode;   /**< Indel mode in use */
  TreeNode *topology;           /**< representative tree from tree
                                   models, used to define topology
                                   with indel model (in this case, all
                                   topologies must be the same) */
  double *alpha, *beta, *tau; /**< category-specific indel
                                   parameters for parameteric indel
                                   model */
  double **T, **t;              /**< branch-length factors used in
                                   parameteric indel model */
  EmData *em_data;              /**< Used in parameter estimation by EM  */
} PhyloHmm;

/* package of data used in estimation of indel parameters */
typedef struct {
  double *u_alpha, *u_beta, *u_tau;
  double **u_self, **T, **fcounts;
  int nfunctional_states, current_dest_cat;         
  GapPatternMap *gpm;
  PhyloHmm *phmm;
} IndelEstimData;

PhyloHmm *phmm_new(HMM *hmm, TreeModel **tree_models, CategoryMap *cm, 
                   List *pivot_cats, indel_mode_type indel_mode);
void phmm_reflect_hmm(PhyloHmm *phmm, List *pivot_cats);
void phmm_create_autocorr_hmm(HMM *hmm, double lambda);
void phmm_rates_cross(PhyloHmm *phmm, int nratecats, double lambda, 
                      int expand_cats);
void phmm_update_cross_prod(PhyloHmm *phmm, double lambda);
void phmm_free(PhyloHmm *phmm);
void phmm_compute_emissions(PhyloHmm *phmm, MSA *msa, int quiet);
double phmm_fit_lambda(PhyloHmm *phmm, double *lambda, FILE *logf);
void phmm_update_cross_prod(PhyloHmm *phmm, double lambda);
GFF_Set* phmm_predict_viterbi(PhyloHmm *phmm, char *seqname, char *grouptag,
                              char *idpref, List *frame);
GFF_Set* phmm_predict_viterbi_cats(PhyloHmm *phmm, List *cats, char *seqname,
                                   char *grouptag, char *idpref, List *frame, 
                                   char *new_type);
double phmm_lnl(PhyloHmm *phmm);
double phmm_postprobs(PhyloHmm *phmm, double **post_probs);
double* phmm_postprobs_cats(PhyloHmm *phmm, List *cats, double *lnl);
void phmm_score_predictions(PhyloHmm *phmm, GFF_Set *preds, List *score_cats, 
                            List *helper_cats, List *null_cats,
                            int score_only_score_cats);
void phmm_add_bias(PhyloHmm *phmm, List *backgd_cat_names, double bias);
void phmm_log_em(FILE *logf, double logl, HMM *hmm, void *data, 
                 int show_header);
void phmm_compute_emissions_copy_em(double **emissions, void **models, 
                                    int nmodels, void *data, int sample, 
                                    int length);
void phmm_compute_emissions_em(double **emissions, void **models, int nmodels,
                               void *data, int sample, int length);
void phmm_estim_mods_em(void **models, int nmodels, void *data, 
                        double **E, int nobs, FILE *logf);
int phmm_get_obs_idx_em(void *data, int sample, int position);
double phmm_fit_em(PhyloHmm *phmm, MSA *msa, int fix_functional,
                   int fix_indel, FILE *logf);
void phmm_set_branch_len_factors(PhyloHmm *phmm);
void phmm_reset(PhyloHmm *phmm);
void phmm_estim_trans_em(HMM *hmm, void *data, double **A);
IndelEstimData *phmm_new_ied(PhyloHmm *phmm, double **A);
void phmm_free_ied(IndelEstimData *ied);
void phmm_em_estim_indels(PhyloHmm *phmm, IndelEstimData *ied);


#endif
