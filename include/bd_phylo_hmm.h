/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/**
   @file bd_phylo_hmm.h
   Birth-death related phylo-HMMs.
   "birth-death" phylo-HMM -- non conserved state, fully conserved
   state, one state per birth event and one state per death event per
   branch of tree 
   @ingroup phylo_hmm
*/

#ifndef BDPHMM
#define BDPHMM

#include <tree_model.h>
#include <phylo_hmm.h>
#include <indel_mod.h>
#include <indel_history.h>
#include <gff.h>

/** Type of event that occurred i.e. birth, death, etc. */
typedef enum {CONS, /**< Constant Event */ 
	     DEATH, /**< Death Event */
	     BIRTH  /**< Birth Event */
} event_t;

/** Birth-Death Phylo-HMM variant */
typedef struct{
  PhyloHmm *phmm;               /**< phylo-HMM */
  int *state_to_branch;         /**< Mapping from states to branches on
                                   which birth/death events occur */
  int *branch_to_state_birth;	/**< Mapping from branches to states
                                   indicating birth events */
  int *branch_to_state_death;	/**< Mapping from branches to states
                                   indicating death events */
  IndelModel **indel_mods;      /**< Indel models, one per state */
  double rho, mu, nu, phi;
  int estim_gamma, /**< Estimated Gamma parameter */
  estim_omega,   /**< Estimated Omega */
  estim_phi;  /**< Estimated Phi */
} BDPhyloHmm;

/** Create a new Birth-Death phylo-HMM based on parameter values. 
    @param source_mod Initial model to build phylo-HMM from 
    @param rho Constant scaling factor
    @param mu State-transition parameter  conserved->nonconserved
    @param nu State-transition parameter  nonconserved->conserved
    @param alpha_c Conserved rate of insertion events per substitution per site
    @param beta_c Conserved rate of deletion events per substitution per site
    @param tau_c Conserved approximate inverse of expected indel length
    @param alpha_n Non-conserved rate of insertion events per substitution per site
    @param beta_n Non-conserved rate of deletion events per substitution per site
    @param tau_n Non-conserved approximate inverse of expected indel length
    @param estim_gamma Whether to estimate Gamma
    @param estim_omega Whether to estimate Omega (expected length of a conserved element)
    @param estim_phi Whether to estimate Phi
*/
BDPhyloHmm *bd_new(TreeModel *source_mod, double rho, double mu, 
                   double nu, double phi, double alpha_c, double beta_c, 
                   double tau_c, double alpha_n, double beta_n, 
                   double tau_n, int estim_gamma, int estim_omega, 
                   int estim_phi);

/** Populate transition_matrix and begin_transitions from mu,nu,phi
    @param bdphmm Birth-Death phylo-HMM to populate transition data for
 */
void bd_set_transitions(BDPhyloHmm *bdphmm);

/** Prevent birth/death events from spanning regions only supported by missing data.
    @param bdphmm Birth-Death phylo-HMM containing regions
    @param msa Alignment with sequence data
*/
void bd_handle_missing_data(BDPhyloHmm *bdphmm, MSA *msa);

/** Compute log odds scores for predictions (wrt neutral), using previously computed emissions. 
   @param bdphmm Birth-Death phylo-HMM 
   @param predictions Feature Set that is scored
*/
void bd_score_predictions(BDPhyloHmm *bdphmm, GFF_Set *predictions);

/** Combine indel emissions with substitution-based emissions 
    @param bdphmm Substitution-based emissions
    @param ih Indel emissions
*/
void bd_add_indel_emissions(BDPhyloHmm *bdphmm, IndelHistory *ih);

/** Estimate free parameters.
    @param bdphmm Birth-Death phylo-HMM to estimate transitions for
    @param msa NOT USED
 */
double bd_estimate_transitions(BDPhyloHmm *bdphmm, MSA *msa);

#endif
