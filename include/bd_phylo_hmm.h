/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* "birth-death" phylo-HMM -- nonconserved state, fully conserved
   state, one state per birth event and one state per death event per
   branch of tree */

#ifndef BDPHMM
#define BDPHMM

#include <tree_model.h>
#include <phylo_hmm.h>
#include <indel_mod.h>
#include <indel_history.h>
#include <gff.h>

typedef enum {CONS, DEATH, BIRTH} event_t;

typedef struct{
  PhyloHmm *phmm;               /* phylo-HMM */
  int *state_to_branch;         /* mapping from states to branches on
                                   which birth/death events occur */
  int *branch_to_state_birth;	/* mapping from branches to states
                                   indicating birth events */
  int *branch_to_state_death;	/* mapping from branches to states
                                   indicating death events */
  IndelModel **indel_mods;      /* indel models, one per state */
  double rho, mu, nu, phi;
  int estim_gamma, estim_omega, estim_phi;
} BDPhyloHmm;


BDPhyloHmm *bd_new(TreeModel *source_mod, double rho, double mu, 
                   double nu, double phi, double alpha_c, double beta_c, 
                   double tau_c, double alpha_n, double beta_n, 
                   double tau_n, int estim_gamma, int estim_omega, 
                   int estim_phi);
void bd_set_transitions(BDPhyloHmm *bdphmm);
void bd_handle_missing_data(BDPhyloHmm *bdphmm, MSA *msa);
void bd_score_predictions(BDPhyloHmm *bdphmm, GFF_Set *predictions);
void bd_add_indel_emissions(BDPhyloHmm *bdphmm, IndelHistory *ih);
double bd_estimate_transitions(BDPhyloHmm *bdphmm, MSA *msa);

#endif
