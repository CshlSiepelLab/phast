/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/**
   @file bgc_hmm.h
   Functions to set up and estimate parameters for a Hidden Markov Model
   representing a  GC-biased gene conversion process.
   @ingroup hmm
*/

#ifndef BGC_HMM_H
#define BGC_HMM_H

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include "misc.h"
#include "tree_model.h"
#include "phast_cons.h"
#include "maf.h"
#include "gff.h"
#include "sufficient_stats.h"
#include "em.h"
#include "tree_likelihoods.h"
#include "list_of_lists.h"

typedef enum {NONE, WIG, FULL} bgchmm_posterior_output_type;


/** Holds all the options for how to set up and estimate parameters in the
    HMM.  Upon return from bgcHmm also may contain the results
 */
struct bgchmm_struct {
  MSA *msa;
  TreeModel *mod;
  double scale;   //overall scale of input neutral model
  double rho;     // conservation tree scale
  double cons_expected_length;
  double cons_target_coverage;
  double bgc;  //bgc parameter B
  double sel;  //sel parameter in conserved state
  double bgc_target_coverage;
  double bgc_expected_length;
  char *foregd_branch;
  char *viterbi_fn;
  char *tract_fn;
  char *informative_fn;
  char *mods_fn;
  int estimate_bgc;
  int estimate_scale;
  int estimate_bgc_target_coverage;
  int estimate_bgc_expected_length;
  int do_bgc;
  int estimate_cons_transitions;
  int null_model;
  int estimate_rho;
  int eqfreqs_from_msa;
  bgchmm_posterior_output_type post_probs;
  int informative_only; // this option implies do nothing except determine which sites are informative for gBGC
  int random_path;
  int get_likelihoods;
  FILE *post_probs_f;
  ListOfLists *results;
};


struct bgchmm_data_struct {
  MSA *msa;
  int *bgc_informative;
  int nsite;
  int *msa_ixd;
  int estimate_cons_transitions;
  int estimate_bgc_target_coverage, estimate_bgc_expected_length;
  double bgc_target_coverage, bgc_expected_length;
  HMM *hmm;
};


int bgcHmm(struct bgchmm_struct *b);

struct bgchmm_struct *bgchmm_struct_new(int rphast);

TreeModel **bgchmm_setup_mods(TreeModel *init_mod, char *foregd_branch,
			      int do_bgc, double bgc, double sel, double rho, 
			      double scale, int estimate_bgc,
			      int estimate_rho, int estimate_scale,
			      int eqfreqs_from_msa, MSA *align, int *npar);

void bgchmm_compute_emissions(double **emissions, void **models, int nmodels,
			      void *data, int sample, int length);

int bgchmm_get_obs_idx(void *data, int i, int j);

void bgchmm_set_hmm(HMM *hmm, double bgc_in, double bgc_out, double cons_in, double cons_out);

void bgchmm_estimate_transitions(HMM* hmm, void * data, double **A);

void bgchmm_get_rates(HMM *hmm, double *bgc_in, double *bgc_out, double *cons_in, double *cons_out);

void bgchmm_estimate_states(TreeModel **mods, int nmod, void *data, double **E, int nobs, FILE *logfile);

int *bgchmm_get_informative(MSA *msa, char *foregd, TreeNode *tree);

void bgchmm_print_informative(MSA *msa, int *bgc_informative, ListOfLists *results, char *informative_fn, int reverse);

char *bgchmm_get_state_name(int state, int do_bgc);

void bgchmm_output_path(int *path, int nsite, MSA *msa, int do_bgc,
         		char *name, char *outfn, ListOfLists *results);

#endif
