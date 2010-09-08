/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: phylo_fit.h,v 1.17 2008-11-12 02:07:59 acs Exp $ */

#ifndef PHAST_CONS_H
#define PHAST_CONS_H


/* default minimum number of informative sites (see -I) */
#define DEFAULT_NSITES_THRESHOLD 50

/* default starting alpha for dgamma */
#define DEFAULT_ALPHA 1

#include <tree_model.h>
#include <msa.h>
#include <category_map.h>
#include <stringsplus.h>
#include <lists.h>
#include <gff.h>
#include "phylo_hmm.h"
#include "list_of_lists.h"

#define DEFAULT_RHO 0.3

struct phastCons_struct {
  MSA *msa;
  int post_probs, score, quiet, gff, FC, estim_lambda,
    estim_transitions, two_state, indels,
    indels_only, estim_indels, estim_trees,
    ignore_missing, estim_rho, set_transitions,
    viterbi, compute_likelihood;
  int nrates, nrates2, refidx, max_micro_indel;
  double lambda, mu, nu, alpha_0, beta_0, tau_0,
    alpha_1, beta_1, tau_1, gc, gamma, rho, omega;
  FILE *viterbi_f, *lnl_f, *log_f, *post_probs_f, *results_f, *progress_f;
  List *states, *pivot_states, *inform_reqd, *not_informative;
  TreeModel **mod;
  int nummod;
  char *seqname, *idpref, *estim_trees_fname_root, 
    *extrapolate_tree_fname;
  HMM *hmm;
  Hashtable *alias_hash;
  TreeNode *extrapolate_tree;
  CategoryMap *cm;
  ListOfLists *results;
};


int phastCons(struct phastCons_struct *p);
struct phastCons_struct* phastCons_struct_new(int rphast);
/* functions implemented below and used internally */
void setup_two_state(HMM **hmm, CategoryMap **cm, double mu, double nu);

double fit_two_state(PhyloHmm *phmm, MSA *msa, int estim_func, int estim_indels,
                     int estim_trees, int estim_rho, double *mu, double *nu, 
                     double *alpha_0, double *beta_0, double *tau_0, 
                     double *alpha_1, double *beta_1, double *tau_1, 
                     double *rho, double gamma, FILE *logf);

void reestimate_trees(void **models, int nmodels, void *data, 
                      double **E, int nobs, FILE *logf);

void reestimate_rho(void **models, int nmodels, void *data, 
		    double **E, int nobs, FILE *logf);

void phmm_estim_trans_em_coverage(HMM *hmm, void *data, double **A);

void collapse_cats(CategoryMap *cm, List *cats_to_merge);

void init_eqfreqs(TreeModel *mod, MSA *msa, double gc);

#endif
