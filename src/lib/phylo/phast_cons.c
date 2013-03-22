/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <phylo_hmm.h>
#include <em.h>
#include <sufficient_stats.h>
#include <bed.h>
#include <dgamma.h>
#include <tree_likelihoods.h>
#include <maf.h>
#include "phast_cons.h"


struct phastCons_struct *phastCons_struct_new(int rphast) {
  struct phastCons_struct *p = smalloc(sizeof(struct phastCons_struct));
  p->post_probs = TRUE; 
  p->score = FALSE; 
  p->gff = FALSE; 
  p->FC = FALSE; 
  p->estim_lambda = TRUE;
  p->estim_transitions = TRUE; 
  p->two_state = TRUE; 
  p->indels = FALSE;
  p->indels_only = FALSE; 
  p->estim_indels = TRUE;
  p->estim_trees = FALSE; 
  p->ignore_missing = FALSE; 
  p->estim_rho = FALSE;
  p->set_transitions = FALSE;
  p->nrates = -1;
  p->nrates2 = -1;
  p->refidx = 1; 
  p->max_micro_indel = 20;
  p->lambda = 0.9; 
  p->mu = 0.01; 
  p->nu = 0.01; 
  p->alpha_0 = 0.05; 
  p->beta_0 = 0.05;
  p->tau_0 = 0.45; 
  p->alpha_1 = 0.05; 
  p->beta_1 = 0.05; 
  p->tau_1 = 0.2;
  p->gc = -1;
  p->gamma = -1; 
  p->rho = DEFAULT_RHO; 
  p->omega = -1;
  p->viterbi_f = NULL;
  p->viterbi = FALSE;
  p->lnl_f = NULL; 
  p->log_f = NULL;
  p->states = NULL; 
  p->pivot_states = NULL; 
  p->inform_reqd = NULL;
  p->mod=NULL;
  p->nummod=0;
  p->not_informative = NULL;
  p->seqname = NULL;
  p->idpref = NULL; 
  p->estim_trees_fname_root = NULL;
  p->extrapolate_tree_fname = NULL;
  p->hmm = NULL;
  p->alias_hash = NULL;
  p->extrapolate_tree = NULL;
  p->cm = NULL;
  p->compute_likelihood = FALSE;
  p->post_probs_f = rphast ? NULL : stdout;
  p->results_f = rphast ? stdout : stderr;
  p->progress_f = rphast ? stdout : stderr;
  p->results = rphast ? lol_new(2) : NULL;
  return p;
}


int phastCons(struct phastCons_struct *p) {
  int post_probs, score, quiet, gff, FC, estim_lambda,
    estim_transitions, two_state, indels, 
    indels_only, estim_indels,
    estim_trees, ignore_missing, estim_rho, set_transitions,
    nummod, viterbi, compute_likelihood;
  int nrates, nrates2, refidx, max_micro_indel, free_cm=0;
  double lambda, mu, nu, alpha_0, beta_0, tau_0, alpha_1, beta_1, tau_1,
    gc, gamma, rho, omega;
  FILE *viterbi_f, *lnl_f, *log_f, *post_probs_f, *results_f;
  List *states, *pivot_states, *inform_reqd,
    *not_informative;
  char *seqname, *idpref, *estim_trees_fname_root,
    *extrapolate_tree_fname;
  HMM *hmm;
  Hashtable *alias_hash;
  TreeNode *extrapolate_tree;
  CategoryMap *cm;
  TreeModel **mod;
  ListOfLists *results;
  MSA *msa;

  /* other vars */
  int i, j, last;
  double lnl = INFTY;
  PhyloHmm *phmm;
  char *newname;
  indel_mode_type indel_mode;

  msa = p->msa;
  post_probs = p->post_probs;
  score = p->score;
  gff = p->gff;
  FC = p->FC;
  estim_lambda = p->estim_lambda;
  estim_transitions = p->estim_transitions;
  two_state = p->two_state; 
  indels = p->indels;
  indels_only = p->indels_only;
  estim_indels = p->estim_indels;
  estim_trees = p->estim_trees;
  ignore_missing = p->ignore_missing;
  estim_rho = p->estim_rho;
  set_transitions = p->set_transitions;
  nrates = p->nrates;
  nrates2 = p->nrates2; 
  refidx = p->refidx;
  max_micro_indel = p->max_micro_indel;
  lambda = p->lambda;
  mu = p->mu;
  nu = p->nu;
  alpha_0 = p->alpha_0;
  beta_0 = p->beta_0;
  tau_0 = p->tau_0;
  alpha_1 = p->alpha_1;
  beta_1 = p->beta_1;
  tau_1 = p->tau_1;
  gc = p->gc;
  gamma = p->gamma;
  rho = p->rho;
  omega = p->omega;
  viterbi_f = p->viterbi_f;
  lnl_f = p->lnl_f;
  log_f = p->log_f;
  states = p->states;
  pivot_states = p->pivot_states;
  inform_reqd = p->inform_reqd; 
  not_informative = p->not_informative;
  mod = p->mod;
  nummod = p->nummod;
  seqname = p->seqname;
  idpref = p->idpref;
  estim_trees_fname_root = p->estim_trees_fname_root;
  extrapolate_tree_fname = p->extrapolate_tree_fname;
  hmm = p->hmm;
  alias_hash = p->alias_hash;
  extrapolate_tree = p->extrapolate_tree;
  cm = p->cm;
  post_probs_f = p->post_probs_f;
  results_f = p->results_f;
  results = p->results;
  viterbi = p->viterbi;
  if (lnl_f != NULL) compute_likelihood=TRUE;
  else if (lnl_f == NULL && results==NULL) compute_likelihood=FALSE;
  else compute_likelihood = p->compute_likelihood;
  if (viterbi_f != NULL) viterbi=TRUE;
  quiet = (results_f == NULL);

  /* enforce usage rules */
  if ((hmm != NULL && FC))
    die("ERROR: --hmm and --FC are mutually exclusive.\n");

  if (indels_only && (hmm != NULL || FC))
    die("ERROR: --indels-only cannot be used with --hmm or --FC.\n");

  if ((estim_trees || gamma != -1 || estim_rho || 
       omega != -1 || set_transitions) && !two_state)
    die("ERROR: --estimate-trees, --target-coverage, --expected-length, --transitions,\nand --estimate-rho can only be used with default two-state HMM.\n");

  if (set_transitions && (gamma != -1 || omega != -1))
    die("ERROR: --transitions and --target-coverage/--expected-length cannot be used together.\n");

  if (omega != -1 && gamma == -1) 
    die("ERROR: --expected-length requires --target-coverage.\n");

  if (cm != NULL && hmm == NULL) 
    die("ERROR: --catmap can only be used with --hmm.\n");
  
  if (indels == TRUE && FC)
    die("ERROR: --indels cannot be used with --FC.\n");

  if (nrates != -1 && hmm != NULL)
    die("ERROR: --nrates currently can't be used with --hmm.\n");

  if (!indels) estim_indels = FALSE;

  if (msa_alph_has_lowercase(msa)) msa_toupper(msa); 
  msa_remove_N_from_alph(msa);  /* for backward compatibility */
  if (msa->ss == NULL) 
    ss_from_msas(msa, nummod==0 ? 1 : mod[0]->order+1, 
		 TRUE, NULL, NULL, NULL, -1, 
		 nummod == 0 ? 0 : subst_mod_is_codon_model(mod[0]->subst_mod));
  if (msa->ss->tuple_idx == NULL)
    die("ERROR: Ordered representation of alignment required.\n");
                                /* SS assumed below */

  /* rename if aliases are defined */
  if (alias_hash != NULL) {
    for (i = 0; i < msa->nseqs; i++) {
      if ((newname = hsh_get(alias_hash, msa->names[i])) != (char*)-1) {
        sfree(msa->names[i]);
        msa->names[i] = copy_charstr(newname);
      }
    }
  }

  /* mask out macro-indels, if necessary */
  if (indels) {
    /* this little hack allows gaps in refseq to be restored before
       output (needed for proper coord conversion) */
    if (msa->seqs == NULL) { ss_to_msa(msa); ss_free(msa->ss); msa->ss = NULL; }
    if (strlen(msa->missing) < 2)
      die("ERROR strlen(msa->missing)=%i\n", strlen(msa->missing));
    for (i = 0; i < msa->length; i++) 
      if (msa->is_missing[(int)msa->seqs[0][i]]) msa->seqs[0][i] = msa->missing[1];
                                /* msa->missing[0] is used in msa_mask_macro_indels */
    msa_mask_macro_indels(msa, max_micro_indel, 0);
  }

  /* Set up array indicating which seqs are informative, if necessary */
  if (not_informative != NULL)
    msa_set_informative(msa, not_informative);

  /* strip missing columns, if necessary */
  if (ignore_missing)
    ss_strip_missing(msa, refidx);

  if ((FC || indels_only) && nummod != 1)
    die("ERROR: only one tree model allowed with --FC and --indels-only.\n");

  if (two_state && nummod > 2)
    die("ERROR: must specify either one or two tree models with default two-state model.\n");

  /* prune/extrapolate tree models, check refidx */
  for (i = 0; i < nummod; i++) {
    int old_nnodes, found;
    List *pruned_names = lst_new_ptr(msa->nseqs);

    old_nnodes = mod[i]->tree->nnodes;

    /* extrapolate tree and/or prune away extra species */
    if (extrapolate_tree != NULL) {
      double scale = tm_extrapolate_and_prune(mod[i], extrapolate_tree, 
                                              msa, pruned_names);
      if (!quiet) 
        fprintf(results_f, "Extrapolating based on %s (scale=%f)...\n", 
                extrapolate_tree_fname, scale);
    }
    else
      tm_prune(mod[i], msa, pruned_names);

    if (lst_size(pruned_names) == (old_nnodes + 1) / 2)
      die("ERROR: no match for leaves of tree in alignment (leaf names must match alignment names).\n");
    if (!quiet && lst_size(pruned_names) > 0) {
      fprintf(results_f, "WARNING: pruned away leaves of tree with no match in alignment (");
      for (j = 0; j < lst_size(pruned_names); j++)
        fprintf(results_f, "%s%s", ((String*)lst_get_ptr(pruned_names, j))->chars, 
                j < lst_size(pruned_names) - 1 ? ", " : ").\n");
    }

    /* also make sure match for reference sequence in tree */
    if (refidx > 0) {
      for (j = 0, found = FALSE; !found && j < mod[i]->tree->nnodes; j++) {
	TreeNode *n = lst_get_ptr(mod[i]->tree->nodes, j);
	if (n->lchild == NULL && n->rchild == NULL && 
	    n->name != NULL && !strcmp(n->name, msa->names[refidx-1]))
	  found = TRUE;
      }
      if (!found) die("ERROR: no match for reference sequence in tree.\n");
    }

    lst_free_strings(pruned_names);
    lst_free(pruned_names);
  }

  /* initial checks and setup of tree models for two-state and FC */
  if (two_state) {
    if (nummod == 2 && (estim_trees || estim_rho))
      die("ERROR: If re-estimating tree models, pass in only one model for initialization.\n");
    if (mod[0]->empirical_rates || 
	(nummod == 2 && mod[1]->empirical_rates))
      die("ERROR: nonparameteric rate variation not allowed with default two-state HMM.\n");

    /* set equilibrium frequencies if estimating tree models */
    if (estim_trees || (gc != -1 && estim_rho)) 
      init_eqfreqs(mod[0], msa, gc);
    
    if (nummod == 1) { /* create 2nd tree model &
			  rescale first */
      mod = srealloc(mod, 2 * sizeof(void*));
      mod[1] = tm_create_copy(mod[0]);
      if (!estim_rho) tm_scale_branchlens(mod[0], rho, TRUE);
    }
    if (nrates != -1 && nrates != mod[0]->nratecats) 
      tm_reinit(mod[0], mod[0]->subst_mod, nrates, mod[0]->alpha, NULL, NULL);
    if (nrates2 != -1 && nrates2 != mod[1]->nratecats) 
      tm_reinit(mod[1], mod[1]->subst_mod, nrates2, mod[1]->alpha, NULL, NULL);
  }
  else if (FC) {
    if (mod[0]->nratecats <= 1)
      die("ERROR: a tree model allowing for rate variation is required.\n");
    if (nrates != -1 && mod[0]->empirical_rates)
      die("ERROR: can't use --nrates with nonparameteric rate model.\n");
    if (nrates == -1) nrates = mod[0]->nratecats;
  }

  if (viterbi && seqname==NULL)
    seqname = "refseq";

  /* set up states */
  if (states == NULL && (two_state || results==NULL)) {
    states = lst_new_ptr(1);
    lst_push_ptr(states, str_new_charstr("0"));
  }

  /* set require-informative to states if null; set to null if "none" */ 
  if (inform_reqd == NULL)
    inform_reqd = states;
  else if (lst_size(inform_reqd) == 1 && 
	   str_equals_charstr(lst_get_ptr(inform_reqd, 0), "none")) 
    inform_reqd = NULL;

  if (two_state) {
    if (!quiet) 
      fprintf(results_f, "Creating 'conserved' and 'nonconserved' states in HMM...\n");
    if (gamma != -1) {
      nu = gamma/(1-gamma) * mu;
      if (nu >= 1) 
        die("ERROR: mu=%f and gamma=%f imply nu >= 1.\n", mu, gamma);
    }
    setup_two_state(&hmm, &cm, mu, nu);
  }
  else if (cm == NULL) {
    cm = cm_create_trivial(nummod-1, NULL);
    free_cm=TRUE;
  }

  /* set up PhyloHmm */
  if (!indels) indel_mode = MISSING_DATA;
  else if (hmm == NULL || hmm->nstates == cm->ncats + 1)
    indel_mode = PARAMETERIC;
  else indel_mode = NONPARAMETERIC;

  phmm = phmm_new(hmm, mod, cm, pivot_states, indel_mode);

  if (FC) {
    if (!quiet) 
      fprintf(results_f, "Creating %d scaled versions of tree model...\n", nrates);
    phmm_rates_cross(phmm, nrates, lambda, TRUE);
  }

  /* set inform_reqd, if necessary.  This has to be done
     *after* the set of models is expanded (two-state or FC) */
  if (inform_reqd != NULL) {
    List *l = cm_get_category_list(cm, inform_reqd, 0);
    for (i = 0; i < lst_size(l); i++) {
      int modno = lst_get_int(l, i);
      if (modno < 0 || modno >= phmm->nmods) 
        die("ERROR: illegal argument to --require-informative.\n");
      phmm->mods[modno]->inform_reqd = TRUE;
    }
    lst_free(l);
  }        
  if (free_cm) cm_free(cm);

  /* compute emissions */
  phmm_compute_emissions(phmm, msa, quiet);

  /* estimate lambda, if necessary */
  if (FC && estim_lambda) {
    if (!quiet) fprintf(results_f, "Finding MLE for lambda...");
    lnl = phmm_fit_lambda(phmm, &lambda, log_f);
    if (!quiet) fprintf(results_f, " (lambda = %g)\n", lambda);
    if (results != NULL) lol_push_dbl(results, &lambda, 1, "lambda");
    phmm_update_cross_prod(phmm, lambda);
  }

  /* estimate mu and nu and indel params, if necessary */
  else if (two_state && 
	   (estim_transitions || estim_indels || estim_trees || estim_rho)) {
    char cons_fname[STR_MED_LEN], noncons_fname[STR_MED_LEN];
    if (!quiet) {
      fprintf(results_f, "Finding MLE for (");
      if (estim_transitions) 
        fprintf(results_f, "mu, nu%s", estim_indels || estim_trees || estim_rho 
		? ", " : "");
      if (estim_indels) 
        fprintf(results_f, "alpha_0, beta_0, tau_0, alpha_1, beta_1, tau_1%s",
                estim_trees || estim_rho ? ", " : "");
      if (estim_trees)
        fprintf(results_f, "[tree models]");
      else if (estim_rho) 
        fprintf(results_f, "rho");
      fprintf(results_f, ")...\n");
    }
    lnl = fit_two_state(phmm, msa, estim_transitions, estim_indels, 
			estim_trees, estim_rho,
                        &mu, &nu, &alpha_0, &beta_0, &tau_0, 
                        &alpha_1, &beta_1, &tau_1, &rho,
                        gamma, log_f);
    if (estim_transitions || estim_indels || estim_rho) {
      if (!quiet) {
	fprintf(results_f, "(");
	if (estim_transitions)
	  fprintf(results_f, "mu = %g. nu = %g%s", mu, nu, 
		  estim_indels || estim_rho ? ", " : "");
	if (estim_indels)
	  fprintf(results_f, 
		  "alpha_0 = %g, beta_0 = %g, tau_0 = %g, alpha_1 = %g, beta_1 = %g, tau_1 = %g%s", 
		  alpha_0, beta_0, tau_0, alpha_1, beta_1, tau_1, 
		  estim_rho ? ", " : "");
	if (estim_rho) 
	  fprintf(results_f, "rho = %g", rho);
	fprintf(results_f, ")\n");
      } 
      if (results != NULL) {
	double *temp;
	if (estim_transitions) {
	  temp = smalloc(2*sizeof(double));
	  temp[0] = mu;
	  temp[1] = nu;
	  lol_push_dbl(results, temp, 2, "transition.rates");
	  sfree(temp);
	}
	if (estim_indels) {
	  temp = smalloc(6*sizeof(double));
	  temp[0] = alpha_0;
	  temp[1] = beta_0;
	  temp[2] = tau_0;
	  temp[3] = alpha_1;
	  temp[4] = beta_1;
	  temp[5] = tau_1;
	  lol_push_dbl(results, temp, 6, "indel.rates");
	  sfree(temp);
	}
	if (estim_rho)
	  lol_push_dbl(results, &rho, 1, "rho");
      }
    }
    if (estim_trees || estim_rho) {
      if (estim_trees_fname_root != NULL) {
	sprintf(cons_fname, "%s.cons.mod", estim_trees_fname_root);
	sprintf(noncons_fname, "%s.noncons.mod", estim_trees_fname_root);
	if (!quiet)
	  fprintf(results_f, "Writing re-estimated tree models to %s and %s...\n", 
		  cons_fname, noncons_fname);
	tm_print(phast_fopen(cons_fname, "w+"), phmm->mods[0]);
	tm_print(phast_fopen(noncons_fname, "w+"), phmm->mods[1]);
      }
      if (results != NULL) {
	ListOfLists *tmplist = lol_new(2);
	lol_push_treeModel(tmplist, phmm->mods[0], "cons.mod");
	lol_push_treeModel(tmplist, phmm->mods[1], "noncons.mod");
	lol_push_lol(results, tmplist, "tree.models");
      }
    }
  }

  /* estimate indel parameters only, if necessary */
  else if (indels_only) {
    if (!quiet) fprintf(results_f, "Estimating parameters for indel model...");
    lnl = phmm_fit_em(phmm, msa, TRUE, FALSE, log_f);
    if (!quiet) fprintf(results_f, "...\n");
  }

  /* still have to set indel params if not estimating */
  else if (indel_mode == PARAMETERIC) {
    phmm->alpha[0] = alpha_0; phmm->beta[0] = beta_0; phmm->tau[0] = tau_0;
    phmm->alpha[1] = alpha_1; phmm->beta[1] = beta_1; phmm->tau[1] = tau_1;
    phmm_reset(phmm);
  }

  /* before output, have to restore gaps in reference sequence, for
     proper coord conversion */
  if (indels && (post_probs || viterbi)) {
    ss_free(msa->ss); msa->ss = NULL; /* msa->seqs must already exist */
    for (i = 0; i < msa->length; i++) 
      if (msa->seqs[0][i] == msa->missing[0]) msa->seqs[0][i] = GAP_CHAR;
  }
    
  /* Viterbi */
  if (viterbi) {
    GFF_Set *predictions;

    if (!quiet) fprintf(results_f, "Running Viterbi algorithm...\n");
    predictions = phmm_predict_viterbi_cats(phmm, states, seqname, NULL,
                                            idpref, NULL, "phastCons_predicted");
    /* note that selected state numbers are also cat numbers  */
   
    /* score predictions, if necessary */
    if (score) { 
      if (!quiet) fprintf(results_f, "Scoring predictions...\n");            
      phmm_score_predictions(phmm, predictions, states, NULL, NULL, FALSE);
    }

    /* convert GFF to coord frame of reference sequence and adjust
       coords by idx_offset, if necessary  */
    if (refidx != 0 || msa->idx_offset != 0)
      msa_map_gff_coords(msa, predictions, 0, refidx, msa->idx_offset);

    if (refidx != 0) 
      gff_flatten(predictions);	
    /* necessary because coord conversion might create overlapping
       features (can happen in deletions in reference sequence) */

    /* now output predictions */
    if (viterbi_f != NULL) {
      if (gff)
	gff_print_set(viterbi_f, predictions);
      else                        /* BED format */
	gff_print_bed(viterbi_f, predictions, FALSE); 
    }
    if (results != NULL)
      lol_push_gff(results, predictions, two_state ? "most.conserved" :
		   (states != NULL ? "in.states" : "viterbi"));
    gff_free_set(predictions);
  }

  /* posterior probs */
  if (post_probs) {
    int *coord=NULL;

    if (!quiet) fprintf(results_f, "Computing posterior probabilities...\n");

    if (states == NULL) {  //this only happens if two_state==FALSE
                           //return posterior probabilites for every state
      double **postprobs = phmm_new_postprobs(phmm), **postprobsNoMissing=NULL;
      int idx=0, j, k, l;
      if (results != NULL) {
	postprobsNoMissing = smalloc(phmm->hmm->nstates * sizeof(double*));
	for (j=0; j < phmm->hmm->nstates; j++)
	  postprobsNoMissing[j] = smalloc(msa->length*sizeof(double));
	coord = smalloc(msa->length*sizeof(int));
      }
      
      /* print to post_probs_f */
      last = -INFTY;
      for (j = 0, k = 0; j < msa->length; j++) {
	checkInterruptN(j, 1000);
	if (refidx == 0 || msa_get_char(msa, refidx-1, j) != GAP_CHAR) {
	  if (!msa_missing_col(msa, refidx, j)) {
	    if (post_probs_f != NULL) {
	      if (k > last + 1) 
		fprintf(post_probs_f, "fixedStep chrom=%s start=%d step=1\n", seqname, 
			k + msa->idx_offset + 1);
	      for (l=0; l < phmm->hmm->nstates; l++) {
		if (l != 0) fprintf(post_probs_f, "\t");
		fprintf(post_probs_f, "%.3f%c", postprobs[l][j],
			l==phmm->hmm->nstates-1 ? '\n' : '\t');
	      }
	    }
	    if (results != NULL) {
	      coord[idx] = k + msa->idx_offset + 1;
	      for (l=0; l < phmm->hmm->nstates; l++) 
		postprobsNoMissing[l][idx] = postprobs[l][j];
	      idx++;
	    }
	    last = k;
	  }
	  k++;
	}
      }
      if (results != NULL) {
	ListOfLists *wigList = lol_new(2);
	char temp[100];
	lol_push_int(wigList, coord, idx, "coord");
	// fix me: can we get actualy state name from category map?  Problem
        //is some states may have same name.  May need to append strand 
        // and/or index?
	for (j=0; j < phmm->hmm->nstates; j++) {
	  //	  sprintf(temp, "%s", cm_get_feature(cm, state_to_cat(phmm->j)));
	  sprintf(temp, "state.%i", j);
	  lol_push_dbl(wigList, postprobsNoMissing[j], idx, temp);
	}
	lol_set_class(wigList, "data.frame");
	lol_push_lol(results, wigList, "post.prob.wig");
	for (j=0; j < phmm->hmm->nstates; j++) 
	  sfree(postprobsNoMissing[j]);
	sfree(postprobsNoMissing);
	sfree(coord);
      }
      for (j=0; j < phmm->hmm->nstates; j++)
	sfree(postprobs[j]);
      sfree(postprobs);
    } else {
      double *postprobs, *postprobsNoMissing=NULL;
      int idx=0, j, k;
      postprobs = phmm_postprobs_cats(phmm, states, &lnl);
      if (results != NULL) {
	postprobsNoMissing = smalloc(msa->length*sizeof(double));
	coord = smalloc(msa->length*sizeof(int));
      }
      
      /* print to post_probs_f */
      last = -INFTY;
      for (j = 0, k = 0; j < msa->length; j++) {
	checkInterruptN(j, 1000);
	if (refidx == 0 || msa_get_char(msa, refidx-1, j) != GAP_CHAR) {
	  if (!msa_missing_col(msa, refidx, j)) {
	    if (post_probs_f != NULL) {
	      if (k > last + 1) 
		fprintf(post_probs_f, "fixedStep chrom=%s start=%d step=1\n", seqname, 
			k + msa->idx_offset + 1);
	      fprintf(post_probs_f, "%.3f\n", postprobs[j]);
	    }
	    if (results != NULL) {
	      coord[idx] = k + msa->idx_offset + 1;
	      postprobsNoMissing[idx++] = postprobs[j];
	    }
	    last = k;
	  }
	  k++;
	}
      }
      if (results != NULL) {
        ListOfLists *wigList = lol_new(2);
        lol_push_int(wigList, coord, idx, "coord");
        lol_push_dbl(wigList, postprobsNoMissing, idx, "post.prob");
        lol_set_class(wigList, "data.frame");
        lol_push_lol(results, wigList, "post.prob.wig");
        sfree(postprobsNoMissing);
        sfree(coord);
      }
      sfree(postprobs);
    }
  }

  if (compute_likelihood) {
    if (lnl > 0) {              /* may have already been computed */
      if (!quiet) fprintf(results_f, "Computing total log likelihood...\n");
      lnl = phmm_lnl(phmm); 
    }
    if (results != NULL)
      lol_push_dbl(results, &lnl, 1, "likelihood");
    if (lnl_f != NULL) {
      fprintf(lnl_f, "lnL = %.4f\n", lnl); 
      if (FC) fprintf(lnl_f, "(lambda = %g)\n", lambda);
      else if (two_state && (estim_transitions || estim_indels)) {
	fprintf(lnl_f, "(");
	if (estim_transitions)
	  fprintf(lnl_f, "mu = %g, nu = %g%s", mu, nu, estim_indels ? ", " : "");
	if (estim_indels)
	  fprintf(lnl_f, "alpha_0 = %g, beta_0 = %g, tau_0 = %g, alpha_1 = %g, beta_1 = %g, tau_1 = %g", alpha_0, beta_0, tau_0, alpha_1, beta_1, tau_1);
	fprintf(lnl_f, ")\n");
      }
    }
  }

  if (!quiet)
    fprintf(results_f, "Done.\n");

  return 0;
}

/* Set up HMM and category map for two-state case */
void setup_two_state(HMM **hmm, CategoryMap **cm, double mu, double nu) {

  *hmm = hmm_new_nstates(2, TRUE, FALSE);

  /* set HMM transitions according to mu and nu */
  mm_set((*hmm)->transition_matrix, 0, 0, 1-mu);
  mm_set((*hmm)->transition_matrix, 0, 1, mu);
  mm_set((*hmm)->transition_matrix, 1, 0, nu);
  mm_set((*hmm)->transition_matrix, 1, 1, 1-nu);

  /* just use stationary distribution for begin transitions */
  vec_set((*hmm)->begin_transitions, 0, nu/(mu+nu));
  vec_set((*hmm)->begin_transitions, 1, mu/(mu+nu));

  hmm_reset(*hmm);

  /* define two-category category map */
  if (cm != NULL)
    *cm = cm_create_trivial(1, "cons_");
}

/* Version of compute_emissions for use when estimating rho only (see
   fit_two_state, below); makes use of fact that emissions for
   nonconserved state need not be recomputed */
void compute_emissions_estim_rho(double **emissions, void **models, 
				 int nmodels, void *data, int sample, 
				 int length) {
  PhyloHmm *phmm = (PhyloHmm*)data;
  tl_compute_log_likelihood(phmm->mods[0], phmm->em_data->msa, 
			    phmm->emissions[0], NULL,  -1, NULL);
}


/* Estimate parameters for the two-state model using an EM algorithm.
   Any or all of the parameters 'mu' and 'nu', the indel parameters, and
   the tree models themselves may be estimated.  Returns ln
   likelihood. */
double fit_two_state(PhyloHmm *phmm, MSA *msa, int estim_func, int estim_indels,
                     int estim_trees, int estim_rho, double *mu, double *nu, 
                     double *alpha_0, double *beta_0, double *tau_0, 
                     double *alpha_1, double *beta_1, double *tau_1, 
                     double *rho, double gamma, FILE *logf) {
  double retval;
  void (*compute_emissions_func)(double **, void **, int, void*, int, int);

  mm_set(phmm->functional_hmm->transition_matrix, 0, 0, 1-*mu);
  mm_set(phmm->functional_hmm->transition_matrix, 0, 1, *mu);
  mm_set(phmm->functional_hmm->transition_matrix, 1, 0, *nu);
  mm_set(phmm->functional_hmm->transition_matrix, 1, 1, 1-*nu);
                                /* note that phmm->functional_hmm ==
                                   phmm->hmm if no indel model */

  phmm->em_data = smalloc(sizeof(EmData));
  phmm->em_data->msa = msa;
  phmm->em_data->fix_functional = !estim_func;
  phmm->em_data->fix_indel = !estim_indels;
  phmm->em_data->rho = *rho;
  phmm->em_data->gamma = gamma;
  phmm->em_data->H = NULL;      /* will be defined as needed */

  if (phmm->indel_mode == PARAMETERIC) {
    phmm->alpha[0] = *alpha_0;
    phmm->beta[0] = *beta_0;
    phmm->tau[0] = *tau_0;
    phmm->alpha[1] = *alpha_1;
    phmm->beta[1] = *beta_1;
    phmm->tau[1] = *tau_1;
  }

  phmm_reset(phmm); 

  if (estim_trees || estim_rho) {
    msa->ncats = phmm->nmods - 1;   /* ?? */
    if (msa->ss == NULL) 
      ss_from_msas(msa, phmm->mods[0]->order+1, TRUE, NULL, NULL, NULL, -1,
		   subst_mod_is_codon_model(phmm->mods[0]->subst_mod));
    else if (msa->ss->cat_counts == NULL)
      ss_realloc(msa, msa->ss->tuple_size, msa->ss->ntuples, TRUE, TRUE);
  }

  if (estim_trees)
    compute_emissions_func = phmm_compute_emissions_em;
  else if (estim_rho) 
    compute_emissions_func = compute_emissions_estim_rho;
  else compute_emissions_func = NULL;

  if (estim_trees) {
    retval = hmm_train_by_em(phmm->hmm, phmm->mods, phmm, 1, &phmm->alloc_len, NULL, 
                             compute_emissions_func, reestimate_trees,
                             gamma > 0 ? 
                             phmm_estim_trans_em_coverage : phmm_estim_trans_em, 
                             phmm_get_obs_idx_em, 
                             phmm_log_em, phmm->emissions, logf) * log(2);


    /* have to do final rescaling of tree models to get units of subst/site */
    if (phmm->mods[0]->subst_mod != JC69 && phmm->mods[0]->subst_mod != F81) {   
                                /* JC69 and F81 are exceptions */
      tm_scale_model(phmm->mods[0], NULL, 1, 0);
      tm_scale_model(phmm->mods[1], NULL, 1, 0);
    }

    phmm->mods[0]->lnL = phmm->mods[1]->lnL = retval;
  }

  else if (estim_rho) {
    phmm->mods[0]->estimate_branchlens = TM_SCALE_ONLY;
    phmm->mods[0]->scale = phmm->em_data->rho;
    tm_set_subst_matrices(phmm->mods[0]);

    retval = hmm_train_by_em(phmm->hmm, phmm->mods, phmm, 1, &phmm->alloc_len, NULL, 
                             compute_emissions_func, reestimate_rho,
                             gamma > 0 ? 
                             phmm_estim_trans_em_coverage : phmm_estim_trans_em, 
                             phmm_get_obs_idx_em, 
                             phmm_log_em, phmm->emissions, logf) * log(2);

    /* do final rescaling of conserved tree */
    tm_scale_branchlens(phmm->mods[0], phmm->em_data->rho, FALSE);  
    phmm->mods[0]->scale = 1;

    phmm->mods[0]->lnL = phmm->mods[1]->lnL = retval;
  }

  else {                        /* not estimating tree models */
    retval = hmm_train_by_em(phmm->hmm, phmm->mods, phmm, 1, 
                             &phmm->alloc_len, NULL, NULL, NULL,
                             gamma > 0 ? 
                             phmm_estim_trans_em_coverage : phmm_estim_trans_em, 
                             NULL, phmm_log_em, phmm->emissions, logf) * log(2);
  }

  *mu = mm_get(phmm->functional_hmm->transition_matrix, 0, 1);
  *nu = mm_get(phmm->functional_hmm->transition_matrix, 1, 0);
  *rho = phmm->em_data->rho;

  if (phmm->indel_mode == PARAMETERIC) {
    *alpha_0 = phmm->alpha[0];
    *beta_0 = phmm->beta[0];
    *tau_0 = phmm->tau[0];
    *alpha_1 = phmm->alpha[1];
    *beta_1 = phmm->beta[1];
    *tau_1 = phmm->tau[1];
  }

  return retval;
}

/* Special-purpose unpack function, adapted from tm_unpack_params */
void unpack_params_mod(TreeModel *mod, Vector *params_in) {
  TreeNode *n;
  int nodeidx, i;
  List *traversal;
  Vector *params = mod->all_params;

  if (!mod->estimate_ratemat)
    die("ERROR unpack_params_mod: mod->estimate_ratemat is FALSE\n");

  /* check parameter values */
  for (i = 0; i < params_in->size; i++) {
    double mu = vec_get(params_in, i);
    if (mu < 0 && fabs(mu) < TM_IMAG_EPS) /* consider close enough to 0 */
      vec_set(params_in, i, mu=0);
    if (mu < 0) die("ERROR: parameter %d has become negative (%g).\n", i, mu);
    if (isinf(mu) || isnan(mu))
      die("ERROR: parameter %d is no longer finite (%g).\n", i, mu);
  }
  for (i = 0; i<params->size; i++) {
    if (mod->param_map[i] >= 0) 
      vec_set(params, i,
	      vec_get(params_in, mod->param_map[i]));
  }

  if (mod->estimate_branchlens == TM_SCALE_ONLY)
    mod->scale = vec_get(params, mod->scale_idx); 
  else if (mod->estimate_branchlens == TM_BRANCHLENS_ALL) {
    traversal = tr_preorder(mod->tree);
    i=0;
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);
      if (n->parent == NULL) continue;

      if ((n == mod->tree->lchild || n == mod->tree->rchild) && 
	  tm_is_reversible(mod))
	n->dparent = vec_get(params, mod->bl_idx+i)/2.0;
      else 
	n->dparent = vec_get(params, mod->bl_idx+i);
      i++;
      if (n->id == mod->root_leaf_id)
	n->dparent = 0.0;
    }
  }

  /* next parameters are for rate variation */
  if (mod->nratecats > 1) 
    mod->alpha = vec_get(params, mod->ratevar_idx);

  tm_set_rate_matrix(mod, params, mod->ratematrix_idx);
  
/* diagonalize, if necessary */
  if (mod->subst_mod != JC69 && mod->subst_mod != F81)
    mm_diagonalize(mod->rate_matrix);
}


/* Unpack all params for two-state HMM.  Used during M step of EM */
void unpack_params_phmm(PhyloHmm *phmm, Vector *params) {
  unpack_params_mod(phmm->mods[0], params);
  unpack_params_mod(phmm->mods[1], params);
  phmm->em_data->rho = vec_get(params, params->size - 1);
  tm_scale_branchlens(phmm->mods[0], phmm->em_data->rho, FALSE);
  
  if (phmm->mods[0]->nratecats > 1) 
    DiscreteGamma(phmm->mods[0]->freqK, phmm->mods[0]->rK, phmm->mods[0]->alpha, 
                  phmm->mods[0]->alpha, phmm->mods[0]->nratecats, 0); 
  if (phmm->mods[1]->nratecats > 1) 
    DiscreteGamma(phmm->mods[1]->freqK, phmm->mods[1]->rK, phmm->mods[1]->alpha, 
                  phmm->mods[1]->alpha, phmm->mods[1]->nratecats, 0); 
  tm_set_subst_matrices(phmm->mods[0]);
  tm_set_subst_matrices(phmm->mods[1]);
}
 
/* Wrapper for computation of likelihood, for use by reestimate_trees (below) */
double likelihood_wrapper(Vector *params, void *data) {
  PhyloHmm *phmm = (PhyloHmm*)data;
  double retval0, retval1;

  unpack_params_phmm(phmm, params);

  retval0 = -tl_compute_log_likelihood(phmm->mods[0], phmm->em_data->msa, NULL, NULL, 0, NULL);
  retval1 = -tl_compute_log_likelihood(phmm->mods[1], phmm->em_data->msa, NULL, NULL, 1, NULL);

  return retval0 + retval1;
                                /* FIXME: what happens when not one to
                                   one cats and mods? */
}

/* Re-estimate phylogenetic model based on expected counts (M step of EM) */
void reestimate_trees(TreeModel **models, int nmodels, void *data, 
                      double **E, int nobs, FILE *logf) {

  PhyloHmm *phmm = (PhyloHmm*)data;
  int k, obsidx, i, npar;
  Vector *params, *lower_bounds, *upper_bounds, *opt_params;
  double ll;
  int haveratevar, orig_nratecats[2];

  /* FIXME: what about when multiple states per model?  Need to
     collapse sufficient stats.  Could probably be done generally...
     need to use state_to_cat, etc. in deciding which categories to
     use */

  for (k = 0; k < phmm->nmods; k++) 
    for (obsidx = 0; obsidx < nobs; obsidx++) 
      phmm->em_data->msa->ss->cat_counts[k][obsidx] = E[k][obsidx];
  
  /* This will set up params in phmm->mods[0] and phmm->mods[1].  The
     tree models should be the same at this point, since only one model
     is allowed for --estimate-trees.  Therefore the parameter setup
     should be the same. */
  phmm->mods[0]->estimate_ratemat = phmm->mods[1]->estimate_ratemat = TRUE;
  phmm->mods[0]->estimate_backgd = phmm->mods[1]->estimate_backgd = FALSE;
  phmm->mods[0]->estimate_branchlens = phmm->mods[1]->estimate_branchlens =
    TM_BRANCHLENS_ALL;

  /* in order to have the parameter mappings the same in both models,
     we need to call tm_setup_params in both models so that they either
     both have rate variation or both not.
   */
  haveratevar=(phmm->mods[0]->nratecats > 1 ||
	       phmm->mods[1]->nratecats > 1);
  for (i=0; i<2; i++) {
    if (phmm->mods[i]->nratecats > 1 && phmm->mods[i]->empirical_rates) 
      die("ERROR: --estimate-trees not implemented with empirical rate model.");
    orig_nratecats[i] = phmm->mods[i]->nratecats;
    if (haveratevar && phmm->mods[i]->nratecats ==1) 
      phmm->mods[i]->nratecats = phmm->mods[!i]->nratecats;
  }

  tm_setup_params(phmm->mods[0], 0);
  params = tm_params_new_init_from_model(phmm->mods[1]);

  //get number of parameters
  npar = 0 ;
  for (i=0; i < tm_get_nparams(phmm->mods[0]); i++)
    if (phmm->mods[0]->param_map[i] >= npar)
      npar = phmm->mods[1]->param_map[i]+1;
  if (orig_nratecats[1] > 1) {
    if (orig_nratecats[0] > 1) 
      phmm->mods[1]->param_map[phmm->mods[1]->ratevar_idx] = npar++;
    else if (orig_nratecats[0]==1) {
      phmm->mods[0]->param_map[phmm->mods[0]->ratevar_idx] = -1;
      phmm->mods[0]->ratevar_idx = -1;
    }
  }
  npar++;  //make room for rho parameter

  phmm->mods[0]->nratecats = orig_nratecats[0];
  phmm->mods[1]->nratecats = orig_nratecats[1];

  opt_params = vec_new(npar);
  for (i=0; i<params->size; i++)
    if (phmm->mods[0]->param_map[i] >= 0) 
      vec_set(opt_params, phmm->mods[0]->param_map[i],
	      vec_get(params, i));
  for (i=0; i<2; i++) {
    if (phmm->mods[i]->nratecats > 1) 
      vec_set(opt_params, phmm->mods[i]->param_map[phmm->mods[i]->ratevar_idx], 
	      phmm->mods[i]->alpha);
  }
  vec_set(opt_params, npar - 1, phmm->em_data->rho);

  lower_bounds = vec_new(npar);
  vec_zero(lower_bounds);
  upper_bounds = vec_new(npar);
  vec_set_all(upper_bounds, INFTY);
  vec_set(upper_bounds, npar - 1, 1); /* 0 < rho < 1 */

  if (logf != NULL)
    fprintf(logf, "\nRE-ESTIMATION OF TREE MODEL:\n");

  /* keep Hessian arround so it can be used from one iteration to the
     next */
  if (phmm->em_data->H == NULL) {
    phmm->em_data->H = mat_new(npar,npar);
    mat_set_identity(phmm->em_data->H);
  }

  vec_copy(phmm->mods[0]->all_params, params);
  vec_copy(phmm->mods[1]->all_params, params);

  if (opt_bfgs(likelihood_wrapper, opt_params, phmm, &ll, lower_bounds, 
               NULL, logf, NULL, OPT_MED_PREC, phmm->em_data->H, NULL) != 0)
    die("ERROR returned by opt_bfgs.\n");

  if (logf != NULL) 
    fprintf(logf, "END RE-ESTIMATION OF TREE MODEL\n\n");

  unpack_params_phmm(phmm, opt_params);

  if (phmm->indel_mode == PARAMETERIC)
    phmm_set_branch_len_factors(phmm);

  vec_free(params); 
  vec_free(opt_params);
  vec_free(lower_bounds);
  vec_free(upper_bounds);
}


/* Wrapper for computation of likelihood, for use by reestimate_rho (below) */
double likelihood_wrapper_rho(double rho, void *data) {
  PhyloHmm *phmm = (PhyloHmm*)data;
  phmm->mods[0]->scale = rho;
  tm_set_subst_matrices(phmm->mods[0]);
  return -tl_compute_log_likelihood(phmm->mods[0], phmm->em_data->msa, 
				    NULL, NULL, 0, NULL);
}

/* Similar to reestimate_trees, but re-estimate only scale parameter
   rho; only the first model (for the conserved state) needs to be
   considered */
void reestimate_rho(TreeModel **models, int nmodels, void *data, 
		    double **E, int nobs, FILE *logf) {

  PhyloHmm *phmm = (PhyloHmm*)data;
  int obsidx;
  double ax, bx, cx, fa, fb, fc;

  for (obsidx = 0; obsidx < nobs; obsidx++) 
    phmm->em_data->msa->ss->cat_counts[0][obsidx] = E[0][obsidx];

  if (logf != NULL)
    fprintf(logf, "\nRE-ESTIMATION OF RHO (BRENT'S METHOD):\n");

  bx = phmm->em_data->rho;
  ax = max(0.1, phmm->em_data->rho - .05);
  mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, likelihood_wrapper_rho, phmm, logf);
  opt_brent(ax, bx, cx, likelihood_wrapper_rho, 5e-3, 
	    &phmm->em_data->rho, phmm, logf);
  //  printf("ll=%f rho=%f\n", ll, phmm->em_data->rho);

  if (logf != NULL) 
    fprintf(logf, "END RE-ESTIMATION OF RHO\n\n");

  if (phmm->indel_mode == PARAMETERIC)
    die("ERROR reestimate:rho: phmm->indel_mode is PARAMETERIC\n");
  /* FIXME: to make work with parameteric indel model, will have to
     propagate scale parameter through phmm_set_branch_len_factors */
}

/* Maximize HMM transition parameters subject to constrain implied by
   target coverage (M step of EM).  For use with two-state HMM.  This
   function is passed to hmm_train_by_em in phmm_fit_em */
void phmm_estim_trans_em_coverage(HMM *hmm, void *data, double **A) {

  PhyloHmm *phmm = data;
  IndelEstimData *ied = NULL;
  double **C;

  if (phmm->em_data->fix_functional && phmm->em_data->fix_indel) return;

  if (phmm->indel_mode == PARAMETERIC) {
    ied = phmm_new_ied(phmm, A);
    C = ied->fcounts;
  }
  else C = A;

  /* estimate transition probs for functional cats subject to
     constraint on coverage */
  if (!phmm->em_data->fix_functional) {
    double a, b, c, mu, nu, nu1, nu2, z, tmp;
    /* if you take the first derivative wrt nu of the expression inside
       the argmax and set it to zero, you get a quadratic eqn which
       can be solved using the quadratic formula */
    z = (1-phmm->em_data->gamma)/phmm->em_data->gamma;
    a = z * (C[0][0] + C[0][1] + C[1][0] + C[1][1]);
    b = -C[0][1] - C[1][0] - C[1][1] - z * (C[0][0] + C[0][1] + C[1][0]);
    c = C[0][1] + C[1][0];

    tmp = b*b - 4*a*c;
    if (tmp < 0)
      die("ERROR phmm_estim_trans_em_coverage: tmp=%e\n", tmp);
    tmp = sqrt(tmp);
    nu1 = (-b + tmp) / (2*a);
    nu2 = (-b - tmp) / (2*a);
    /* only one root can be valid */
    if (nu1 < 1e-10 || z * nu1 > 1 - 1e-10)                                 
      nu = nu2;                   /* (allow for rounding errors) */
    else nu = nu1;

    /* double check that derivative is really zero */
//    if (!(fabs(-z*C[0][0]/(1-z*nu) + (C[0][1] + C[1][0])/nu - C[1][1]/(1-nu)) < 1e-4))
//     die("ERROR phmm_estim_trans_em_coverage: derivative not zero?\n");

    mu = z * nu;
    if (!(nu >= 0 && nu <= 1 && mu >= 0 && mu <= 1))
      die("ERROR phmm_estim_trans_em_coverage: mu=%e, nu=%e\n", mu, nu);

    mm_set(phmm->functional_hmm->transition_matrix, 0, 0, 1-mu);
    mm_set(phmm->functional_hmm->transition_matrix, 0, 1, mu);
    mm_set(phmm->functional_hmm->transition_matrix, 1, 0, nu);
    mm_set(phmm->functional_hmm->transition_matrix, 1, 1, 1-nu);

    /* use stationary distribution for begin transitions */
    vec_set(phmm->functional_hmm->begin_transitions, 0, nu/(mu+nu));
    vec_set(phmm->functional_hmm->begin_transitions, 1, mu/(mu+nu));
  }

  if (phmm->indel_mode == PARAMETERIC) {
    if (!phmm->em_data->fix_indel) phmm_em_estim_indels(phmm, ied);
    phmm_free_ied(ied);
  }

  phmm_reset(phmm);
}

/* initialize equilibrium freqs for tree model; either make consistent
   with given G+C content or estimate from alignment */
void init_eqfreqs(TreeModel *mod, MSA *msa, double gc) {
  if (gc != -1) {               /* gc specified */
    if (strlen(mod->rate_matrix->states) != 4 || 
        mod->rate_matrix->inv_states[(int)'A'] < 0 ||
        mod->rate_matrix->inv_states[(int)'C'] < 0 ||
        mod->rate_matrix->inv_states[(int)'G'] < 0 ||
        mod->rate_matrix->inv_states[(int)'T'] < 0)
      die("ERROR: Four-character DNA alphabet required with --gc.\n");
    if (! (gc > 0 && gc < 1))
      die("ERROR init_eqfreqs got gc=%e\n", gc);
    vec_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'G'], gc/2);
    vec_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'C'], gc/2);
    vec_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'A'], (1-gc)/2);
    vec_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'T'], (1-gc)/2);
  }
  else {                        /* estimate from alignment */
    if (subst_mod_is_codon_model(mod->subst_mod))
      msa_get_backgd_3x4(mod->backgd_freqs, msa);
    else if (mod->subst_mod == JC69 || mod->subst_mod == K80)
      vec_set_all(mod->backgd_freqs, 1.0/mod->backgd_freqs->size);
    else
      msa_get_base_freqs_tuples(msa, mod->backgd_freqs, mod->order+1, -1);
  }
}
