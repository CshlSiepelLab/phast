/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include "bgc_hmm.h"

/* 
   Like phastCons, but with two versions of each state: with and without
   gBGC acting on a particular subtree.

   Parameters that describe this HMM:
   tree scale- hold at 1? (ie, assume input is neutral model)
   conserved tree scale (hold at 0.3?)
   transition rate cons->neutral (use phastCons estimates?)
   transition rate neutral->cons (use phastCons estimates?)
   rate bgc->no_bgc
   rate no_bgc->bgc
   bgc  parameter

   We can use EM to estimate rates between bgc and no bgc states, and
   estimate the rest of the parameters with tm_fit_multi (which can
   optimize parameters shared across multiple models).
 */


/* set up the tree models representing the 4 HMM states, initialize
   parameters so that they are appropriately shared between the states
   and determine which ones to optimize.
 */
TreeModel **bgchmm_setup_mods(TreeModel *init_mod, 
			      char *foregd_branch,
			      int do_bgc,
			      double bgc,
			      double rho,
			      double init_scale,
			      int estimate_bgc,
			      int estimate_rho,
			      int estimate_scale,
			      int eqfreqs_from_msa,
			      MSA *align,
			      int *npar_rv) {
  int npar=0, nmod, i;
  TreeModel **mods;
  subst_mod_type subst_mod;
  Vector *subst_mod_params;
  String *altmod_str=NULL;
  char tempstr[1000];
  AltSubstMod *altmod;
  int scale_pos=-1, cons_scale_pos=-1, bgc_pos=-1;
  double scale;

  //set npar to the number of parameters that we want to optimize.  
  // All others held constant.

  npar=0;
  if (estimate_scale) {
    fprintf(stderr, "estimating tree scale\n");
    scale_pos = npar++;
  } else scale_pos = -1;
  if (estimate_rho) {
    fprintf(stderr, "estimating rho (conserved state scale)\n");
    cons_scale_pos = npar++;
  } else cons_scale_pos = -1;
  if (estimate_bgc && do_bgc) {
    fprintf(stderr, "estimating bgc strength parameter B\n");
    bgc_pos = npar++;
  } else bgc_pos = -1;
  if (npar_rv != NULL) 
    *npar_rv = npar;

  if (do_bgc) {
    sprintf(tempstr, "%s:%s", foregd_branch, "bgc[0,]");
    altmod_str = str_new_charstr(tempstr);
    nmod=4;
  } else nmod=2;

  mods = smalloc(nmod * sizeof(TreeModel*));

  nmod = 0;
  mods[nmod] = tm_create_copy(init_mod);
  if (mods[nmod]->alt_subst_mods != NULL)
    tm_free_alt_subst_mods(mods[nmod]);
  subst_mod = mods[nmod]->subst_mod;
  subst_mod_params = vec_new(tm_get_nratematparams(mods[nmod]));
  if (mods[nmod]->selection != 0.0) {
    tm_unapply_selection_bgc(mods[nmod]->rate_matrix, mods[nmod]->selection, 0.0);
    mods[nmod]->selection = 0.0;
  }
  if (eqfreqs_from_msa) 
    init_eqfreqs(mods[nmod], align, -1);
  tm_rate_params_init_from_model(mods[nmod], subst_mod_params, 0, 0.0, 0.0);
  tm_set_rate_matrix(mods[nmod], subst_mod_params, 0);
  mods[nmod]->scale_during_opt = 1;
  if (mods[nmod]->nratecats > 1) {
    tm_reinit(mods[nmod], subst_mod, 1, 0.0, NULL, NULL);
    tm_set_rate_matrix(mods[nmod], subst_mod_params, 0);
  }
  scale = tm_scale_rate_matrix(mods[nmod]) * init_scale;
  tr_scale(mods[nmod]->tree, scale);

  //mods[0] is neutral, no bgc, noncoding
  mods[nmod]->estimate_branchlens = TM_SCALE_ONLY;  //note scale parameter will not be optimized unless estimate_scale=TRUE
  mods[nmod]->estimate_backgd = FALSE;
  mods[nmod]->estimate_ratemat = FALSE;
  tm_setup_params(mods[nmod], 0);
  mods[nmod]->param_map[mods[nmod]->scale_idx] = scale_pos;
  tm_params_init_from_model(mods[nmod], mods[nmod]->all_params, FALSE);
  if (estimate_scale) {
    mods[nmod]->bound_arg = lst_new_ptr(1);
    lst_push_ptr(mods[nmod]->bound_arg, str_new_charstr("scale[0,200]"));
  }
  nmod++;

  //mods[1] is conserved, no bgc.  Share all params with mod[0] except scale_sub
  mods[nmod] = tm_create_copy(mods[0]);
  mods[nmod]->scale_sub = rho;
  vec_set(mods[nmod]->all_params, mods[nmod]->scale_idx+1, rho);
  mods[nmod]->in_subtree = smalloc(mods[nmod]->tree->nnodes * sizeof(int));
  for (i=0; i < mods[nmod]->tree->nnodes; i++)
    mods[nmod]->in_subtree[i] = 1;
  mods[nmod]->param_map[mods[nmod]->scale_idx+1] = cons_scale_pos;
  nmod++;

  if (do_bgc) {
    //next mod is neutral, bgc.  Share all params with mod[0] except bgc
    mods[nmod] = tm_create_copy(mods[0]);
    tm_add_alt_mod(mods[nmod], altmod_str);
    altmod = lst_get_ptr(mods[nmod]->alt_subst_mods, 0);
    altmod->bgc = bgc;
    tm_apply_selection_bgc(altmod->rate_matrix, 0.0, bgc);
    tm_setup_params(mods[nmod], 0);
    tm_params_init_from_model(mods[nmod], mods[nmod]->all_params, FALSE);
    mods[nmod]->param_map[mods[nmod]->scale_idx] = scale_pos;
    mods[nmod]->param_map[altmod->bgc_idx] = bgc_pos;
    nmod++;
    
    //mods[3] is conserved, bgc.  Share scale with mods[1] and bgc with mods[2]
    mods[nmod] = tm_create_copy(mods[2]);
    mods[nmod]->scale_sub = rho;
    mods[nmod]->in_subtree = smalloc(mods[nmod]->tree->nnodes * sizeof(int));
    for (i=0; i < mods[nmod]->tree->nnodes; i++)
      mods[nmod]->in_subtree[i] = 1;
    vec_set(mods[nmod]->all_params, mods[nmod]->scale_idx + 1, rho);
    nmod++;
  }

  if (altmod_str != NULL) str_free(altmod_str);
  vec_free(subst_mod_params);

  for (i=0; i < nmod; i++) {
    if (mods[i] != NULL) {
      vec_set(mods[i]->all_params, mods[i]->scale_idx, 1.0);
      mods[i]->scale = 1.0;
    }
  }
  return mods;
}


struct bgchmm_struct *bgchmm_struct_new(int rphast) {
  struct bgchmm_struct *rv = smalloc(sizeof(struct bgchmm_struct));
  rv->msa = NULL;  //msa and mod need to be set before bgcHmm can be called
  rv->mod = NULL;
  rv->scale = 1.0;
  rv->rho = 0.31;
  rv->cons_expected_length = 45;
  rv->cons_target_coverage = 0.3;
  rv->bgc_target_coverage = 0.01;
  rv->bgc_expected_length = 1000;
  rv->bgc = 3.0;
  rv->foregd_branch = NULL;  //this needs to be set too unless gff==NULL and do_bgc==FALSE
  rv->estimate_bgc_target_coverage = TRUE;
  rv->estimate_bgc_expected_length = FALSE;
  rv->estimate_scale = FALSE;
  rv->do_bgc = TRUE;
  rv->estimate_bgc = FALSE;
  rv->estimate_cons_transitions = FALSE;
  rv->estimate_rho = FALSE;
  rv->eqfreqs_from_msa = TRUE;
  rv->results = rphast ? lol_new(15) : NULL;
  rv->viterbi_fn = NULL;  //need to set if want viterbi path output to file
  rv->post_probs = WIG;
  rv->post_probs_f = rphast ? NULL : stdout;
  rv->random_path = 0;
  rv->get_likelihoods = FALSE;
  rv->informative_only = FALSE;
  rv->non_informative_fn = NULL;
  rv->mods_fn = NULL;
  rv->tract_fn = NULL;
  return rv;
}


int bgcHmm(struct bgchmm_struct *b) {
  MSA *msa=b->msa;
  TreeModel **mods;
  int i, j, numstate, nsite, *path, npar;
  double mu, nu,likelihood,
    **emissions, path_likelihood,
    bgc_in_rate, bgc_out_rate;
  HMM *hmm;
  struct bgchmm_data_struct *data;
  int do_bgc;
  ListOfLists *results = b->results;
  FILE *post_probs_f = b->post_probs_f;
  void bgchmm_print_non_informative(struct bgchmm_data_struct *data, ListOfLists *results,
				    char *non_informative_fn);

  /*  TODO: prune tree so it has only sequences in msa; prune MSA so it only has sequences in tree, make warning. 
   */
  List *pruned_names = lst_new_ptr(msa->nseqs);
  int old_nnodes = b->mod->tree->nnodes;

  tm_prune(b->mod, msa, pruned_names);
  if (lst_size(pruned_names) != 0) {
    char warning_str[10000];
    if (lst_size(pruned_names) == (old_nnodes + 1)/2) 
      die("ERROR: no match for leaves of tree in alignment (leaf names must match alignment names");
    sprintf(warning_str, "WARNING: Pruned away leaves of tree with no match in alignment (%s", 
	    ((String*)lst_get_ptr(pruned_names, 0))->chars);
    for (j=1; j < lst_size(pruned_names); j++) 
      sprintf(warning_str, "%s, %s", warning_str, ((String*)lst_get_ptr(pruned_names, j))->chars);
    sprintf(warning_str, "%s)\n", warning_str);
    phast_warning(warning_str);
  }


  do_bgc = b->do_bgc || b->informative_only;
  //make sure we have ordered sufficient stats with tuple_size 1, or get them
  if (msa->ss == NULL || msa->ss->tuple_size != 1 || msa->ss->tuple_idx == NULL) {
    if (msa->seqs == NULL) {
      if (msa->ss->tuple_idx == NULL)
	die("Error: Need ordered alignment to run phylo-HMM");
      ss_to_msa(msa);
    }
    if (msa->ss != NULL)
      ss_free(msa->ss);
    ss_from_msas(msa, 1, 1, NULL, NULL, NULL, -1, 0);
  }

  mu = 1.0/b->cons_expected_length;
  nu = b->cons_target_coverage/(1.0-b->cons_target_coverage) * mu;

  bgc_out_rate = 1.0/b->bgc_expected_length;
  bgc_in_rate = b->bgc_target_coverage/(1.0-b->bgc_target_coverage) * bgc_out_rate;

  data = smalloc(sizeof(struct bgchmm_data_struct));
  data->nsite = msa->length;
  data->msa = msa;
  data->bgc_informative = NULL;
  data->estimate_cons_transitions = b->estimate_cons_transitions;
  data->bgc_target_coverage = b->bgc_target_coverage;
  data->bgc_expected_length = b->bgc_expected_length;
  data->estimate_bgc_target_coverage = b->estimate_bgc_target_coverage;
  data->estimate_bgc_expected_length = b->estimate_bgc_expected_length;
  if (data->estimate_bgc_target_coverage) 
    fprintf(stderr, "estimating gBGC target coverage by EM\n");
  if (data->estimate_bgc_expected_length)
    fprintf(stderr, "estimating gBGC expected length by EM\n");
  if (data->estimate_cons_transitions)
    fprintf(stderr, "estimating transition rates between neutral/conserved states by EM\n");

  if (do_bgc)
    numstate = 4;
  else numstate = 2;

  //create hmm and initialize transition rates
  hmm = hmm_new_nstates(numstate, TRUE, FALSE);
  bgchmm_set_hmm(hmm, bgc_in_rate, bgc_out_rate, nu, mu);
  data->hmm = hmm;

  if (do_bgc) {
    data->bgc_informative = bgchmm_get_informative(data->msa, b->foregd_branch, b->mod->tree);
    if (b->informative_only || results != NULL || b->non_informative_fn != NULL)
      bgchmm_print_non_informative(data, results, b->non_informative_fn);
    if (b->informative_only) return 0;
  } else data->bgc_informative=NULL;

  //set up treeModels and setup param_map so that relevant parameters are optimized
  mods = bgchmm_setup_mods(b->mod, b->foregd_branch, do_bgc, 
			   b->bgc, b->rho, b->scale, 
			   b->estimate_bgc, b->estimate_rho, b->estimate_scale,
			   b->eqfreqs_from_msa, msa, &npar);
  
  nsite = (int)msa->length;
  emissions = smalloc(hmm->nstates * sizeof(double*));
  for (i=0; i < hmm->nstates; i++)
    emissions[i] = smalloc(nsite * sizeof(double));

  //Here is the main function call: optimize the HMM!
  fprintf(stderr, "calling hmm_train_by_em...");
  likelihood=hmm_train_by_em(hmm, mods, data, 1, &nsite, NULL, 
			     bgchmm_compute_emissions,
			     npar > 0 ? bgchmm_estimate_states : NULL, 
			     bgchmm_estimate_transitions, 
			     bgchmm_get_obs_idx, 
			     NULL, emissions, NULL);
  fprintf(stderr, "Done.\n\n");
  bgchmm_get_rates(hmm, &bgc_in_rate, &bgc_out_rate, &nu, &mu);
  
  //output results
  fprintf(stderr, "likelihood=%f\n", likelihood);
  if (do_bgc) {
    fprintf(stderr, "bgc_in_rate=%g\n", bgc_in_rate);
    fprintf(stderr, "bgc_out_rate=%g\n", bgc_out_rate);
  }
  fprintf(stderr, "mu=%g\n", mu);
  fprintf(stderr, "nu=%g\n", nu);
  fprintf(stderr, "scale=%g\n", mods[0]->scale);
  fprintf(stderr, "rho=%g\n", mods[1]->scale_sub);
  if (do_bgc) fprintf(stderr, "bgc=%g\n", ((AltSubstMod*)lst_get_ptr(mods[2]->alt_subst_mods, 0))->bgc);
  
  if (results != NULL) {
    lol_push_dbl(results, &likelihood, 1, "likelihood");
    if (do_bgc) {
      lol_push_dbl(results, &bgc_in_rate, 1, "bgc.in");
      lol_push_dbl(results, &bgc_out_rate, 1, "bgc.out");
    }
    lol_push_dbl(results, &mu, 1, "mu");
    lol_push_dbl(results, &nu, 1, "nu");
    lol_push_dbl(results, &mods[0]->scale, 1, "scale");
    lol_push_dbl(results, &mods[1]->scale_sub, 1, "rho");
    if (do_bgc) 
      lol_push_dbl(results, &(((AltSubstMod*)lst_get_ptr(mods[2]->alt_subst_mods, 0))->bgc), 1, "bgc");
  }

  if (b->mods_fn != NULL) {
    FILE *outfile;
    char outfn[1000];
    for (i = 0 ; i < hmm->nstates; i++) {
      sprintf(outfn, "%s.%s.mod", b->mods_fn, bgchmm_get_state_name(i, do_bgc));
      outfile = phast_fopen(outfn, "w");
      tm_print(outfile, mods[i]);
      fclose(outfile);
    }
  }

  //output viterbi path - don't bother with this anymore since viterbi results don't correspond to posterior probs very well - threshold posteriors at 0.5 instead
  /*  if (b->viterbi_fn != NULL || results != NULL) {  //need to make a GFF_Set out of viterbi elements and push them onto results and/or output them to viterbi_fn
    ListOfLists *viterbi_lol = lol_new(2);
    
    path = smalloc(nsite * sizeof(int));
    hmm_viterbi(hmm, emissions, nsite, path);
    if (b->get_likelihoods) {
      path_likelihood = hmm_path_likelihood(hmm, emissions, nsite, path);
      lol_push_dbl(viterbi_lol, &path_likelihood, 1, "likelihood");
    }
    bgchmm_output_path(path, nsite, msa, do_bgc, "features", 
		       b->viterbi_fn, viterbi_lol);
    if (results != NULL)
      lol_push_lol(results, viterbi_lol, "viterbi");
    sfree(path);
    }*/
  
  //output random paths
  if (b->random_path > 0 && results != NULL) {
    ListOfLists *random_paths = lol_new(b->random_path);
    ListOfLists *currpath;
    double **forward_scores = smalloc(hmm->nstates * sizeof(double*));
    char tempname[1000];
    path = smalloc(nsite * sizeof(int));
    for (i=0; i< hmm->nstates; i++)
      forward_scores[i] = smalloc(nsite * sizeof(double));
    hmm_forward(hmm, emissions, nsite, forward_scores);
    for (i=0; i < b->random_path; i++) {
      currpath = lol_new(2);
      sprintf(tempname, "random.path.%i", i+1);
      hmm_stochastic_traceback(hmm, forward_scores, nsite, path);
      if (b->get_likelihoods) {
	path_likelihood = hmm_path_likelihood(hmm, emissions, nsite, path);
	lol_push_dbl(currpath, &path_likelihood, 1, "likelihood");
      }
      bgchmm_output_path(path, nsite, msa, do_bgc, "features", NULL, currpath);
      lol_push_lol(random_paths, currpath, tempname);
    }
    lol_push_lol(results, random_paths, "random.path");
    for (i=0; i < hmm->nstates; i++) 
      sfree(forward_scores[i]);
    sfree(forward_scores);
    sfree(path);
  }
  
  //now look at posteriors
  if (b->post_probs != NONE || b->tract_fn) {
    double **postprobs, *prob_bgc;
    int k, *coord=NULL, reflen;
    
    postprobs = smalloc(hmm->nstates * sizeof(double*));
    for (i=0; i < hmm->nstates; i++)
      postprobs[i] = smalloc(msa->length * sizeof(double));
    
    hmm_posterior_probs(hmm, emissions, msa->length, postprobs);
    
    /*  For now, we always want to do things in the frame of reference of the
	first species.  So condense the probs.
    */
    for (i=0,j=0; i < msa->length; i++) {
      if (msa_get_char(msa, 0, i) != GAP_CHAR) {
	if (i != j) {
	  for (k=0; k < hmm->nstates; k++) 
	    postprobs[k][j] = postprobs[k][i];
	}
	j++;
      }
    }
    reflen = j;
    prob_bgc = smalloc(reflen * sizeof(double));
    for (i=0; i < reflen; i++)
      prob_bgc[i] = do_bgc ? postprobs[2][i] + postprobs[3][i] : postprobs[1][i];
    
    if (results != NULL) {
      coord = smalloc(reflen*sizeof(int));
      for (i=0; i < reflen; i++) coord[i] = msa->idx_offset + i + 1;
    }
    
    /* print to post_probs_f */
    if (post_probs_f != NULL) {
      if (b->post_probs == WIG) {
	fprintf(post_probs_f, "fixedStep chrom=%s start=%d step=1\n",
		msa->names[0], msa->idx_offset + 1);
	for (i=0; i < reflen; i++)
	  fprintf(post_probs_f, "%.5g\n", prob_bgc[i]);
      } else if (b->post_probs == FULL) {
	if (post_probs_f != NULL) {
	  fprintf(post_probs_f, "#coord");
	  for (j=0; j < hmm->nstates; j++)
	    fprintf(post_probs_f, "\t%s", bgchmm_get_state_name(j, do_bgc));
	  fprintf(post_probs_f, "\n");
	  for (j=0; j < reflen; j++) {
	    fprintf(post_probs_f, "%i", msa->idx_offset + 1 + j);
	    for (k=0; k < hmm->nstates; k++)
	      fprintf(post_probs_f, "\t%0.5g", postprobs[k][j]);
	    fprintf(post_probs_f, "\n");
	  }
	}
      }
    }
    
    if (results != NULL) {
      ListOfLists *wigList = lol_new(2);
      lol_push_int(wigList, coord, reflen, "coord");
      // fix me: can we get actualy state name from category map?  Problem
      //is some states may have same name.  May need to append strand 
      // and/or index?
      if (b->post_probs == WIG)
	lol_push_dbl(wigList, prob_bgc, reflen, do_bgc ? "prob_gBGC" : "prob_cons");
      else {
	for (j=0; j < hmm->nstates; j++) 
	  lol_push_dbl(wigList, postprobs[j], reflen, bgchmm_get_state_name(j, do_bgc));
      }
      lol_set_class(wigList, "data.frame");
      lol_push_lol(results, wigList, "post.prob");
      sfree(coord);
    }
    
    if (b->tract_fn != NULL || results != NULL) {
      FILE *outfile;
      GFF_Set *bgc_tracts = gff_from_wig_threshold(msa->names[0], 
						   msa->idx_offset+1,
						   prob_bgc, reflen, 0.5, "phastBias",
						   "gBGC_tract");
      if (b->tract_fn != NULL) {
	outfile = phast_fopen(b->tract_fn, "w");
	gff_print_set(outfile, bgc_tracts);
	fclose(outfile);
      }
      if (results != NULL)
	lol_push_gff_ptr(results, bgc_tracts, "tracts");
      else gff_free_set(bgc_tracts);
    }
    for (j=0; j < hmm->nstates; j++)
      sfree(postprobs[j]);
    sfree(postprobs);
    sfree(prob_bgc);
  }
  return 0;
}

void bgchmm_print_non_informative(struct bgchmm_data_struct *data, ListOfLists *results,
				  char *non_informative_fn) {
  //create a GFF containing non-informative regions of alignment (if any)
  GFF_Set *non_informative_gff;
  GFF_Feature *feat;
  int *non_informative, i, j;
  MSA *msa = data->msa;

  non_informative = smalloc(msa->length * sizeof(int));
  for (i=0; i < msa->length; i++) 
    non_informative[i] = !data->bgc_informative[msa->ss->tuple_idx[i]];
  
  non_informative_gff = gff_new_set();
  for (i=0; i < msa->length; i++) {
    if (msa_get_char(msa, 0, i) == GAP_CHAR) continue;
    if (non_informative[i]) {
      for (j=i+1; j < msa->length; j++)
	if (!non_informative[j]) break;
      feat = gff_new_feature_copy_chars(msa->names[0], "bgcHmm", "non_informative",
					i+1, j,
					0, '+', GFF_NULL_FRAME, ".", TRUE);
      lst_push_ptr(non_informative_gff->features, feat);
      i = j-1;
    }
  }
  msa_map_gff_coords(msa, non_informative_gff, 0, 1, msa->idx_offset);
  if (results != NULL)
    lol_push_gff(results, non_informative_gff, "not.informative");
  if (non_informative_fn != NULL) {
    FILE *outfile = phast_fopen(non_informative_fn, "w");
    gff_print_set(outfile, non_informative_gff);
    fclose(outfile);
  }
}

char *bgchmm_get_state_name(int state, int do_bgc) {
  switch (state) {
  case 0: 
    return("neutral");
  case 1: 
    return("conserved");
  case 2: 
    if (!do_bgc) die("Got state=%i but do_bgc==FALSE\n", state);
    return("gBGC_neutral");
  case 3: 
    if (!do_bgc) die("Got state=%i but do_bgc==FALSE\n", state);
    return("gBGC_conserved");
  }
  die("ERROR bgchmm_get_state_name got state=%i do_bgc=%i\n", state, do_bgc);
  return("error");
}


void bgchmm_output_path(int *path, int nsite, MSA *msa, int do_bgc,
			char *name, char *outfn, ListOfLists *results) {
  GFF_Set *gff = gff_new_set();
  GFF_Feature *feat;
  int start=-1, state=-1, coord = msa->idx_offset, i;
  FILE *outfile;
  char *feat_type=NULL;
  for (i=0; i < nsite; i++) {
    if (msa_get_char(msa, 0, i) != GAP_CHAR) coord++;  //these are one-based coords
    if (path[i] != state) {
      if (state != -1) {
	feat = gff_new_feature_copy_chars(msa->names[0],
					  "bgcHmm", feat_type, 
					  start, coord-1, 0, '+', 
					  GFF_NULL_FRAME, ".", TRUE);
	lst_push_ptr(gff->features, feat);
      }
      state = path[i];
      start = coord;
      feat_type = bgchmm_get_state_name(state, do_bgc);
    }
  }
  feat = gff_new_feature_copy_chars(msa->names[0],
				    "bgcHmm", feat_type, 
				    start, coord, 0, '+', 
				    GFF_NULL_FRAME, ".", TRUE);
  lst_push_ptr(gff->features, feat);
  if (results != NULL)
    lol_push_gff(results, gff, name);
  if (outfn != NULL) {
    outfile = phast_fopen(outfn, "w");
    gff_print_set(outfile, gff);
    fclose(outfile);
  }
  gff_free_set(gff);
}


//return an array of coordinates containing all tuples in MSA that are not 
// informative enough to detect bgc on a particular branch.  We require at least 
// one piece of non-missing data above the sister node and an additional outgroup 
// node.  If the node is a leaf, there must be data at this leaf.  Otherwise, 
// there must be data above each of the children of the node.
//FIXME: this is a bit over-zealous if bgc acts on several branches;
// currently it combines the requirements for all branches.
int *bgchmm_get_informative(MSA *msa, char *foregd, TreeNode *tree) {
  List *foregd_nodes = lst_new_ptr(1), *leaf_names=NULL, *informative_nodes = lst_new_ptr(4), *inside=lst_new_ptr(tree->nnodes/2), *outside = lst_new_ptr(tree->nnodes/2);
  TreeNode *node, *parent, *sib;
  int i, j, k, tupleidx, col, *spec, *rv, numleaf, have_informative;
  char c;
  //use msa_get_informative_feats(msa, 1, List *specList, 0, 0)

  //first need to find all foregd branches
  node = tr_get_node(tree, foregd);
  if (node != NULL) lst_push_ptr(foregd_nodes, node);
  else tr_get_labelled_nodes(tree, foregd, foregd_nodes);
  if (lst_size(foregd_nodes)==0)
    die("No nodes with name or label %s\n", foregd);
  rv = smalloc(msa->ss->ntuples * sizeof(int));
  for (i=0; i < msa->ss->ntuples; i++) rv[i] = 1;

  for (i=0; i < lst_size(foregd_nodes); i++) {
    node = lst_get_ptr(foregd_nodes, i);
    parent = node->parent;
    if (parent ==  NULL) {
      if (lst_size(foregd_nodes) == 1) die("root node cannot be bgc branch");
      continue;
    }
    if (parent->lchild == node) sib = parent->rchild;
    else sib = parent->lchild;
    
    if (node->lchild == NULL)
      lst_push_ptr(informative_nodes, node);
    else {
      lst_push_ptr(informative_nodes, node->lchild);
      lst_push_ptr(informative_nodes, node->rchild);
    }
    lst_push_ptr(informative_nodes, sib);

    tr_partition_leaves(tree, parent, inside, outside);
    if (lst_size(outside) != 0)
      leaf_names = lst_new_ptr(lst_size(outside));
    for (j=0; j < lst_size(outside); j++) 
      lst_push_ptr(leaf_names, str_new_charstr(((TreeNode*)lst_get_ptr(outside, j))->name));

    //first check leaf_names defined above, then the descendants of each of
    //the nodes in informative_nodes.  Each check must have at least one
    //leaf with information, otherwise this column is not informative.
    for (j=-1; j < lst_size(informative_nodes); j++) {
      if (j != -1)  //first loop leaf_names already defined from outside
	leaf_names = tr_leaf_names((TreeNode*)lst_get_ptr(informative_nodes, j));
      if (leaf_names == NULL || lst_size(leaf_names) == 0) continue;
      numleaf = lst_size(leaf_names);
      spec = smalloc(numleaf*sizeof(int));
      for (k=0; k < numleaf; k++)
	spec[k] = msa_get_seq_idx(msa, ((String*)lst_get_ptr(leaf_names, k))->chars);
      for (tupleidx = 0; tupleidx < msa->ss->ntuples; tupleidx++) {
	have_informative = 0;  //set to 1 when we find a species from leaf_names that is informative in all cols
	for (k=0; k < lst_size(leaf_names); k++) {
	  if (spec[k] == -1) continue;
	  for (col= 1-msa->ss->tuple_size; col <= 0; col++)  {
	    c = ss_get_char_tuple(msa, tupleidx, spec[k], col);
	    if (c == GAP_CHAR || msa->is_missing[(int)c]) break;
	  }
	  if (col==1) {
	    have_informative = 1;
	    break;
	  }
	}
	if (!have_informative) {
	  //one of the 3 or four nodes is not informative; set to 0 and go on
	  rv[tupleidx] = 0;
	  continue;
	}
      }
      sfree(spec);
      lst_free_strings(leaf_names);
      lst_free(leaf_names);
    }
    lst_clear(informative_nodes);
  }
  lst_free(informative_nodes);
  lst_free(foregd_nodes);
  lst_free(inside);
  lst_free(outside);

  return rv;
}


void bgchmm_estimate_states(TreeModel **mods, int nmod, void *data0, double **E, 
			    int nobs, FILE *logfile) {
  struct bgchmm_data_struct *data = (struct bgchmm_data_struct*)data0;
  MSA *msa;
  int use_nmod, i, j;
  TreeModel **usemods;

  msa = data->msa;
  if (msa->ncats != nmod && msa->ss->cat_counts != NULL) {
    for (i=0; i < msa->ncats; i++) sfree(msa->ss->cat_counts[i]);
    sfree(msa->ss->cat_counts);
    msa->ss->cat_counts = NULL;
  }
  if (msa->ss->cat_counts == NULL) {
    msa->ss->cat_counts = smalloc(nmod * sizeof(double*));
    for (i=0; i < nmod; i++)
      msa->ss->cat_counts[i] = smalloc(msa->ss->ntuples * sizeof(double));
  }
  for (i=0; i < nmod; i++)
    for (j=0; j < msa->ss->ntuples; j++)
      msa->ss->cat_counts[i][j] = E[i][j];
  if (data->bgc_informative != NULL && nmod==4) {
    for (i=0; i < msa->ss->ntuples; i++) {
      if (data->bgc_informative[i]==0) {
	for (j=2; j<=3; j++) {
	  msa->ss->cat_counts[j-2][i] += msa->ss->cat_counts[j][i];
	  msa->ss->cat_counts[j][i] = 0;
	}
      }
    }
  }
  msa->ncats = nmod;
  usemods = smalloc(nmod * sizeof(TreeModel*));
  use_nmod = nmod;
  for (i=0; i < nmod; i++) usemods[i] = mods[i];
  
  tm_fit_multi(usemods, use_nmod, &msa, 1, OPT_VERY_HIGH_PREC, NULL, 1);
  sfree(usemods);
}


void bgchmm_compute_emissions(double **emissions, void **models, int nmodels,
			      void *data0, int sample, int length) {
  struct bgchmm_data_struct *data = (struct bgchmm_data_struct*)data0;
  double *temp_emissions;
  int state, i, j, sspos;
  MSA *msa;
  if (sample != 0) 
    die("bgchmm_compute_emissions got sample=%i (should always be 0)\n", sample);
  msa = data->msa;

  temp_emissions = smalloc(msa->ss->ntuples * sizeof(double));
  
  for (state=0; state < nmodels; state++) {
    for (i=0; i < length; i++) emissions[state][i] = NEGINFTY;
    tl_compute_log_likelihood(models[state], msa,
			      NULL, temp_emissions, -1, NULL);
    
    for (j=0; j < msa->length; j++) {
      sspos = msa->ss->tuple_idx[j];
      if (nmodels==4 && (state==2 || state==3) && 
	  data->bgc_informative != NULL && 
	  data->bgc_informative[sspos]==0)
	emissions[state][j] = emissions[state-2][j];
      else emissions[state][j] = temp_emissions[sspos];
    }
  }
  sfree(temp_emissions);
}


int bgchmm_get_obs_idx(void *data0, int i, int j) {
  struct bgchmm_data_struct *data = (struct bgchmm_data_struct*)data0;
  if (i==-1 || j== -1) {
    return data->msa->ss->ntuples;
  }
  if (i != 0) die("bgchmm_get_obs_idx got i=%i (should always be 0)\n", i);
  return data->msa->ss->tuple_idx[j];
}


void bgchmm_set_hmm(HMM *hmm, double bgc_in, double bgc_out, double cons_in, double cons_out) {
  int i, j, from_cons, from_bgc, to_cons, to_bgc, have_bgc, num_states;
  double rate, denom;

  for (i=0; i < hmm->nstates; i++)
    for (j=0; j < hmm->nstates; j++)
      mm_set(hmm->transition_matrix, i, j, 0.0);
  
  have_bgc = (hmm->nstates == 4);
  if (have_bgc) 
    num_states = 4;
  else num_states = 2;

  for (i=0; i < num_states; i++) {
    from_cons = (i%2==1);
    from_bgc = (i >= 2);
    for (j=0; j < num_states; j++) {
      to_cons = (j%2 == 1);
      to_bgc = (j >= 2);
      rate = 1.0;
      if (from_cons) {
	if (to_cons) rate *= (1.0-cons_out);
	else rate *= cons_out;
      } else {  //from is not conserved
	if (to_cons) rate *= cons_in;
	else rate *= (1.0-cons_in);
      }
      
      if (have_bgc) {
	if (from_bgc) {
	  if (to_bgc) rate *= (1.0-bgc_out);
	  else rate *= bgc_out;
	} else { //from is no bgc
	  if (to_bgc) rate *= bgc_in;
	  else rate *= (1.0-bgc_in);
	}
      }
      mm_set(hmm->transition_matrix, i, j, rate);
    }
  }
  if (have_bgc) {
    vec_set(hmm->eq_freqs, 0, cons_out/(cons_out+cons_in)*bgc_out/(bgc_out+bgc_in));
    vec_set(hmm->eq_freqs, 1, cons_in/(cons_out+cons_in)*bgc_out/(bgc_out+bgc_in));
    vec_set(hmm->eq_freqs, 2, cons_out/(cons_out+cons_in)*bgc_in/(bgc_out+bgc_in));
    vec_set(hmm->eq_freqs, 3, cons_in/(cons_out+cons_in)*bgc_in/(bgc_out+bgc_in));
  } else {
    vec_set(hmm->eq_freqs, 0, cons_out/(cons_out+cons_in));
    vec_set(hmm->eq_freqs, 1, cons_in/(cons_out+cons_in));
  }
  vec_copy(hmm->begin_transitions, hmm->eq_freqs);

  if (1) {
    //check and see if transition probs sum to 1
    for (i=0; i < hmm->nstates; i++) {
      denom = 0.0;
      for (j=0; j < hmm->nstates; j++)
	denom += mm_get(hmm->transition_matrix, i, j);
      if (fabs(denom-1.0) > 1.0e-6) 
	die("hmm row %i sum=%e\n", i, denom);
    }
  }

  hmm_reset(hmm);
}


void bgchmm_get_rates(HMM *hmm, double *bgc_in, double *bgc_out, double *cons_in, double *cons_out) {
  double denom;
  int i, have_bgc=-1;
  
  if (hmm->nstates == 4) have_bgc = 1;
  else if (hmm->nstates == 2) have_bgc=0;
  else die("bgc_hmm_get_rates got nstates=%i\n", hmm->nstates);

  if (! have_bgc) {
    if (bgc_in != NULL)  *bgc_in = 0;
    if (bgc_out != NULL) *bgc_out = 0;
    if (cons_in != NULL)
      *cons_in = mm_get(hmm->transition_matrix, 0, 1)/(mm_get(hmm->transition_matrix, 0, 0) + mm_get(hmm->transition_matrix, 0, 1));
    if (cons_out != NULL)
      *cons_out = mm_get(hmm->transition_matrix, 1, 0)/(mm_get(hmm->transition_matrix, 1, 0)+mm_get(hmm->transition_matrix, 1, 1));
  } else {  //have bgc
    if (cons_in != NULL) {
      denom = 0.0;
      for (i=0; i < 4; i++)
	denom += mm_get(hmm->transition_matrix, 0, i);
      *cons_in = (mm_get(hmm->transition_matrix, 0, 1) + mm_get(hmm->transition_matrix, 0, 3))/denom;
    }
    if (cons_out != NULL) {
      denom = 0.0;
      for (i=0; i < 4; i++)
	denom += mm_get(hmm->transition_matrix, 1, i);
      *cons_out = (mm_get(hmm->transition_matrix, 1, 0) + mm_get(hmm->transition_matrix, 1, 2))/denom;
    }
    if (bgc_in != NULL) {
      denom = 0.0;
      for (i=0; i < 4; i++)
	denom += mm_get(hmm->transition_matrix, 0, i);
      *bgc_in = (mm_get(hmm->transition_matrix, 0, 2) + mm_get(hmm->transition_matrix, 0, 3))/denom;
    }
    if (bgc_out != NULL) {
      denom = 0.0;
      for (i=0; i < 4; i++)
	denom += mm_get(hmm->transition_matrix, 2, i);
      *bgc_out = (mm_get(hmm->transition_matrix, 2, 0) + mm_get(hmm->transition_matrix, 2, 1))/denom;
    }
  }

  //use the bgc versions of conserved/nonconserved states to check if we get the same results
  if (have_bgc) {
    double cons_in_test, cons_out_test;
    if (cons_in != NULL) {
      denom = 0.0;
      for (i=0; i < 4; i++)
	denom += mm_get(hmm->transition_matrix, 2, i);
      cons_in_test = (mm_get(hmm->transition_matrix, 2, 1) + mm_get(hmm->transition_matrix,2, 3))/denom;
      if (fabs(cons_in_test - *cons_in) > 1.0e-6) 
	die("Got cons_in=%e, cons_in_test=%e\n", *cons_in, cons_in_test);
    }
    
    if (cons_out != NULL) {
      denom = 0.0;
      for (i=0; i < 4; i++)
	denom += mm_get(hmm->transition_matrix, 3, i);
      cons_out_test = (mm_get(hmm->transition_matrix, 3, 0) + mm_get(hmm->transition_matrix, 3, 2))/denom;
      if (fabs(cons_out_test - *cons_out) > 1.0e-6)
	die("Got cons_out=%e, cons_out_test=%e\n", *cons_out, cons_out_test);
    }
  }
}


// here is estimate transition function.  For now we assume that the 
// only rates we want to estimate are those in/out of bgc states.  
// Everything else is held constant
void bgchmm_estimate_transitions(HMM* hmm, void *data0, double **A) {
  //A holds the counts of transitions from one state to the next
  //need four counts: to bgc from bgc, to bgc from non-bgc, to non-bgc from bgc,
  //to non-bgc from non-bgc

  double counts[2][2], cons_in, cons_out, bgc_in, bgc_out;
  int i, j, from_bgc, to_bgc, from_cons, to_cons, have_bgc=-1,
    numstate;
  struct bgchmm_data_struct *data = (struct bgchmm_data_struct*)data0;
  
  if (! (data->estimate_bgc_target_coverage ||
	 data->estimate_bgc_expected_length ||
	 data->estimate_cons_transitions))
    return;

  if (hmm->nstates == 4) have_bgc = 1;
  else if (hmm->nstates == 2) have_bgc=0;
  else die("bgc_hmm_get_rates got nstates=%i\n", hmm->nstates);

  if (have_bgc) 
    numstate = 4;
  else 
    numstate = 2;

  bgc_in = bgc_out = -1.0;
  if (data->estimate_bgc_target_coverage || data->estimate_bgc_expected_length) {
    for (i=0; i < 2; i++)
      for (j=0; j < 2; j++)
	counts[i][j] = 0.0;
    
    for (i=0; i < numstate; i++) {  //from state
      from_bgc = (i >= 2);
      for (j=0; j < numstate; j++) {
	to_bgc = (j >= 2);
	counts[from_bgc][to_bgc] += A[i][j];
      }
    }
    
    if (data->estimate_bgc_expected_length && !data->estimate_bgc_target_coverage) {
      //use same algorithm as phmm_estim_trans_em_coverage.  But counts are flipped, since here we want expected coverage of bgc state (which is stored in index 1), whereas in phastCons conserved state is state 0
      double a, b, c, mu, nu, nu1, nu2, z, tmp;
      double gamma = data->bgc_target_coverage;
      /* if you take the first derivative wrt nu of the expression inside
	 the argmax and set it to zero, you get a quadratic eqn which
	 can be solved using the quadratic formula */
      z = (1-gamma)/gamma;
      a = z * (counts[0][0] + counts[0][1] + counts[1][0] + counts[1][1]);
      b = -counts[0][1] - counts[1][0] - counts[0][0] - z * (counts[1][1] + counts[0][1] + counts[1][0]);
      c = counts[0][1] + counts[1][0];
      
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
      
      /* double check that derivative is really zero.  NOTE I couldn't get the equivalent expression in phmm_estim_trans_em_coverage to come out to zero all the time.  It was usually close but sometimes had an absolute value > 1.0e-4 (the biggest threshold I tried)  The expression that doesn't work is commented out.  The expression that is used I derived myself... */
      //if (!(fabs(-z*counts[1][1]/(1-z*nu) + (counts[1][0] + counts[0][1])/nu - counts[0][0]/(1-nu)) < 1e-4))
      if (fabs(-counts[0][0]*nu*(gamma - nu + gamma * nu) + (counts[0][1]+counts[1][0])*(1.0-nu)*(gamma - nu + gamma * nu) + counts[1][1]*(gamma-1.0)*(1.0-nu)*nu) >= 1.0e-6)
	die("ERROR phmm_estim_trans_em_coverage: derivate not zero?\n");
      
      mu = z * nu;
      if (!(nu >= 0 && nu <= 1 && mu >= 0 && mu <= 1))
	die("ERROR phmm_estim_trans_em_coverage: mu=%e, nu=%e\n", mu, nu);
      bgc_in = nu;
      bgc_out = mu;
    } else {
      bgc_in = counts[0][1]/(counts[0][0]+counts[0][1]);
      if (!data->estimate_bgc_expected_length)
	bgc_out = 1.0/data->bgc_expected_length;
      else
	bgc_out = counts[1][0]/(counts[1][1]+counts[1][0]);
    }
  } else {
    bgchmm_get_rates(hmm, 
		     bgc_in < 0.0 ? &bgc_in : NULL,
		     bgc_out < 0.0 ? &bgc_out : NULL,
		     NULL, NULL);
  }

  if (data->estimate_cons_transitions) {
    for (i=0; i < 2; i++)
      for (j=0; j < 2; j++)
	counts[i][j] = 0.0;
    for (i=0; i < numstate; i++) {
      from_cons = (i%2==1);
      for (j=0; j < numstate; j++) {
	to_cons = (j%2==1);
	counts[from_cons][to_cons] += A[i][j];
      }
    }
    cons_in = counts[0][1]/(counts[0][0]+counts[0][1]);
    cons_out = counts[1][0]/(counts[1][1]+counts[1][0]);
  }
  else bgchmm_get_rates(hmm, NULL, NULL, &cons_in, &cons_out);

  bgchmm_set_hmm(hmm, bgc_in, bgc_out, cons_in, cons_out);
}

