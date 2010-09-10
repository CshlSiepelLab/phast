/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: phylo_hmm.c,v 1.31 2008-11-12 02:07:59 acs Exp $ */

/** \file phylo_hmm.c
   Code for phylo-HMMs.  Allows for automatic expansion of the state
   space to accommodate features on the reverse strand, and for the
   indel model described in Siepel & Haussler, RECOMB '04.  Also
   allows for cross-product constructions involving functional states
   and rate categories (Siepel & Haussler, RECOMB '03).  */

#include <phylo_hmm.h>
#include <dgamma.h>
#include <sufficient_stats.h>
#include <gap_patterns.h>
#include <tree_likelihoods.h>
#include <em.h>

/* initial values for alpha, beta, tau; possibly should be passed in instead */
#define ALPHA_INIT 0.05
#define BETA_INIT 0.05
#define TAU_INIT 0.15
#define COMPLEX_EPSILON 1e-5

/** Create a new PhyloHmm object. Optionally expands original HMM to
    allow for features on both the positive and negative strands. */
PhyloHmm *phmm_new(HMM *hmm,    /**< HMM.  If indel_mode ==
                                   MISSING_DATA or indel_mode ==
                                   PARAMETERIC, then an HMM for
                                   functional cateogories (with
                                   nstates == cm->ncats + 1) should be
                                   passed in; otherwise (indel_mode ==
                                   NONPARAMETERIC), an hmm for the
                                   full gapped (but not 'reflected')
                                   phylo-HMM should be passed in.
                                   NULL may also be given, in which
                                   case a trivial, single-state HMM
                                   will be created.  Any HMM that is
                                   passed in will be left unchanged (a
                                   copy will be created) */
                   TreeModel **tree_models, 
                                /**< Array of TreeModel objects.
                                   Number of elements is assumed to
                                   equal number of categories
                                   (cm->ncats+1)  */
                   CategoryMap *cm, 
                                /**< CategoryMap.  A copy is created,
                                   orig. isn't touched.  If NULL, a
                                   trivial map is created, with a
                                   separate category for every HMM
                                   state.  Must be non-NULL if
                                   indel_mode != MISSING_DATA. */
                   List *pivot_cats,  
                                /**< Categories (by name or number)
                                   about which to "reflect" the HMM.
                                   Allows for prediction on both
                                   strands (see hmm_reverse_compl).
                                   Pass NULL for no reflection */
                   indel_mode_type indel_mode 
                                /**< How to model indels.  Allowable
                                     values are MISSING_DATA,
                                     PARAMETERIC, and
                                     NONPARAMETERIC. */
                   ) { 
  int max_nstates, s, cat, j, nunspooled_cats;

  PhyloHmm *phmm = smalloc(sizeof(PhyloHmm));
  TreeNode *topology = NULL;

  if (hmm == NULL) 
    hmm = hmm_create_trivial();

  phmm->functional_hmm = hmm;
  if (indel_mode == PARAMETERIC) phmm->hmm = NULL; /* will be created below */
  else phmm->hmm = hmm;         /* hmm == functional_hmm unless
                                   parameteric indel model (or rates
                                   cross) */

  hmm = NULL;                   /* be sure this doesn't get used by
                                   mistake */

  if (cm != NULL)
    phmm->cm = cm_create_copy(cm);
  else {
    if (indel_mode != MISSING_DATA) 
      die("ERROR: must pass non-NULL category map if using indel model");
                                /* would get a little tricky without a cm */
    phmm->cm = cm_create_trivial(phmm->hmm->nstates-1, "model_");
  }

  cm = NULL;                    /* be sure doesn't get used by mistake */
  
  phmm->mods = tree_models;
  phmm->nmods = phmm->cm->ncats + 1;
  phmm->nratecats = 1;
  phmm->indel_mode = indel_mode;
  phmm->autocorr_hmm = NULL;
  phmm->reflected = pivot_cats != NULL;
  phmm->emissions = NULL;
  phmm->forward = NULL;
  phmm->alloc_len = -1;
  phmm->state_pos = phmm->state_neg = NULL;
  phmm->gpm = NULL;
  phmm->T = phmm->t = NULL;
  phmm->em_data = NULL;
  phmm->alpha = NULL;
  phmm->beta = NULL;
  phmm->tau = NULL;
  

  /* make sure tree models all have trees and all have the same number
     of leaves; keep a pointer to a representative tree for use with
     indel model (in this case, topologies must be the same) */
  for (s = 0; s < phmm->nmods; s++) {
    if (phmm->mods[s]->tree == NULL) 
      die("ERROR: tree model #%d for phylo-HMM has no tree.\n", s+1);
    else if (topology == NULL)
      topology = phmm->mods[s]->tree;
    else if (phmm->mods[s]->tree->nnodes != topology->nnodes) 
      die("ERROR: tree models for phylo-HMM have different numbers of nodes.\n");
  }

  /* expand category map for indel model, if necessary.  Also obtain
     mappings between cats and gapcats */
  if (indel_mode != MISSING_DATA) {
    /* start by enlarging catmap to allows for gap cats, and creating
       mappings between cats and gapcats */
    /* need to pass gp_create_gapcats the feature types in question --
       we'll use all except the ones for which gaps are completely
       prohibited (i.e., in all categories in range) */
    List *indel_types = lst_new_ptr(phmm->cm->ncats+1);
    for (cat = 0; cat <= phmm->cm->ncats; ) {
      int allow_gaps = FALSE;
      for (j = cat; !allow_gaps && j <= phmm->cm->ranges[cat]->end_cat_no; j++)
        if (phmm->mods[j]->allow_gaps) allow_gaps = TRUE;
      if (allow_gaps) lst_push_ptr(indel_types, cm_get_feature(phmm->cm, cat));
      cat = phmm->cm->ranges[cat]->end_cat_no + 1;
    }
    phmm->gpm = gp_create_gapcats(phmm->cm, indel_types, topology, 
                                  indel_mode == PARAMETERIC);
    lst_free(indel_types);
  }

  nunspooled_cats = phmm->cm->unspooler == NULL ? phmm->cm->ncats + 1 : 
    phmm->cm->unspooler->nstates_unspooled;

  if (phmm->hmm == NULL)        /* parameteric indel model -- create hmm */
    phmm->hmm = hmm_new_nstates(nunspooled_cats, TRUE, FALSE);

  /* now the number of unspooled categories should equal the number of
     states in the HMM */
  if (phmm->hmm->nstates != nunspooled_cats)
    die("ERROR: number of states in HMM must equal number of site categories (unspooled).\n");

  /* initialize mappings */
  max_nstates = phmm->reflected ? phmm->hmm->nstates * 2 : 
    phmm->hmm->nstates;
  phmm->state_to_mod = smalloc(max_nstates * sizeof(int));
  phmm->state_to_cat = smalloc(max_nstates * sizeof(int));
  phmm->state_to_pattern = smalloc(max_nstates * sizeof(int));
  phmm->reverse_compl = smalloc(max_nstates * sizeof(int));
  phmm->cat_to_states = smalloc((phmm->cm->ncats+1) * sizeof(List*));
  for (cat = 0; cat <= phmm->cm->ncats; cat++) 
    phmm->cat_to_states[cat] = lst_new_int(5);

  /* we need separate initialization logic for the cases with and
     without indel cats.  This is related to the fact that we're using
     the "unspooling" mechanism of the category map to help define the
     indel model.  When there are no indel cats, the HMM states
     correspond to ordinary ("unspooled") categories, but when there
     are indel cats, they correspond to (unspooled) *gap* categories.
     Without indel cats, there is a direct correspondence between
     "spooled" categories and models; with indel cats, multiple
     spooled gap categories map to the same model */

  for (s = 0; s < max_nstates; s++) {
    phmm->reverse_compl[s] = 0;
    phmm->state_to_pattern[s] = -1;

    if (s >= phmm->hmm->nstates) {
      phmm->state_to_cat[s] = phmm->state_to_mod[s] = -1;
      continue; 
    }

    /* in both cases below, 'cat' is the spooled ordinary (non-gap)
       category corresponding to state s */
    if (indel_mode == MISSING_DATA) 
      cat = (phmm->cm->unspooler == NULL ? s : 
             phmm->cm->unspooler->unspooled_to_spooled[s]);
    else {
      int sp_gapcat = (phmm->cm->unspooler == NULL ? s : 
                       phmm->cm->unspooler->unspooled_to_spooled[s]);
      cat = phmm->gpm->gapcat_to_cat[sp_gapcat];
      phmm->state_to_pattern[s] = phmm->gpm->gapcat_to_pattern[sp_gapcat];
    }

    phmm->state_to_cat[s] = phmm->state_to_mod[s] = cat;
    lst_push_int(phmm->cat_to_states[cat], s); 
                                /* FIXME: I don't think this gets
                                   propagated through phmm_reflect_hmm
                                   and phmm_rate_cats yet */
  }

  /* reflect if necessary */
  if (pivot_cats != NULL) 
    phmm_reflect_hmm(phmm, pivot_cats);

  /* initialize alpha, beta, and tau if parameteric indel model */
  if (indel_mode == PARAMETERIC) {
    phmm->alpha = smalloc(phmm->functional_hmm->nstates * sizeof(double));
    phmm->beta = smalloc(phmm->functional_hmm->nstates * sizeof(double));
    phmm->tau = smalloc(phmm->functional_hmm->nstates * sizeof(double));
    for (j = 0; j < phmm->functional_hmm->nstates; j++) {
      phmm->alpha[j] = ALPHA_INIT;
      phmm->beta[j] = BETA_INIT;
      phmm->tau[j] = TAU_INIT;
    }
    phmm_reset(phmm);
  }

  return phmm;
}

/* Reflect PhyloHmm about pivot states corresponding to specified list of
   category names; update all mappings accordingly */
void phmm_reflect_hmm(PhyloHmm *phmm, List *pivot_cats) {
  int ncats = phmm->cm->ncats + 1;
  int mark[ncats];
  int i;
  int *new_to_old;
  HMM *tmp_hmm;
  List *pivot_states = lst_new_int(50);

  for (i = 0; i < ncats; i++) mark[i] = 0;
  mark[0] = 1;                  /* background is implicit */

  if (pivot_cats != NULL) {
    List *cats = cm_get_category_list(phmm->cm, pivot_cats, TRUE);
    for (i = 0; i < lst_size(cats); i++) {
      mark[lst_get_int(cats, i)] = 1;
      //      printf("pivot state %i\n", lst_get_int(cats, i));
    }
    lst_free(cats);
  }

  /* enumerate states corresponding to pivot categories.  This will
     take care of indel categories as well as ordinary unspooled
     categories */
  for (i = 1; i < phmm->hmm->nstates; i++) /* skip 0; already taken care of */
    if (mark[phmm->state_to_cat[i]])
      lst_push_int(pivot_states, i);

  /* now reflect HMM */
  new_to_old = smalloc(phmm->hmm->nstates * 2 * sizeof(int));
  tmp_hmm = hmm_reverse_compl(phmm->hmm, pivot_states, new_to_old); 
  hmm_free(phmm->hmm);
  phmm->hmm = tmp_hmm;

  /* finally, revise mappings */
  for (i = phmm->hmm->nstates-1; i >= 0; i--) {
    /* we use the fact that abs(new_to_old[i]) <= i to do the
       remapping in place */
    phmm->state_to_mod[i] = phmm->state_to_mod[abs(new_to_old[i])];
    phmm->state_to_cat[i] = phmm->state_to_cat[abs(new_to_old[i])];
    phmm->state_to_pattern[i] = phmm->state_to_pattern[abs(new_to_old[i])];
    phmm->reverse_compl[i] = (new_to_old[i] < 0);
  }

  free(new_to_old);
  lst_free(pivot_states);
}

/* create an HMM representing the autocorrelation model of Felsenstein
   and Churchill; HMM object must already be allocated with desired
   number of states */
void phmm_create_autocorr_hmm(HMM *hmm, double lambda) {
  int i, j;

  /* fill markov matrix */
  for (i = 0; i < hmm->nstates; i++) 
    for (j = 0; j < hmm->nstates; j++) 
      mm_set(hmm->transition_matrix, i, j, 
             (1-lambda)/hmm->nstates + (i == j ? lambda : 0));

  /* set eqfreqs and begin transitions (uniform) */
  vec_set_all(hmm->eq_freqs, (double)1/hmm->nstates);
  vec_set_all(hmm->begin_transitions, (double)1/hmm->nstates);

  hmm_reset(hmm);
}


/** Scale tree models and replace HMM by its cross-product with an
   autocorr HMM defined by 'lambda'.  Update all mappings
   appropriately.  Weight-matrix tree models are not scaled.  */
void phmm_rates_cross(PhyloHmm *phmm, 
                                /**< PhyloHmm to be altered */
                      int nratecats, 
                                /**< Number of rate categories. */
                      double lambda,
                                /**< Autocorrelation parameter  */
                      int expand_cats
                                /**< whether to expand the category
                                   map to reflect the rate categories
                                   (useful, e.g., in prediction of
                                   rate categories).  Currently not
                                   allowed with the indel model, or
                                   with any category map having
                                   "conditioned_on" dependencies */
                    ) {

  TreeModel **new_mods;
  double *sconsts;
  int max_nstates, mod, state, ratecat, thismod_new = 0, i, cat;
  int *old_mod_to_new_base;     /* mapping from old model index to
                                   base index in new space */
  int *old_cat_to_new;           /* mapping of old to new category number */

  if (nratecats <= 1) 
    die("ERROR: phmm_rate_cats requires nratecats > 1.\n");

  if (phmm->indel_mode != MISSING_DATA)
    die("ERROR: phmm_rates_cross cannot be used with indel model (parameteric or nonparameteric).\n");

  phmm->nratecats = nratecats;

  /* realloc and init */
  max_nstates = phmm->hmm->nstates * phmm->nratecats;
  phmm->state_to_mod = srealloc(phmm->state_to_mod, max_nstates * sizeof(int));
  phmm->state_to_cat = srealloc(phmm->state_to_cat, max_nstates * sizeof(int));
  phmm->state_to_pattern = srealloc(phmm->state_to_pattern, max_nstates * sizeof(int));
  phmm->reverse_compl = srealloc(phmm->reverse_compl, max_nstates * sizeof(int));

  for (state = phmm->hmm->nstates; state < max_nstates; state++) {
    phmm->state_to_mod[state] = -1;
    phmm->state_to_cat[state] = -1;
    phmm->state_to_pattern[state] = -1;
    phmm->reverse_compl[state] = 0;
  }

  new_mods = smalloc(sizeof(TreeModel*) * phmm->nmods * phmm->nratecats);  
  sconsts = smalloc(phmm->nratecats * sizeof(double));
  old_mod_to_new_base = smalloc(phmm->nmods * sizeof(int));

  /* copy and scale tree models */
  for (mod = 0; mod < phmm->nmods; mod++) {
    if (phmm->mods[mod]->tree == NULL) { /* weight matrix */
      old_mod_to_new_base[mod] = thismod_new;
      new_mods[thismod_new++] = phmm->mods[mod]; /* just reuse same object */
    }
    else {      /* full tree model */
      /* in this case, the scaling constants must be determined from
         the model (and may be different for each model) */
      double freqs[phmm->nratecats];
      if (phmm->mods[mod]->alpha <= 0)
	die("ERROR phmm_rates_cross: alpha=%f, should be > 0\n", 
	    phmm->mods[mod]->alpha);
      DiscreteGamma(freqs, sconsts, phmm->mods[mod]->alpha, 
                    phmm->mods[mod]->alpha, phmm->nratecats, 0);
          
      tm_reinit(phmm->mods[mod], phmm->mods[mod]->subst_mod, 1, 0, NULL, NULL);
      /* we don't want to treat the scaled versions as
         models with rate variation (would be "counting twice") */

      /* now make nratecats scaled copies of model */
      old_mod_to_new_base[mod] = thismod_new;
      for (ratecat = 0; ratecat < phmm->nratecats; ratecat++) {
        new_mods[thismod_new] = tm_create_copy(phmm->mods[mod]);

/*         fprintf(stderr, "Scaling tree model %d by factor %f ...\n", */
/*                 mod+1, sconsts[ratecat]); */ /* remove? */

        tm_scale_branchlens(new_mods[thismod_new], sconsts[ratecat], 1);
        thismod_new++;
      }
    }
  }

  /* create new categories, if necessary */
  if (expand_cats) {
    int old_ncats = phmm->cm->ncats + 1;
    int newcat = old_ncats;
    old_cat_to_new = smalloc(old_ncats * phmm->nratecats * sizeof(int));
    cm_realloc(phmm->cm, old_ncats * phmm->nratecats - 1);
    for (cat = 0; cat < old_ncats; cat++) {
      CategoryRange *r = phmm->cm->ranges[cat];
      int range_size = r->end_cat_no - r->start_cat_no + 1;
      String *oldtype = str_dup(lst_get_ptr(r->feature_types, 0)), *newtype;

      /* update name of orig category with rate 0 */
      if (old_ncats == 1) 
        str_cpy_charstr(lst_get_ptr(r->feature_types, 0), "rate_0"); 
                                /* just replace */
      else str_append_charstr(oldtype, ",rate_0");
      if (phmm->cm->conditioned_on[cat] != NULL)
        die("ERROR: cannot expand cats if conditioned_on != NULL (phmm_rate_cross)\n");
      old_cat_to_new[phmm->nratecats * cat + 0] = cat; /* mapping */

      /* now add nratecats new categories, with names indicating rate cat */
      for (ratecat = 1; ratecat < phmm->nratecats; ratecat++) {
        if (old_ncats == 1) newtype = str_new_charstr("rate_");
        else {
          newtype = str_dup(oldtype);
          str_append_charstr(newtype, ",rate_");
        }
        str_append_int(newtype, ratecat);
        phmm->cm->ranges[newcat] = cm_new_category_range(newtype, newcat, 
                                                         newcat + range_size - 1);
        for (i = 0; i < range_size; i++) {
          phmm->cm->conditioned_on[newcat + i] = NULL;
          if (i > 0) 
            phmm->cm->ranges[newcat + i] = phmm->cm->ranges[newcat];
        }
        old_cat_to_new[phmm->nratecats * cat + ratecat] = newcat;
        newcat += range_size;
      }
      str_free(oldtype);
    }
  }
  else {                        /* just use identity mapping */
    old_cat_to_new = smalloc((phmm->cm->ncats + 1) * sizeof(int));
    for (cat = 0; cat <= phmm->cm->ncats; cat++) old_cat_to_new[cat] = cat;
  }

  /* revise mappings */
  for (state = phmm->hmm->nstates-1; state >= 0; state--) { 
    for (ratecat = phmm->nratecats-1; ratecat >= 0; ratecat--) {
      int new_state = state*phmm->nratecats + ratecat;
      phmm->state_to_mod[new_state] = 
        old_mod_to_new_base[phmm->state_to_mod[state]] + 
        (phmm->mods[phmm->state_to_mod[state]]->tree == NULL ? 0 : ratecat); 
      phmm->state_to_cat[new_state] = 
        old_cat_to_new[phmm->nratecats * phmm->state_to_cat[state] + ratecat];
      phmm->state_to_pattern[new_state] = phmm->state_to_pattern[state];
      phmm->reverse_compl[new_state] = phmm->reverse_compl[state];
      /* we use the fact that new_state > old_state
         except when state == 0 && ratecat == 0 (the last case
         considered) to do the remapping in place */
    }
  }

  /* replace HMM with cross product of (functional) HMM and autocorr HMM */
  phmm->autocorr_hmm = hmm_new_nstates(phmm->nratecats, 1, 0);
  phmm_create_autocorr_hmm(phmm->autocorr_hmm, lambda);
  phmm->hmm = hmm_new_nstates(phmm->nratecats * phmm->functional_hmm->nstates, 1, 0); 
  hmm_cross_product(phmm->hmm, phmm->functional_hmm, phmm->autocorr_hmm);

  /* free memory */
  for (mod = 0; mod < phmm->nmods; mod++)
    if (phmm->mods[mod]->tree != NULL) tm_free(phmm->mods[mod]);
  free(phmm->mods);
  phmm->mods = new_mods;
  phmm->nmods = thismod_new;
  free(old_mod_to_new_base);
  free(old_cat_to_new);
}

/* update cross product HMM using new value of lambda */
void phmm_update_cross_prod(PhyloHmm *phmm, double lambda) {
  phmm_create_autocorr_hmm(phmm->autocorr_hmm, lambda);
  hmm_cross_product(phmm->hmm, phmm->functional_hmm, phmm->autocorr_hmm);
}

void phmm_free(PhyloHmm *phmm) {
  int i;
  for (i = 0; i < phmm->nmods; i++) tm_free(phmm->mods[i]);
  free(phmm->mods);

  if (phmm->emissions != NULL) {
    for (i = 0; i < phmm->hmm->nstates; i++) 
      if (phmm->state_pos[phmm->state_to_mod[i]] == i ||
          phmm->state_neg[phmm->state_to_mod[i]] == i || 
          phmm->state_to_pattern[i] >= 0)
        free(phmm->emissions[i]);
    free(phmm->emissions); free(phmm->state_pos); free(phmm->state_neg);
  }

  if (phmm->forward != NULL) {
    for (i = 0; i < phmm->hmm->nstates; i++) free(phmm->forward[i]);
    free(phmm->forward);
  }

  free(phmm->state_to_mod);
  free(phmm->state_to_cat);
  free(phmm->state_to_pattern);
  free(phmm->reverse_compl);
  for (i = 0; i <= phmm->cm->ncats; i++)
    lst_free(phmm->cat_to_states[i]);
  free(phmm->cat_to_states);
  cm_free(phmm->cm);
  if (phmm->T != NULL) {
    for (i = 0; i < phmm->functional_hmm->nstates; i++) {
      free(phmm->T[i]);
      free(phmm->t[i]);
    }
    free(phmm->T);
    free(phmm->t);
  }
  if (phmm->gpm != NULL) gp_free_map(phmm->gpm);
  if (phmm->functional_hmm != phmm->hmm) hmm_free(phmm->functional_hmm);
  if (phmm->autocorr_hmm != NULL) hmm_free(phmm->autocorr_hmm);
  if (phmm->alpha != NULL) free(phmm->alpha);
  if (phmm->beta != NULL) free(phmm->beta);
  if (phmm->tau != NULL) free(phmm->tau);
  if (phmm->em_data != NULL) {
    if (phmm->em_data->H != NULL) 
      mat_free(phmm->em_data->H);
    free(phmm->em_data);
  }
  hmm_free(phmm->hmm);
  free(phmm);
}

/** Compute emissions for given PhyloHmm and MSA.  Preprocessor for
    phmm_viterbi_features, phmm_posterior_probs, and phmm_lnl
    (often only needs to be run once). */
void phmm_compute_emissions(PhyloHmm *phmm,
                                /**< Initialized PhyloHmm */
                            MSA *msa,
                                /**< Source alignment */
                            int quiet
                                /**< Determins whether progress is
                                   reported to stderr */
                            ) {

  int i, mod, j;
  MSA *msa_compl = NULL;
  int new_alloc = (phmm->emissions == NULL); 
  /* allocate new memory if emissions is NULL; otherwise reuse */ 

  if (new_alloc) {
    if (phmm->emissions == NULL)
      phmm->emissions = smalloc(phmm->hmm->nstates * sizeof(double*));  
    phmm->alloc_len = msa->length;
  }
  if (phmm->alloc_len < msa->length)
    die("ERROR phmm_compute_emissions: phmm->alloc_len (%i) < msa->length (%i)\n",
	phmm->alloc_len, msa->length);

  /* if HMM is reflected, we need the reverse complement of the
     alignment as well */
  if (phmm->reflected) {          
    int idx1, idx2;
    msa_compl = msa_create_copy(msa, 0);
    msa_reverse_compl(msa_compl);

    /* we actually want to keep the indexing of the forward strand; to
       save code, we'll just *reverse* the reverse complement */
    for (idx1 = 0, idx2 = msa_compl->length-1; 
         idx1 < idx2; idx1++, idx2--) {
      int tmp = msa_compl->ss->tuple_idx[idx2];
      msa_compl->ss->tuple_idx[idx2] = msa_compl->ss->tuple_idx[idx1];
      msa_compl->ss->tuple_idx[idx1] = tmp;
    }

    /* get rid of the sequences! they'll be wrong! */
    if (msa_compl->seqs != NULL) {
      for (i = 0; i < msa_compl->nseqs; i++) free(msa_compl->seqs[i]);
      free(msa_compl->seqs);
      msa_compl->seqs = NULL;
    }
  }

  /* set up mapping from model/strand to first associated state
     (allows phmm->emissions to be computed only once for each
     model/strand pair) */
  if (new_alloc) {
    phmm->state_pos = smalloc(phmm->nmods * sizeof(int));
    phmm->state_neg = smalloc(phmm->nmods * sizeof(int));
  }
  for (i = 0; i < phmm->nmods; i++) 
    phmm->state_pos[i] = phmm->state_neg[i] = -1;

  for (i = 0; i < phmm->hmm->nstates; i++) {
    if (!quiet) {
      fprintf(stderr, "Computing emission probs (state %d, cat %d, mod %d",
              i, phmm->state_to_cat[i], phmm->state_to_mod[i]);
      if (phmm->state_to_pattern[i] != -1) 
        fprintf(stderr, ", pattern %d", phmm->state_to_pattern[i]);
      if (phmm->reflected) 
        fprintf(stderr, ", strand %c", phmm->reverse_compl[i] ? '-' : '+');
      fprintf(stderr, ")...\n");
    }
            
    /* reuse already computed values if possible */
    mod = phmm->state_to_mod[i];
    if (!phmm->reverse_compl[i] && phmm->state_pos[mod] != -1)
      phmm->emissions[i] = phmm->emissions[phmm->state_pos[mod]]; /* saves memory */
    else if (phmm->reverse_compl[i] && phmm->state_neg[mod] != -1)
      phmm->emissions[i] = phmm->emissions[phmm->state_neg[mod]];
    else {
      if (new_alloc)
	phmm->emissions[i] = smalloc(msa->length * sizeof(double));

      tl_compute_log_likelihood(phmm->mods[mod], 
                                phmm->reverse_compl[i] ? msa_compl : msa,
                                phmm->emissions[i], -1, NULL);
      if (!phmm->reverse_compl[i]) phmm->state_pos[mod] = i;
      else phmm->state_neg[mod] = i;            
    }
  }
  if (msa_compl != NULL) msa_free(msa_compl);

  /* finally, adjust for indel model, if necessary */
  if (phmm->indel_mode != MISSING_DATA) {
    int *matches = smalloc(msa->ss->ntuples * sizeof(int));
                                /* msa->ss should exist
                                   (tl_compute_log_likelihood) */

    if (!quiet)
      fprintf(stderr, "Adjusting emission probs according to gap patterns...\n");
    for (i = phmm->hmm->nstates - 1; i >= 0; i--) {
                                /* by going backwards, we ensure that
                                   the "base" state (with gap pattern
                                   == 0) is visited last */
      if (phmm->state_to_pattern[i] >= 0) {
        double *orig_emissions = phmm->emissions[i];

        if (phmm->state_to_pattern[i] > 0)
          phmm->emissions[i] = smalloc(msa->length * sizeof(double));
                                /* otherwise, use the array already
                                   allocated */

        gp_tuple_matches_pattern(phmm->gpm, msa, phmm->state_to_pattern[i],
                                 matches);

        for (j = 0; j < msa->length; j++) 
          phmm->emissions[i][j] = 
            (matches[msa->ss->tuple_idx[j]] ? orig_emissions[j] : NEGINFTY);
      }
    }
    free(matches);
  }
}

/** Run the Viterbi algorithm and return a set of predictions.
    Emissions must have already been computed (see
    phmm_compute_emissions) */
GFF_Set* phmm_predict_viterbi(PhyloHmm *phmm, 
                                /**< PhyloHmm object */
                              char *seqname,
                                /**< seqname for feature set (e.g.,
                                   "chr1") */
                              char *grouptag,
                                /**< tag to use for groups (e.g.,
                                   "exon_id", "transcript_id"); if
                                   NULL, a default will be used */
                              char *idpref,
                                /**< prefix for assigned ids (e.g.,
                                   "chr1.15") (may be NULL) */
                              List *frame
                                /**< names of features for which to obtain
                                   frame (NULL to ignore) */
                               ) {
  int *path = (int*)smalloc(phmm->alloc_len * sizeof(int));
  GFF_Set *retval;

  if (phmm->emissions == NULL)
    die("ERROR: emissions required for phmm_viterbi_features.\n");
          
  hmm_viterbi(phmm->hmm, phmm->emissions, phmm->alloc_len, path);

  retval = cm_labeling_as_gff(phmm->cm, path, phmm->alloc_len, 
                              phmm->state_to_cat, 
                              phmm->reverse_compl, seqname, "PHAST",  
                              frame, grouptag, idpref);
  free(path);

  return retval;
}

/** Run the Viterbi algorithm and return a set of features indicating
    regions in which the Viterbi path remains in states corresponding
    to any of the specified categories.  Emissions must have already
    been computed (see phmm_compute_emissions) */
GFF_Set* phmm_predict_viterbi_cats(PhyloHmm *phmm, 
                                   /**< PhyloHmm object */
                                   List *cats,
                                   /**< categories of interest, by
                                      name or number */
                                   char *seqname,
                                   /**< seqname for feature set (e.g.,
                                      "chr1") */
                                   char *grouptag,
                                   /**< tag to use for groups (e.g.,
                                      "exon_id", "transcript_id"); if
                                      NULL, a default will be used */
                                   char *idpref,
                                   /**< prefix for assigned ids (e.g.,
                                      "chr1.15") (may be NULL) */
                                   List *frame,
                                   /**< features for which to obtain
                                      frame (NULL to ignore) */
                                   char *new_type
                                   /**< replace type of each retained
                                      feature with this string if
                                      non-NULL (old types may no
                                      longer make sense, because of
                                      merging) */
                                   ) {
  int i;
  GFF_Feature *lastkeeper = NULL;
  List *types, *keepers, *catnos;
  GFF_Set *retval = phmm_predict_viterbi(phmm, seqname, grouptag, 
                                         idpref, frame);

  /* do this way to allow input to be numbers or names */
  catnos = cm_get_category_list(phmm->cm, cats, 1);
  types = cm_get_features(phmm->cm, catnos);
  lst_free(catnos);

  /* filter out unwanted types */
  gff_filter_by_type(retval, types, FALSE, NULL);
  lst_free(types);

  /* now merge adjacent features */
  keepers = lst_new_ptr(lst_size(retval->features));
  for (i = 0; i < lst_size(retval->features); i++) {
    GFF_Feature *f = lst_get_ptr(retval->features, i);
    if (lastkeeper != NULL && f->start == lastkeeper->end+1) {
      lastkeeper->end = f->end;
      gff_free_feature(f);
    }
    else {
      lst_push_ptr(keepers, f);
      lastkeeper = f;
      if (new_type != NULL) str_cpy_charstr(f->feature, new_type);
    }
  }
  lst_free(retval->features);
  retval->features = keepers;
  return retval;
}

/** Compute and return log likelihood.  Uses forward algorithm.
    Emissions must have already been computed (see
    phmm_compute_emissions) */
double phmm_lnl(PhyloHmm *phmm) {
  double **forward = smalloc(phmm->hmm->nstates * sizeof(double*));
  int i;
  double logl;

  if (phmm->emissions == NULL)
    die("ERROR: emissions required for phmm_lnl.\n");
          
  for (i = 0; i < phmm->hmm->nstates; i++)
    forward[i] = (double*)smalloc(phmm->alloc_len * sizeof(double));
  logl = hmm_forward(phmm->hmm, phmm->emissions, 
                     phmm->alloc_len, forward);
  for (i = 0; i < phmm->hmm->nstates; i++) free(forward[i]);
  free(forward);
  return logl * log(2); /* convert to natural log */
}

/** Computes posterior probabilities for a PhyloHmm.  Emissions must
    have already been computed (see phmm_compute_emissions).  Returns 
    log likelihood.  */
double phmm_postprobs(PhyloHmm *phmm, double **post_probs) {
  if (phmm->emissions == NULL)
    die("ERROR: emissions required for phmm_posterior_probs.\n");

  return hmm_posterior_probs(phmm->hmm, phmm->emissions, phmm->alloc_len, 
                             post_probs) * log(2);
                                /* convert to natural log */          
}

/** Computes and returns an array of length phmm->alloc_len
    representing the marginal posterior prob at every site, summed
    over states corresponding to categories in the specified list.
    Emissions must have already been computed (see
    phmm_compute_emissions).  Computes log likelihood as a side
    effect. */
double* phmm_postprobs_cats(PhyloHmm *phmm, 
                                /**< PhyloHmm object */
                              List *cats, 
                                /**< Categories of interest, by name
                                   or number */
                              double *lnl
                                /**< if non-NULL, will point to log
                                   likelihood on return */
                              ) {
  int i, j, state;
  double **pp = smalloc(phmm->hmm->nstates * sizeof(double*));
  double *retval = smalloc(phmm->alloc_len * sizeof(double));
  double l;
  List *states = lst_new_int(phmm->hmm->nstates), *catnos;
  int docat[phmm->cm->ncats+1];

  /* get states corresponding to specified cats */
  catnos = cm_get_category_list(phmm->cm, cats, 1);
  for (i = 0; i <= phmm->cm->ncats; i++) docat[i] = 0;
  for (i = 0; i < lst_size(catnos); i++) docat[lst_get_int(catnos, i)] = 1;
  for (i = 0; i < phmm->hmm->nstates; i++)
    if (docat[phmm->state_to_cat[i]]) lst_push_int(states, i);
                                               
  /* only allocate memory for states of interest; NULLs for the others
     will cause hmm_postprobs to ignore them */
  for (i = 0; i < phmm->hmm->nstates; i++) pp[i] = NULL;
  for (i = 0; i < lst_size(states); i++) {
    state = lst_get_int(states, i);
    if (pp[state] == NULL)
      pp[state] = smalloc(phmm->alloc_len * sizeof(double));
  }

  l = phmm_postprobs(phmm, pp);

  for (j = 0; j < phmm->alloc_len; j++) retval[j] = 0;
  for (i = 0; i < lst_size(states); i++) {
    state = lst_get_int(states, i);
    for (j = 0; j < phmm->alloc_len; j++) retval[j] += pp[state][j];
  }
      
  if (lnl != NULL) *lnl = l;

  return retval;
}

/* wrapper for hmm_forward used by phmm_fit_lambda */
double fit_lambda_lnl(double lambda, void *data) {
  PhyloHmm *phmm = data;
  if (lambda < 0 || lambda > 1) return INFTY;
  phmm_update_cross_prod(phmm, lambda);
  return log(2) * -hmm_forward(phmm->hmm, phmm->emissions, 
                               phmm->alloc_len, phmm->forward);
}

/* returns log likelihood */
double phmm_fit_lambda(PhyloHmm *phmm, double *lambda, FILE *logf) {
  double neglnl, ax, bx, cx, fa, fb, fc;
  int i;

  /* allocate memory for forward alg */
  if (phmm->forward == NULL) {  /* otherwise assume already alloc */
    phmm->forward = smalloc(phmm->hmm->nstates * sizeof(double*));
    for (i = 0; i < phmm->hmm->nstates; i++)
      phmm->forward[i] = smalloc(phmm->alloc_len * sizeof(double));
  }

  /* start with a range of 0.2 around the starting value of lambda;
     seems to be a reasonable heuristic (mnbrak can go outside this
     range) */
  bx = min(*lambda + 0.1, .97);  /* don't get too close to boundary */
  ax = bx - 0.2; 
  mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, fit_lambda_lnl, phmm, logf);
  neglnl = opt_brent(ax, bx, cx, fit_lambda_lnl, 5e-3, lambda, phmm, logf);

  return -neglnl;
}

/** Score a set of predicted features using log odds scoring. */
void phmm_score_predictions(PhyloHmm *phmm, 
                                /**< PhyloHmm object */
                            GFF_Set *preds, 
                                /**< Predictions to score */
                            List *score_cats, 
                                /**< Categories to score (specify by
                                   name or number) */
                            List *helper_cats,  
                                /**< Secondary categories to be
                                   included in scoring if adjacent to
                                   main cats.  Pass NULL if none. */
                            List *null_cats,
                                /**< Categories in null model.  Pass
                                   NULL to use all categories not in
                                   score_cats or helper_cats */
                            int score_only_score_cats
                                /**< Restrict scoring to features
                                   matching score_cats; if FALSE,
                                   all features are scored */
                            )  {

  int i, j, cat, ncats, nscore, state;
  List *cats, *score_states, *null_states, *score_types;
  List **cat_to_states;
  int *is_scored;

  /* convert lists of cat names to lists of states */

  /* first create inverse mapping from category numbers to lists of
     states.  FIXME: should have this available (partially
     implemented)  */
  ncats = phmm->cm->ncats + 1;
  cat_to_states = smalloc(ncats * sizeof(List*));
  for (i = 0; i < ncats; i++) 
    cat_to_states[i] = lst_new_int(2*phmm->hmm->nstates/ncats);
  for (i = 0; i < phmm->hmm->nstates; i++) 
    lst_push_int(cat_to_states[phmm->state_to_cat[i]], i);

  /* now fill in state lists */
  nscore = lst_size(score_cats) + 
    (helper_cats != NULL ? lst_size(helper_cats) : 0);
  score_states = lst_new_int(nscore * 10);
  null_states = lst_new_int(20);
  is_scored = smalloc(phmm->hmm->nstates * sizeof(int));
  for (i = 0; i < phmm->hmm->nstates; i++) is_scored[i] = 0;

  cats = cm_get_category_list(phmm->cm, score_cats, 1);
  for (i = 0; i < lst_size(cats); i++) {
    cat = lst_get_int(cats, i);
    for (j = 0; j < lst_size(cat_to_states[cat]); j++) {
      state = lst_get_int(cat_to_states[cat], j);
      lst_push_int(score_states, state);
      is_scored[state] = 1;
    }
  }
  score_types = cm_get_features(phmm->cm, cats);
  lst_free(cats);

  if (helper_cats != NULL) {
    cats = cm_get_category_list(phmm->cm, helper_cats, 1);
    for (i = 0; i < lst_size(cats); i++) {
      cat = lst_get_int(cats, i);
      for (j = 0; j < lst_size(cat_to_states[cat]); j++) {
        state = lst_get_int(cat_to_states[cat], j);
        lst_push_int(score_states, state);
        is_scored[state] = 1;
      }
    }
    lst_free(cats);
  }

  /* null states */
  if (null_cats == NULL) {      /* assume everything not scored */
    for (i = 0; i < phmm->hmm->nstates; i++)
      if (!is_scored[i]) lst_push_int(null_states, i);
  }
  else {
    cats = cm_get_category_list(phmm->cm, null_cats, 1);
    for (i = 0; i < lst_size(cats); i++) {
      cat = lst_get_int(cats, i);
      for (j = 0; j < lst_size(cat_to_states[cat]); j++)
        lst_push_int(null_states, lst_get_int(cat_to_states[cat], j));
    }
    lst_free(cats);
  }

  /* now score each feature */
  for (i = 0; i < lst_size(preds->features); i++) {
    GFF_Feature *feat = lst_get_ptr(preds->features, i);
    int is_score_cat = str_in_list(feat->feature, score_types);
    if (!score_only_score_cats || is_score_cat) {
      int start = feat->start;
      int end = feat->end;

      /* if score cat, extend range as far as possible using helper
         cats */
      if (is_score_cat && helper_cats != NULL) {
        for (j = i-1; j >= 0; j--) {
          GFF_Feature *prev_feat = lst_get_ptr(preds->features, j);
          if (str_in_list(prev_feat->feature, helper_cats) && 
              prev_feat->end == start - 1) 
            start = prev_feat->start;
          else break;
        }
        for (j = i+1; j < lst_size(preds->features); j++) {
          GFF_Feature *next_feat = lst_get_ptr(preds->features, j);
          if (str_in_list(next_feat->feature, helper_cats) && 
              next_feat->start == end + 1) 
            end = next_feat->end;
          else break;
        }
      }

      /* score from start to end */
      feat->score = 
        hmm_log_odds_subset(phmm->hmm, phmm->emissions, score_states, 
                            null_states, start - 1, end - start + 1);

      feat->score_is_null = 0;
    }
  }

  lst_free(score_states);
  lst_free(null_states);
  lst_free(score_types);
  for (i = 0; i < ncats; i++) lst_free(cat_to_states[i]);
  free(cat_to_states);
  free(is_scored);
}

/** Add specified "bias" to log transition probabilities from
   designated background categories to non-background categories, then
   renormalize.  Provides a simple "knob" for controlling the
   sensitivity/specificity tradeoff.  Bias may be negative (increases
   specificity) or positive (increases sensitivity). */
void phmm_add_bias(PhyloHmm *phmm, List *backgd_cat_names, double bias) {
  double is_backgd_cat[phmm->cm->ncats+1];
  int i, j;
  List *backgd_cat_nos;
  double bias_factor = exp(bias); /* multiplicative factor for
                                     transition probs implied by
                                     bias */

  for (i = 0; i <= phmm->cm->ncats; i++) is_backgd_cat[i] = 0;
  backgd_cat_nos = cm_get_category_list(phmm->cm, backgd_cat_names, 0); 
  for (i = 0; i < lst_size(backgd_cat_nos); i++) 
    is_backgd_cat[lst_get_int(backgd_cat_nos, i)] = 1;
  lst_free(backgd_cat_nos);

  for (i = 0; i < phmm->hmm->nstates; i++)
    if (is_backgd_cat[phmm->state_to_cat[i]])
      for (j = 0; j < phmm->hmm->nstates; j++)
        if (!is_backgd_cat[phmm->state_to_cat[j]] && 
            mm_get(phmm->hmm->transition_matrix, i, j) != 0) 
          mm_set(phmm->hmm->transition_matrix, i, j, 
                 bias_factor * mm_get(phmm->hmm->transition_matrix, i, j));
  
  hmm_renormalize(phmm->hmm);
}

/* special-purpose logging function for phmm_fit_em */
void phmm_log_em(FILE *logf, double logl, HMM *hmm, void *data, 
                 int show_header) {
  PhyloHmm *phmm = data;
  int i, j;

  if (show_header) {
    fprintf(logf, "\nlnl\t");
    for (i = 0; i < phmm->functional_hmm->nstates; i++) {
      for (j = 0; j < phmm->functional_hmm->nstates; j++) {
        if (i == j) continue;
        fprintf(logf, "(%d,%d)\t", i, j);
      }
    }
    if (phmm->indel_mode == PARAMETERIC) 
      for (i = 0; i < phmm->functional_hmm->nstates; i++) 
        fprintf(logf, "alpha[%d]\tbeta[%d]\ttau[%d]\t", i, i, i);
    fprintf(logf, "\n");
  }

  fprintf(logf, "%f\t", logl * log(2));
  for (i = 0; i < phmm->functional_hmm->nstates; i++) {
    for (j = 0; j < phmm->functional_hmm->nstates; j++) {
      if (i == j) continue;
      fprintf(logf, "%f\t", mm_get(phmm->functional_hmm->transition_matrix, i, j));
    }
  }
  if (phmm->indel_mode == PARAMETERIC) 
    for (i = 0; i < phmm->functional_hmm->nstates; i++) 
      fprintf(logf, "%f\t%f\t%f\t", phmm->alpha[i], phmm->beta[i], phmm->tau[i]);
  fprintf(logf, "\n");
  fflush(logf);
}

/* The functions below are used by phmm_fit_em; they are phylo-HMM
   specific implementations of generic routines required by
   hmm_train_by_em.  Currently, a single training alignment is assumed
   (i.e., no PooledMSA)  */

/* wrapper for phmm_cmopute_emissions for use in EM */
void phmm_compute_emissions_em(double **emissions, void **models, int nmodels,
                               void *data, int sample, int length) {
  PhyloHmm *phmm = (PhyloHmm*)data;
  phmm_compute_emissions(phmm, phmm->em_data->msa, TRUE);
}

/* re-estimate phylogenetic models based on expected counts */
void phmm_estim_mods_em(void **models, int nmodels, void *data, 
                        double **E, int nobs, FILE *logf) {

  /* FIXME: what about when multiple states per model?  Need to
     collapse sufficient stats.  Could probably be done
     generally... */

  int k, obsidx;
  Vector *params;
  PhyloHmm *phmm = (PhyloHmm*)data;

  if (phmm->em_data->msa->ss == NULL) {
    phmm->em_data->msa->ncats = phmm->nmods - 1;   /* ?? */
    ss_from_msas(phmm->em_data->msa, phmm->mods[0]->order+1, TRUE, 
                 NULL, NULL, NULL, -1);
  }
  else if (phmm->em_data->msa->ncats != phmm->nmods - 1 ||
           phmm->em_data->msa->ss->cat_counts == NULL) {
    phmm->em_data->msa->ncats = phmm->nmods - 1;   /* ?? */
    ss_realloc(phmm->em_data->msa, phmm->em_data->msa->ss->tuple_size, 
               phmm->em_data->msa->ss->ntuples, TRUE, TRUE);
  }

  for (k = 0; k < phmm->nmods; k++) {
    params = tm_params_new_init_from_model(phmm->mods[k]);
    for (obsidx = 0; obsidx < nobs; obsidx++) 
      phmm->em_data->msa->ss->cat_counts[k][obsidx] = E[k][obsidx];
    msa_get_base_freqs_tuples(phmm->em_data->msa, phmm->mods[k]->backgd_freqs, 
                              phmm->mods[k]->order+1, k);
                                /* need to reestimate background
                                   freqs, using new category counts */

    /* FIXME: need to use state_to_cat, etc. in deciding which categories to use */

    tm_fit(phmm->mods[k], phmm->em_data->msa, params, k, OPT_HIGH_PREC, logf);
    vec_free(params); 
  }

  if (phmm->indel_mode == PARAMETERIC)
    phmm_set_branch_len_factors(phmm);
}

/* return observation index associated with given position, here a tuple index */
int phmm_get_obs_idx_em(void *data, int sample, int position) {
  MSA *msa = ((PhyloHmm*)data)->em_data->msa;
  if (sample == -1 || position == -1) 
    return msa->ss->ntuples;
  return msa->ss->tuple_idx[position];
}

/** General routine to estimate the parameters of a phylo-HMM by EM.
   Can be used with or without the indel model, and for estimation of
   transition params only or transition params and tree models.
   Returns ln likelihood. */
double phmm_fit_em(PhyloHmm *phmm, 
                   MSA *msa,     /** NULL means not to estimate tree
                                    models (emissions must be
                                    precomputed) */
                   int fix_functional,
                                /**< Whether to fix transition
                                   parameters of functional HMM */
                   int fix_indel,
                                /**< Whether to fix indel parameters */
                   FILE *logf
                   ) {
  double retval;

  if (msa == NULL && phmm->emissions == NULL)
    die("ERROR (phmm_fit_em): emissions must be precomputed if not estimating tree models.\n");

  phmm->em_data = smalloc(sizeof(EmData));
  phmm->em_data->msa = msa;
  phmm->em_data->fix_functional = fix_functional;
  phmm->em_data->fix_indel = fix_indel;
  phmm->em_data->H = NULL;

  if (msa != NULL)              /* estimating tree models */
    retval = hmm_train_by_em(phmm->hmm, phmm->mods, phmm, 1, 
                             &phmm->alloc_len, NULL, 
                             phmm_compute_emissions_em, phmm_estim_mods_em,
                             phmm_estim_trans_em, phmm_get_obs_idx_em, 
                             phmm_log_em, phmm->emissions, logf);

  else                          /* not estimating tree models */
    retval = hmm_train_by_em(phmm->hmm, phmm->mods, phmm, 1, 
                             &phmm->alloc_len, NULL, NULL, NULL,
                             phmm_estim_trans_em, NULL,
                             phmm_log_em, phmm->emissions, logf);

  return log(2) * retval;
}

/* Set up phmm->T and phmm->t, the arrays of branch lengths and branch
   length sums that are used in the parameteric indel model */
void phmm_set_branch_len_factors(PhyloHmm *phmm) {
  int i, j;

  if (phmm->T == NULL) {
    phmm->T = smalloc(phmm->functional_hmm->nstates * sizeof(void*));
    phmm->t = smalloc(phmm->functional_hmm->nstates * sizeof(void*));
    for (i = 0; i < phmm->functional_hmm->nstates; i++) {
      phmm->T[i] = smalloc(phmm->gpm->ngap_patterns * sizeof(double));
      phmm->t[i] = smalloc(phmm->gpm->ngap_patterns * sizeof(double));
    }
  }

  for (i = 0; i < phmm->functional_hmm->nstates; i++) {
    phmm->T[i][0] = tr_total_len(phmm->mods[i]->tree);
    phmm->T[i][phmm->gpm->ngap_patterns - 1] = -1;
    phmm->t[i][0] = -1;
    phmm->t[i][phmm->gpm->ngap_patterns - 1] = -1;

    for (j = 1; j <= 2*phmm->gpm->nbranches; j++) { /* exclude null and
                                                       complex patterns */
      TreeNode *n = lst_get_ptr(phmm->mods[i]->tree->nodes, 
                                phmm->gpm->pattern_to_node[j]);
                                /* n is the node below the branch
                                   associated with gap pattern j */
      phmm->t[i][j] = n->dparent;
      phmm->T[i][j] = phmm->T[i][0] - phmm->t[i][j];
    }
  }
}

/* Reset HMM of phylo-HMM, allowing for possible changes in indel
   parameters */ 
void phmm_reset(PhyloHmm *phmm) {
  int i, j, cat_i, pat_i, cat_j, pat_j;
  pattern_type pat_i_type, pat_j_type;
  double val;

  if (phmm->indel_mode == PARAMETERIC) {

    phmm_set_branch_len_factors(phmm);

    for (i = 0; i < phmm->hmm->nstates; i++) {
      cat_i = phmm->state_to_cat[i];
      pat_i = phmm->state_to_pattern[i];
      pat_i_type = gp_pattern_type(phmm->gpm, pat_i);

      for (j = 0; j < phmm->hmm->nstates; j++) {
        cat_j = phmm->state_to_cat[j];
        pat_j = phmm->state_to_pattern[j];
        pat_j_type = gp_pattern_type(phmm->gpm, pat_j);

        val = mm_get(phmm->functional_hmm->transition_matrix, cat_i, cat_j);

        if (pat_i_type == COMPLEX_PATTERN && pat_j_type == COMPLEX_PATTERN)
          val *= (1 - (phmm->gpm->ngap_patterns - 1) * COMPLEX_EPSILON);

        else if (pat_i_type == COMPLEX_PATTERN || pat_j_type == COMPLEX_PATTERN) 
          val *= COMPLEX_EPSILON;

        else if (pat_i == pat_j) {
          if (pat_i == 0)
            val *= (1 - COMPLEX_EPSILON - phmm->T[cat_j][0] * 
                    (phmm->alpha[cat_j] + phmm->beta[cat_j]));
          else if (pat_i_type == INSERTION_PATTERN)
            val *= (1 - COMPLEX_EPSILON 
                    - phmm->tau[cat_j] * phmm->alpha[cat_j] * 
                    phmm->T[cat_j][pat_j] 
                    - phmm->tau[cat_j] * phmm->beta[cat_j] * 
                    phmm->T[cat_j][0] 
                    - phmm->tau[cat_j]);
          else if (pat_i_type == DELETION_PATTERN)
            val *= (1 - COMPLEX_EPSILON 
                    - phmm->tau[cat_j] * phmm->alpha[cat_j] * 
                    phmm->T[cat_j][0] 
                    - phmm->tau[cat_j] * phmm->beta[cat_j] * 
                    phmm->T[cat_j][pat_j] 
                    - phmm->tau[cat_j]);
        }

        else {                  /* pat_i != pat_j, neither complex */
          if (pat_i_type != NULL_PATTERN) 
            val *= phmm->tau[cat_j];

          if (pat_j_type == INSERTION_PATTERN) 
            val *= phmm->alpha[cat_j] * phmm->t[cat_j][pat_j];

          else if (pat_j_type == DELETION_PATTERN) 
            val *= phmm->beta[cat_j] * phmm->t[cat_j][pat_j];
        }

        mm_set(phmm->hmm->transition_matrix, i, j, val);
      }

      vec_set(phmm->hmm->begin_transitions, i, 
                     vec_get(phmm->functional_hmm->begin_transitions, 
                                    cat_i) * 
                     1.0/phmm->gpm->ngap_patterns);
                                /* assume uniform distrib. for gap
                                   patterns */
    }
  }

  hmm_reset(phmm->hmm);
}

/* Function to be maximized when estimating indel params by EM --
   portion of expected log likelihood that depends on the indel
   parameters for a given functional category (defined by
   phmm->current_dest_cat) */
double indel_max_function(Vector *params, void *data) {
  IndelEstimData *ied = data;
  int pat, j = ied->current_dest_cat;
  double alpha_j = vec_get(params, 0), beta_j = vec_get(params, 1),
    tau_j = vec_get(params, 2);
  double retval = ied->u_alpha[j] * log(alpha_j) + 
    ied->u_beta[j] * log(beta_j) + ied->u_tau[j] * log(tau_j) +
    ied->u_self[j][0] * log(1 - COMPLEX_EPSILON -
                            (alpha_j + beta_j) * ied->T[j][0]);

  for (pat = 1; pat < ied->gpm->ngap_patterns - 1; pat++) {
    pattern_type type = gp_pattern_type(ied->gpm, pat);
    if (type == INSERTION_PATTERN)
      retval += ied->u_self[j][pat] * 
        log(1 - COMPLEX_EPSILON - tau_j * alpha_j * ied->T[j][pat] 
            - tau_j * beta_j * ied->T[j][0] - tau_j);
    else {
      if (type != DELETION_PATTERN)
	die("ERROR indel_max_function, got unknown type (%i)\n", type);
      retval += ied->u_self[j][pat] * 
        log(1 - COMPLEX_EPSILON - tau_j * alpha_j * ied->T[j][0] 
            - tau_j * beta_j * ied->T[j][pat] - tau_j);
    }
  }
  return -retval;               /* negate for minimization */
}

/* Gradient of indel_max_function */
void indel_max_gradient(Vector *grad, Vector *params,void *data, 
                          Vector *lb, Vector *ub) {
  IndelEstimData *ied = data;
  int pat, j = ied->current_dest_cat;
  double alpha_j = vec_get(params, 0), beta_j = vec_get(params, 1),
    tau_j = vec_get(params, 2);
  double grad_alpha_j = ied->u_alpha[j] / alpha_j 
    - ied->u_self[j][0] * ied->T[j][0] / 
    (1 - COMPLEX_EPSILON - ied->T[j][0] * (alpha_j + beta_j));
  double grad_beta_j = ied->u_beta[j] / beta_j 
    - ied->u_self[j][0] * ied->T[j][0] / 
    (1 - COMPLEX_EPSILON - ied->T[j][0] * (alpha_j + beta_j));
  double grad_tau_j = ied->u_tau[j] / tau_j;

  for (pat = 1; pat < ied->gpm->ngap_patterns - 1; pat++) {
    pattern_type type = gp_pattern_type(ied->gpm, pat);
    if (type == INSERTION_PATTERN) {
      grad_alpha_j -= 
        (ied->u_self[j][pat] * tau_j * ied->T[j][pat] /
         (1 - COMPLEX_EPSILON - tau_j * alpha_j * ied->T[j][pat]
          - tau_j * beta_j * ied->T[j][0] - tau_j));
      grad_beta_j -= 
        (ied->u_self[j][pat] * tau_j * ied->T[j][0] /
         (1 - COMPLEX_EPSILON - tau_j * alpha_j * ied->T[j][pat]
          - tau_j * beta_j * ied->T[j][0] - tau_j));
      grad_tau_j -= 
        (ied->u_self[j][pat] * 
         (alpha_j * ied->T[j][pat] + beta_j * ied->T[j][0] + 1) /
         (1 - COMPLEX_EPSILON - tau_j * alpha_j * ied->T[j][pat]
          - tau_j * beta_j * ied->T[j][0] - tau_j));
    }
    else {
      if (type != DELETION_PATTERN)
	die("ERROR indel_max_gradient, got unknown type %i\n", type);
      grad_alpha_j -= 
        (ied->u_self[j][pat] * tau_j * ied->T[j][0] /
         (1 - COMPLEX_EPSILON - tau_j * alpha_j * ied->T[j][0]
          - tau_j * beta_j * ied->T[j][pat] - tau_j));
      grad_beta_j -= 
        (ied->u_self[j][pat] * tau_j * ied->T[j][pat] /
         (1 - COMPLEX_EPSILON - tau_j * alpha_j * ied->T[j][0]
          - tau_j * beta_j * ied->T[j][pat] - tau_j));
      grad_tau_j -= 
        (ied->u_self[j][pat] * 
         (alpha_j * ied->T[j][0] + beta_j * ied->T[j][pat] + 1) /
         (1 - COMPLEX_EPSILON - tau_j * alpha_j * ied->T[j][0]
          - tau_j * beta_j * ied->T[j][pat] - tau_j));
    }
  }
  
  vec_set(grad, 0, -grad_alpha_j); /* negate each one because function is negated */
  vec_set(grad, 1, -grad_beta_j);
  vec_set(grad, 2, -grad_tau_j);
}

/* Maximize all params for state transition (M step of EM).  This
   function is passed to hmm_train_by_em in phmm_fit_em;  */
void phmm_estim_trans_em(HMM *hmm, void *data, double **A) {

  /* NOTE: if re-estimating state models, be sure to call
     phmm_set_branch_len_factors; it's not called here */

  PhyloHmm *phmm = data;
  int i, j;
  IndelEstimData *ied = NULL;
  double **fcounts;

  if (phmm->em_data->fix_functional && phmm->em_data->fix_indel) return;

  if (phmm->indel_mode == PARAMETERIC) {
    ied = phmm_new_ied(phmm, A);
    fcounts = ied->fcounts;
  }
  else fcounts = A;

  /* estimate transition probs for functional cats in ordinary way but
     using marginal counts */
  if (!phmm->em_data->fix_functional) {
    for (i = 0; i < phmm->functional_hmm->nstates; i++) {
      double rowsum = 0;
      for (j = 0; j < phmm->functional_hmm->nstates; j++)
        rowsum += fcounts[i][j];
      for (j = 0; j < phmm->functional_hmm->nstates; j++)
        mm_set(phmm->functional_hmm->transition_matrix, i, j, 
               fcounts[i][j] / rowsum);
    }
  }

  if (phmm->indel_mode == PARAMETERIC) {
    if (!phmm->em_data->fix_indel) phmm_em_estim_indels(phmm, ied);
    phmm_free_ied(ied);
  }

  phmm_reset(phmm);
}

IndelEstimData *phmm_new_ied(PhyloHmm *phmm, double **A) {
  int i, j;
  IndelEstimData *ied = smalloc(sizeof(IndelEstimData));
  ied->phmm = phmm;

  /* initialize marginal counts for functional categories */
  if (!phmm->em_data->fix_functional) {
    ied->fcounts = smalloc(phmm->functional_hmm->nstates * sizeof(void*));
    for (i = 0; i < phmm->functional_hmm->nstates; i++) {
      ied->fcounts[i] = smalloc(phmm->functional_hmm->nstates * 
                                sizeof(double));
      for (j = 0; j < phmm->functional_hmm->nstates; j++) 
        ied->fcounts[i][j] = 0;
    } 
  }

  /* initialize indel-related marginal counts */
  if (!phmm->em_data->fix_indel) {     /* this part we'll skip if possible */
    ied->nfunctional_states = phmm->functional_hmm->nstates;
    ied->gpm = phmm->gpm;
    ied->u_alpha = smalloc(ied->nfunctional_states * sizeof(double));
    ied->u_beta = smalloc(ied->nfunctional_states * sizeof(double));
    ied->u_tau = smalloc(ied->nfunctional_states * sizeof(double));
    ied->u_self = smalloc(ied->nfunctional_states * sizeof(void*));
    for (i = 0; i < ied->nfunctional_states; i++) {
      ied->u_alpha[i] = ied->u_beta[i] = ied->u_tau[i] = 0;
      ied->u_self[i] = smalloc(ied->gpm->ngap_patterns * sizeof(double));
      for (j = 0; j < ied->gpm->ngap_patterns; j++) ied->u_self[i][j] = 0;
    }
  }

  /* compute marginal counts */
  for (i = 0; i < phmm->hmm->nstates; i++) { 
    int cat_i = phmm->state_to_cat[i];
    int pat_i = phmm->state_to_pattern[i];
    pattern_type pat_i_type = gp_pattern_type(phmm->gpm, pat_i);

    if (cat_i < 0 || cat_i >= phmm->functional_hmm->nstates)
      die("ERROR phmm_new_ied: cat=%i, should be in [0,%i)\n",
	  cat_i, phmm->functional_hmm->nstates);
    if (pat_i < 0 || pat_i >= phmm->gpm->ngap_patterns)
      die("ERROR phmm_new_ied: pat=%i, should bein [0, %i)\n",
	  pat_i, phmm->gpm->ngap_patterns);

    for (j = 0; j < phmm->hmm->nstates; j++) {
      int cat_j = phmm->state_to_cat[j];
      int pat_j = phmm->state_to_pattern[j];
      pattern_type pat_j_type = gp_pattern_type(phmm->gpm, pat_j);

      if (cat_j < 0 || cat_j >= phmm->functional_hmm->nstates)
	die("ERROR phmm_new_ied: cat_j=%i, should be in [0,%i)\n",
	    cat_j, phmm->functional_hmm->nstates);
      if (pat_j < 0 || pat_j >= phmm->gpm->ngap_patterns)
	die("ERROR phmm_new_ied: pat_j=%i, should bein [0, %i)\n",
	    pat_j, phmm->gpm->ngap_patterns);

      if (!phmm->em_data->fix_functional)
        ied->fcounts[cat_i][cat_j] += A[i][j];

      if (phmm->em_data->fix_indel) continue;

      if (pat_i_type == COMPLEX_PATTERN || pat_j_type == COMPLEX_PATTERN) 
        continue;             /* include these for ied->fcounts but not for
                                 the indel-related marginals */

      if (pat_i == pat_j) 
        ied->u_self[cat_j][pat_j] += A[i][j]; /* c_.kjk from notes */
      else {
        if (pat_j_type == INSERTION_PATTERN)
          ied->u_alpha[cat_j] += A[i][j];
        else if (pat_j_type == DELETION_PATTERN)
          ied->u_beta[cat_j] += A[i][j];

        if (pat_i_type != NULL_PATTERN)
          ied->u_tau[cat_j] += A[i][j];
      }
    }
  }

  return ied;
}

void phmm_free_ied(IndelEstimData *ied) {
  int i;
  if (!ied->phmm->em_data->fix_functional) {
    for (i = 0; i < ied->phmm->functional_hmm->nstates; i++) 
      free(ied->fcounts[i]);
    free(ied->fcounts);
  }
  if (!ied->phmm->em_data->fix_indel) {
    for (i = 0; i < ied->phmm->functional_hmm->nstates; i++) 
      free(ied->u_self[i]);
    free(ied->u_alpha);
    free(ied->u_beta);
    free(ied->u_tau);
    free(ied->u_self);
  }
  free(ied);
}

/* for E step of EM: estimate indel params using a multi-dimensional
   optimization routine.  This can be done separately for each
   functional category  */
void phmm_em_estim_indels(PhyloHmm *phmm, IndelEstimData *ied) {
  int i;
  Vector *params = vec_new(3);
  Vector *lb = vec_new(3);
  vec_zero(lb);
  ied->T = phmm->T;

  for (i = 0; i < phmm->functional_hmm->nstates; i++) {
    double retval;
    ied->current_dest_cat = i;

    vec_set(params, 0, phmm->alpha[i]);
    vec_set(params, 1, phmm->beta[i]);
    vec_set(params, 2, phmm->tau[i]);
    opt_bfgs(indel_max_function, params, ied, &retval, lb, NULL, NULL, 
             indel_max_gradient, OPT_HIGH_PREC, NULL); 
    phmm->alpha[i] = vec_get(params, 0);
    phmm->beta[i] = vec_get(params, 1);
    phmm->tau[i] = vec_get(params, 2);
  }
  vec_free(params);
  vec_free(lb);
}

