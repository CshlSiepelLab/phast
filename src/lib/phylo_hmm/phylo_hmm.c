/* $Id: phylo_hmm.c,v 1.10 2004-08-05 07:15:04 acs Exp $
   Written by Adam Siepel, 2003
   Copyright 2003, Adam Siepel, University of California */

/* Code for phylo-HMMs.  Allows for automatic expansion of the state
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

/** Create a new PhyloHmm object. Optionally expands original HMM to
    allow for features on both the positive and negative strands.
    Also optionally expands original HMM for indel modeling. */
PhyloHmm *phmm_new(HMM *hmm,    /**< HMM.  A copy is created,
                                   orig. isn't touched.  If NULL, a
                                   trivial, single-state HMM will be
                                   created.  */
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
                                   indel_cats != NULL. */
                   List *pivot_cats,  
                                /**< Categories (by name or number)
                                   about which to "reflect" the HMM.
                                   Allows for prediction on both
                                   strands (see hmm_reverse_compl).
                                   Pass NULL for no reflection */
                   int indels, 
                                /**< Whether to use indel model. */
                   int nseqs    /**< Number of sequences to be used.
                                   Required for indel modeling, ignored
                                   if indels == FALSE. */
                   ) { 
  int max_nstates, s, cat, j;

  PhyloHmm *phmm = smalloc(sizeof(PhyloHmm));
  GapPatternMap *gpm = NULL;

  if (hmm != NULL)
    phmm->hmm = hmm_create_copy(hmm);
  else 
    phmm->hmm = hmm_create_trivial();

  hmm = NULL;                   /* be sure this doesn't get used by
                                   mistake */

  if (cm != NULL)
    phmm->cm = cm_create_copy(cm);
  else {
    if (indels) 
      die("ERROR: must pass non-NULL category map if using indel model");
                                /* would get a little tricky without a cm */
    phmm->cm = cm_create_trivial(phmm->hmm->nstates-1, "model_");
  }

  cm = NULL;                    /* be sure doesn't get used by mistake */

  phmm->mods = tree_models;
  phmm->nmods = phmm->cm->ncats + 1;
  phmm->nratecats = 1;
  phmm->functional_hmm = phmm->autocorr_hmm = NULL;
  phmm->reflected = pivot_cats != NULL;
  phmm->emissions = NULL;
  phmm->forward = NULL;
  phmm->alloc_len = -1;
  phmm->state_pos = phmm->state_neg = NULL;
  phmm->do_indels = indels;
  phmm->topology = NULL;

  /* make sure tree models all have trees and all have the same number
     of leaves; keep a pointer to a representative tree for use with
     indel model (in this case, topologies must be the same) */
  for (s = 0; s < phmm->nmods; s++) {
    if (phmm->mods[s]->tree == NULL) 
      die("ERROR: tree model #%d for phylo-HMM has no tree.\n", s+1);
    else if (phmm->topology == NULL)
      phmm->topology = phmm->mods[s]->tree;
    else if (phmm->mods[s]->tree->nnodes != phmm->topology->nnodes) 
      die("ERROR: tree models for phylo-HMM have different numbers of nodes.\n");
  }

  /* now initialize mappings */
  max_nstates = phmm->reflected ? phmm->hmm->nstates * 2 : phmm->hmm->nstates;
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

  if (indels) {
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
    gpm = gp_create_gapcats(phmm->cm, indel_types, nseqs);
    lst_free(indel_types);
  }

  /* now the number of unspooled categories should equal the number of
     states in the HMM */
  if (phmm->hmm->nstates != 
      (phmm->cm->unspooler == NULL ? phmm->cm->ncats + 1 : 
       phmm->cm->unspooler->nstates_unspooled)) 
    die("ERROR: number of states in HMM must equal number of site categories (unspooled).\n");

  for (s = 0; s < max_nstates; s++) {
    phmm->reverse_compl[s] = 0;
    phmm->state_to_pattern[s] = -1;

    if (s >= phmm->hmm->nstates) {
      phmm->state_to_cat[s] = phmm->state_to_mod[s] = -1;
      continue; 
    }

    /* in both cases below, 'cat' is the spooled ordinary (non-gap)
       category corresponding to state s */
    if (!indels) 
      cat = (phmm->cm->unspooler == NULL ? s : 
             phmm->cm->unspooler->unspooled_to_spooled[s]);
    else {
      int sp_gapcat = phmm->cm->unspooler->unspooled_to_spooled[s];
      cat = gpm->gapcat_to_cat[sp_gapcat];
      phmm->state_to_pattern[s] = gpm->gapcat_to_pattern[sp_gapcat];
    }

    phmm->state_to_cat[s] = phmm->state_to_mod[s] = cat;
    lst_push_int(phmm->cat_to_states[cat], s); 
                                /* FIXME: I don't think this gets
                                   propagated through phmm_reflect_hmm
                                   and phmm_rate_cats yet */
  }

  if (pivot_cats != NULL) 
    phmm_reflect_hmm(phmm, pivot_cats);

  if (gpm != NULL) gp_free_map(gpm);

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
    for (i = 0; i < lst_size(cats); i++)
      mark[lst_get_int(cats, i)] = 1;
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
  gsl_vector_set_all(hmm->eq_freqs, (double)1/hmm->nstates);
  gsl_vector_set_all(hmm->begin_transitions, (double)1/hmm->nstates);

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
      assert(phmm->mods[mod]->alpha > 0); 
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

        tm_scale(new_mods[thismod_new], sconsts[ratecat], 1);
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
  phmm->functional_hmm = phmm->hmm;
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
  hmm_free(phmm->hmm);
  cm_free(phmm->cm);
  if (phmm->functional_hmm != NULL) hmm_free(phmm->functional_hmm);
  if (phmm->autocorr_hmm != NULL) hmm_free(phmm->autocorr_hmm);

  free(phmm);
}

/** Compute emissions for given PhyloHmm and MSA.  Preprocessor for
    phmm_viterbi_features, phmm_posterior_probs, and phmm_lnl
    (typically only needs to be run once). */
void phmm_compute_emissions(PhyloHmm *phmm,
                                /**< Initialized PhyloHmm */
                            MSA *msa,
                                /**< Source alignment */
                            int quiet
                                /**< Whether to report progress to stderr */
                            ) {

  int i, mod, j;
  MSA *msa_compl = NULL;

  phmm->emissions = smalloc(phmm->hmm->nstates * sizeof(double*));  
  phmm->alloc_len = msa->length;

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
    for (i = 0; i < msa_compl->nseqs; i++) free(msa_compl->seqs[i]);
    free(msa_compl->seqs);
    msa_compl->seqs = NULL;
  }

  /* set up mapping from model/strand to first associated state
     (allows phmm->emissions to be computed only once for each
     model/strand pair) */
  phmm->state_pos = smalloc(phmm->nmods * sizeof(int));
  phmm->state_neg = smalloc(phmm->nmods * sizeof(int));
  for (i = 0; i < phmm->nmods; i++) 
    phmm->state_pos[i] = phmm->state_neg[i] = -1;

  for (i = 0; i < phmm->hmm->nstates; i++) {
    if (!quiet) {
      fprintf(stderr, "Computing emission probs (state %d, cat %d, mod %d",
              i+1, phmm->state_to_cat[i]+1, phmm->state_to_mod[i]+1);
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
      phmm->emissions[i] = smalloc(msa->length * sizeof(double));
      tl_compute_log_likelihood(phmm->mods[mod], 
                                phmm->reverse_compl[i] ? msa_compl : msa,
                                phmm->emissions[i], -1, NULL);
      if (!phmm->reverse_compl[i]) phmm->state_pos[mod] = i;
      else phmm->state_neg[mod] = i;            
    }
  }

  /* finally, adjust for indel model, if necessary */
  if (phmm->do_indels) {
    int *msa_gap_patterns = smalloc(msa->length * sizeof(int));

    if (!quiet)
      fprintf(stderr, "Obtaining gap patterns...\n");

    gp_set_phylo_patterns(msa_gap_patterns, msa, phmm->topology);

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
        for (j = 0; j < msa->length; j++) 
          phmm->emissions[i][j] = 
            (msa_gap_patterns[j] == phmm->state_to_pattern[i] ? 
             orig_emissions[j] : NEGINFTY);                                 
      }
    }
    free(msa_gap_patterns);
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

  /* filter out unwanted types */
  gff_filter_by_type(retval, types, FALSE, NULL);

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

/* assumes 1-state HMM, 1-cat category map, and single tree model.
   Currently not compatible with 'reflected' HMMs or indel model.
   Assumes category map is to be updated */
void phmm_rates_cut(PhyloHmm *phmm, 
                    int nrates, /* may be different from value in tree
                                   model if dgamma */
                    int cut_idx, 
                    double p, 
                    double q) {
  double freq1 = 0;
  int i;
  int dgamma = !phmm->mods[0]->empirical_rates; /* whether using
                                                   discrete gamma
                                                   model */

  String *newtype;
  List *rconsts, *rweights;
  double tmp_rates[nrates], tmp_weights[nrates];

  assert(phmm->hmm != NULL && phmm->hmm->nstates == 1 && phmm->nmods == 1);
  assert(phmm->cm->ncats == 0);
  assert(nrates > 1);
  assert(cut_idx >= 1 && cut_idx <= nrates);

  hmm_free(phmm->hmm);
  phmm->hmm = hmm_new_nstates(2, TRUE, FALSE);

  if (dgamma) 
    /* if using dgamma, need to compute rate consts and weights -- may
       not have been computed yet */
    DiscreteGamma(tmp_weights, tmp_rates, phmm->mods[0]->alpha, 
                  phmm->mods[0]->alpha, nrates, 0); 
  else {
    for (i = 0; i < nrates; i++) {
      tmp_weights[i] = phmm->mods[0]->freqK[i];
      tmp_rates[i] = phmm->mods[0]->rK[i];
    }
  }

  /* set HMM transitions according to p and q */
  mm_set(phmm->hmm->transition_matrix, 0, 0, 1-p);
  mm_set(phmm->hmm->transition_matrix, 0, 1, p);
  mm_set(phmm->hmm->transition_matrix, 1, 0, q);
  mm_set(phmm->hmm->transition_matrix, 1, 1, 1-q);

  /* set HMM begin transitions according to weights */
  for (i = 0; i < cut_idx; i++) freq1 += tmp_weights[i];
  gsl_vector_set(phmm->hmm->begin_transitions, 0, freq1);
  gsl_vector_set(phmm->hmm->begin_transitions, 1, 1 - freq1);

  hmm_reset(phmm->hmm);

  /* create 2nd tree model, then update rate categories in both */
  phmm->mods = srealloc(phmm->mods, 2 * sizeof(void*));
  phmm->nmods = 2;
  phmm->mods[1] = tm_create_copy(phmm->mods[0]);

  rconsts = lst_new_dbl(nrates);
  rweights = lst_new_dbl(nrates);
  for (i = 0; i < cut_idx; i++) {
    lst_push_dbl(rweights, tmp_weights[i]);
    lst_push_dbl(rconsts, tmp_rates[i]);
  }

  if (cut_idx == 1) {
    tm_reinit(phmm->mods[0], phmm->mods[0]->subst_mod, 1, 0, NULL, NULL);
    tm_scale(phmm->mods[0], lst_get_dbl(rconsts, 0), TRUE);
                                /* in this case, have to by-pass rate
                                   variation machinery and just scale
                                   tree; tree model code won't do
                                   rate variation with single rate
                                   category */
  }
  else
    tm_reinit(phmm->mods[0], phmm->mods[0]->subst_mod, cut_idx, 
              phmm->mods[0]->alpha, rconsts, rweights);
                                /* note that dgamma model will be
                                   redefined as empirical rates mod */

  lst_clear(rconsts); lst_clear(rweights);
  for (i = cut_idx; i < nrates; i++) {
    lst_push_dbl(rweights, tmp_weights[i]);
    lst_push_dbl(rconsts, tmp_rates[i]);
  }

  if (cut_idx == nrates-1) {    /* unlikely but possible */
    tm_reinit(phmm->mods[1], phmm->mods[1]->subst_mod, 1, 0, NULL, NULL);
    tm_scale(phmm->mods[1], lst_get_dbl(rconsts, nrates-1), TRUE);
  }
  else
    tm_reinit(phmm->mods[1], phmm->mods[1]->subst_mod, nrates - cut_idx, 
              phmm->mods[1]->subst_mod, rconsts, rweights);
  lst_free(rconsts); lst_free(rweights);

  /* expand category map */
  cm_realloc(phmm->cm, 1);
  str_cpy_charstr(lst_get_ptr(phmm->cm->ranges[0]->feature_types, 0), 
                  "conserved");
  newtype = str_new_charstr("nonconserved");
  phmm->cm->ranges[1] = cm_new_category_range(newtype, 1, 1);
  assert(phmm->cm->conditioned_on[0] == NULL);
  
  /* update mappings */
  phmm->state_to_mod = srealloc(phmm->state_to_mod, 2 * sizeof(int));
  phmm->state_to_cat = srealloc(phmm->state_to_cat, 2 * sizeof(int));
  phmm->state_to_pattern = srealloc(phmm->state_to_pattern, 2 * sizeof(int));
  phmm->reverse_compl = srealloc(phmm->reverse_compl, 2 * sizeof(int));
  phmm->state_to_mod[1] = 1; 
  phmm->state_to_cat[1] = 1;
  phmm->state_to_pattern[1] = -1;
  phmm->reverse_compl[1] = 0;
}

/* function used by phmm_fit_rates_cut */
void compute_emissions(double **emissions, void **models, int nmodels,
                       void *data, int sample, int length) {
  /* just copy emissions; they're already computed */
  PhyloHmm *phmm = (PhyloHmm*)data;
  int i, j;
  for (i = 0; i < phmm->hmm->nstates; i++)
    for (j = 0; j < phmm->alloc_len; j++)
      emissions[i][j] = phmm->emissions[i][j];
}

/** Estimate the parameters 'p' and 'q' that define the two-state
    "rates-cut" model using an EM algorithm.  Returns ln likelihood. */
double phmm_fit_rates_cut(PhyloHmm *phmm, 
                          double *p, 
                          double *q, 
                          FILE *logf
                          ) {
  double retval;

  if (phmm->emissions == NULL)
    die("ERROR: emissions required for phmm_fit_rates_cut.\n");

  mm_set(phmm->hmm->transition_matrix, 0, 0, 1-*p);
  mm_set(phmm->hmm->transition_matrix, 0, 1, *p);
  mm_set(phmm->hmm->transition_matrix, 1, 0, *q);
  mm_set(phmm->hmm->transition_matrix, 1, 1, 1-*q);
  hmm_reset(phmm->hmm);

  retval = hmm_train_by_em(phmm->hmm, phmm->mods, phmm, 1, &phmm->alloc_len, 
                           NULL, compute_emissions, NULL, NULL, logf);

  *p = mm_get(phmm->hmm->transition_matrix, 0, 1);
  *q = mm_get(phmm->hmm->transition_matrix, 1, 0);

  return log(2) * retval;
}
 
/* wrapper function used by fit_rates_cut_bfgs */
double fit_rates_cut_lnl(gsl_vector *params, void *data) {
  PhyloHmm *phmm = data;
  double p = exp(gsl_vector_get(params, 0)), 
    q = exp(gsl_vector_get(params, 1));
  mm_set(phmm->hmm->transition_matrix, 0, 0, 1-p);
  mm_set(phmm->hmm->transition_matrix, 0, 1, p);
  mm_set(phmm->hmm->transition_matrix, 1, 0, q);
  mm_set(phmm->hmm->transition_matrix, 1, 1, 1-q);
  hmm_reset(phmm->hmm);
  return log(2) * -hmm_forward(phmm->hmm, phmm->emissions, 
                               phmm->alloc_len, phmm->forward);
}

/* returns ln likelihood */
double phmm_fit_rates_cut_bfgs(PhyloHmm *phmm, double *p, double *q, 
                               FILE *logf) {
  int i;
  double neglnl = INFTY;
  gsl_vector *params = gsl_vector_alloc(2),  *lbounds = NULL, 
    *ubounds = gsl_vector_calloc(2);

  /* initialize to given values */
  gsl_vector_set(params, 0, log(*p)); /* use log parameterization */
  gsl_vector_set(params, 1, log(*q));

  /* allocate memory for forward alg */
  if (phmm->forward == NULL) {  /* otherwise assume already alloc */
    phmm->forward = smalloc(phmm->hmm->nstates * sizeof(double*));
    for (i = 0; i < phmm->hmm->nstates; i++)
      phmm->forward[i] = smalloc(phmm->alloc_len * sizeof(double));
  }

  opt_bfgs(fit_rates_cut_lnl, params, phmm, &neglnl, lbounds, ubounds, 
           logf, NULL, OPT_HIGH_PREC, NULL);

  *p = exp(gsl_vector_get(params, 0));
  *q = exp(gsl_vector_get(params, 1));

  gsl_vector_free(params);
  gsl_vector_free(ubounds);

  return -neglnl;
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
