/* $Id: phylo_hmm.c,v 1.1 2004-06-09 17:12:15 acs Exp $
   Written by Adam Siepel, 2003
   Copyright 2003, Adam Siepel, University of California */

/* Code for phylo-HMMs.  Allows for automatic expansion of the state
   space to accommodate features on the reverse strand, and for the
   indel model described in Siepel & Haussler (RECOMB '04).  Also
   allows for cross-product constructions involving functional states
   and rate categories (Siepel & Haussler, RECOMB '03).  */

#include <phylo_hmm.h>
#include <dgamma.h>
#include <sufficient_stats.h>
#include <gap_patterns.h>
#include <tree_likelihoods.h>

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
                   List *indel_cats, 
                                /**< Categories (by name or number)
                                   for which to model indels.  Pass
                                   NULL for no indel model.  */
                   int nseqs    /**< Number of sequences to be used.
                                   Requied for indel modeling, ignored
                                   if indel_cats == NULL. */
                   ) { 
  int max_nstates, s, cat;

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
    if (indel_cats != NULL) 
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
  phmm->emissions_len = -1;
  phmm->state_pos = phmm->state_neg = NULL;

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

  if (indel_cats != NULL)
    /* start by enlarging catmap to allows for gap cats, and creating
       mappings between cats and gapcats */
    gpm = gp_create_gapcats(phmm->cm, indel_cats, nseqs);

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
    if (indel_cats == NULL) 
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
    List *cats = cm_get_category_list(phmm->cm, pivot_cats, 0);
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
                                /**< Number of rate categories.  Use
                                   -1 to obtain the number from the
                                   individual tree models in the
                                   PhyloHmm (all must be the same) */
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

  phmm->nratecats = nratecats;

  if (nratecats < 0) {          /* obtain nratecats from tree models;
                                   all must be the same */
    for (mod = 0; mod < phmm->nmods; mod++) {
      if (phmm->mods[mod]->tree == NULL) continue; /* ignore weight-matrix 
                                                      tree models */
      if (phmm->nratecats < 0) phmm->nratecats = phmm->mods[mod]->nratecats;
      else if (phmm->nratecats != phmm->mods[mod]->nratecats) 
        die("ERROR: tree models specify different numbers of rate categories.\n");
    }
  }

  if (phmm->nratecats <= 1) 
    die("ERROR: phmm_rate_cats requires nratecats > 1.\n");

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
      assert(phmm->mods[mod]->alpha >= 0); /* FIXME: check in label.c */
      DiscreteGamma(freqs, sconsts, phmm->mods[mod]->alpha, 
                    phmm->mods[mod]->alpha, phmm->nratecats, 0);
          
      tm_reinit(phmm->mods[mod], phmm->mods[mod]->subst_mod, 1, 0, NULL);
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
        phmm->cm->feat_ext_lst[newcat] = lst_new_ptr(1);
        lst_push_ptr(phmm->cm->feat_ext_lst[newcat], newtype);
        for (i = 0; i < range_size; i++) {
          phmm->cm->conditioned_on[newcat + i] = NULL;
          if (i > 0) {
            phmm->cm->ranges[newcat + i] = phmm->cm->ranges[newcat];
            phmm->cm->feat_ext_lst[newcat + i] = phmm->cm->feat_ext_lst[newcat];

          }
        }
        /* skip labelling precedence, fill precedence, feat_ext_list */
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

  int i, mod;
  MSA *msa_compl = NULL;

  phmm->emissions = smalloc(phmm->hmm->nstates * sizeof(double*));  
  phmm->emissions_len = msa->length;

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
    if (!quiet)
      fprintf(stderr, "Computing emission probs (state %d,  cat %d, mod %d, pattern %d, strand %c) ...\n", 
              i, phmm->state_to_cat[i], phmm->state_to_mod[i], 
              phmm->state_to_pattern[i], phmm->reverse_compl[i] ? '-' : '+');
    
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
}

/** Run the Viterbi algorithm and return a set of predictions.
    Emissions must have already been computed (see
    phmm_compute_emissions) */
GFF_Set* phmm_predict_viterbi(PhyloHmm *phmm, 
                                /**< PhyloHmm object */
                               char *seqname,
                                /**< seqname for feature set (e.g.,
                                   "chr1") */
                               List *frame 
                                /**< features for which to obtain
                                   frame (NULL to ignore) */
                               ) {
  int *path = (int*)smalloc(phmm->emissions_len * sizeof(int));
  GFF_Set *retval;

  if (phmm->emissions == NULL)
    die("ERROR: emissions required for phmm_viterbi_features.\n");
          
  hmm_viterbi(phmm->hmm, phmm->emissions, phmm->emissions_len, path);

  retval = cm_labeling_as_gff(phmm->cm, path, phmm->emissions_len, 
                              phmm->state_to_cat, 
                              phmm->reverse_compl, seqname, "PHAST",  
                              phmm->reflected ? '-' : '+', frame, 
                              NULL);

                                /* FIXME: gff groups (group_root)
                                   should pass tag instead */
  free(path);

  return retval;
}

GFF_Set* phmm_predict_viterbi_cats(PhyloHmm *phmm, 
                                     /**< PhyloHmm object */
                                   List *cats,
                                     /**< states of interest
                                        (integers >= 0 and <
                                        phmm->hmm->nstates) */
                                   char *seqname,
                                     /**< seqname for feature set (e.g.,
                                        "chr1") */
                                   List *frame,
                                     /**< features for which to obtain
                                        frame (NULL to ignore) */
                                   char *new_type
                                     /**< replace type of each
                                        retained feature with this
                                        string if non-NULL */
                                   ) {
  int i, cat, lastone = -1;
  GFF_Set *retval = phmm_predict_viterbi(phmm, seqname, frame);
  List *types = lst_new_ptr(lst_size(cats));
  GFF_Feature *lastkeeper = NULL;
  List *keepers, *catnos;

  /* do this way to allow input to be numbers or names */
  catnos = cm_get_category_list(phmm->cm, cats, 1);
  for (i = 0; i < lst_size(catnos); i++) {
    cat = lst_get_int(catnos, i);
    if (phmm->cm->ranges[cat]->start_cat_no == lastone) continue;
    lst_push_ptr(types, cm_get_feature(phmm->cm, cat));
    lastone = phmm->cm->ranges[cat]->start_cat_no;
  }
  lst_free(catnos);

  /* filter out unwanted types */
  gff_filter_by_type(retval, types);

  /* now merge adjacent features */
  keepers = lst_new_ptr(lst_size(retval->features));
  for (i = 0; i < lst_size(retval->features); i++) {
    GFF_Feature *f = lst_get_ptr(retval->features, i);
    if (f->start == lastkeeper->end+1) {
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
    forward[i] = (double*)smalloc(phmm->emissions_len * sizeof(double));
  logl = hmm_forward(phmm->hmm, phmm->emissions, 
                     phmm->emissions_len, forward);
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

  return hmm_posterior_probs(phmm->hmm, phmm->emissions, phmm->emissions_len, 
                             post_probs) * log(2);
                                /* convert to natural log */          
}

/** Computes and returns an array of length phmm->emissions_len
    representing the marginal posterior prob at every site, summed
    over all states in the specified list.  Emissions must have
    already been computed (see phmm_compute_emissions) Computes log
    likelihood as a side effect. */
double* phmm_postprobs_states(PhyloHmm *phmm, 
                                /**< PhyloHmm object */
                              List *states, 
                                /**< List of states to consider
                                   (integers >= 0 and <
                                   phmm->hmm->nstates) */
                              double *lnl
                                /**< if non-NULL, will point to log
                                   likelihood on return */
                              ) {
  int i, j, state;
  double **pp = smalloc(phmm->hmm->nstates * sizeof(double*));
  double *retval = smalloc(phmm->emissions_len * sizeof(double));
  double l;

  /* only allocate memory for states of interest; NULLs for the others
     will cause hmm_postprobs to ignore them */
  for (i = 0; i < phmm->hmm->nstates; i++) pp[i] = NULL;
  for (i = 0; i < lst_size(states); i++) {
    state = lst_get_int(states, i);
    if (pp[state] == NULL)
      pp[state] = smalloc(phmm->emissions_len * sizeof(double));
  }

  l = phmm_postprobs(phmm, pp);

  for (j = 0; j < phmm->emissions_len; j++) retval[j] = 0;
  for (i = 0; i < lst_size(states); i++) {
    state = lst_get_int(states, i);
    for (j = 0; j < phmm->emissions_len; j++) retval[j] += pp[state][j];
  }
      
  if (lnl != NULL) *lnl = l;

  return retval;
}

typedef struct {                /* "package" for data needed by
                                   log_likelihood_wrapper (see below) */
  double **forward;
  int msa_len;
  PhyloHmm *phmm;
} AutoratesData;

double log_likelihood_wrapper(double lambda, void *data) {
  AutoratesData *ad = (AutoratesData*)data;
  if (lambda < 0 || lambda > 1) return INFTY;
  phmm_update_cross_prod(ad->phmm, lambda);
  return -hmm_forward(ad->phmm->hmm, ad->phmm->emissions, ad->msa_len, 
                      ad->forward);
}

double phmm_fit_lambda(PhyloHmm *phmm, int msa_len) {
  AutoratesData ad;
  double lambda, final_score, ax, bx, cx, fa, fb, fc;
  int i;

  /* allocate memory for forward alg */
  ad.forward = smalloc(phmm->hmm->nstates * sizeof(double*));
  for (i = 0; i < phmm->hmm->nstates; i++)
    ad.forward[i] = smalloc(msa_len * sizeof(double));

  ad.phmm = phmm;
  ad.msa_len = msa_len;

  ax = .80; bx = .97;           /* FIXME -- parameterize */
  mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, log_likelihood_wrapper, &ad);
  final_score = opt_brent(ax, bx, cx, log_likelihood_wrapper, 5e-3, &lambda, &ad);
/*   fprintf(stderr, "Returned from opt_brent; freeing forward ...\n"); */

  for (i = 0; i < phmm->hmm->nstates; i++) free(ad.forward[i]);
  free(ad.forward);

  return lambda;
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
                            List *null_cats
                                /**< Categories in null model.  Pass
                                   NULL to use all categories not in
                                   score_cats or helper_cats */
                            )  {

  int i, j, cat, ncats, nscore, state;
  List *cats, *score_states, *null_states;
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
    helper_cats != NULL ? lst_size(helper_cats) : 0;
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
    if (str_in_list(feat->feature, score_cats)) {
      int start = feat->start;
      int end = feat->end;

      /* extend range as far as possible using helper cats */
      if (helper_cats != NULL) {
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
  for (i = 0; i < ncats; i++) lst_free(cat_to_states[i]);
  free(cat_to_states);
  free(is_scored);
}
