/* $Id: puzzler.c,v 1.2 2004-06-11 05:58:51 acs Exp $
   Written by Adam Siepel, 2003
   Copyright 2003, Adam Siepel, University of California */

/* Functions to puzzle through the somewhat complex relationships
   among HMM states, tree models, and functional/rate categories.
   Allows for "unspooling" of categories, tied categories, cross
   product of functional and rate categories, "reflection" of HMM
   about pivot states to accommodate reverse strand, and special
   states for indel modeling.  */

/* TO DO: general state tying not implemented yet.  Unspooling
   seems to be sufficient for now */

#include <puzzler.h>
#include <dgamma.h>
#include <sufficient_stats.h>
#include <gap_patterns.h>

/* create a new puzzler.  Inputs are a category map (which
   encapsulates unspooling and category-tying information), a list of
   tree models, an HMM, whether to reflect the HMM, an optional list
   of additional pivot category names for reflection (pass NULL to
   ignore), a number of rate categories (pass 1 for no cross product
   or -1 to use number from tree models), an optional list of explicit
   scaling constants for tree models (if NULL, "auto-rates" will be
   assumed -- see label.c), an optional autocorrelation parameter
   (ignored if nratecats == 1), an optional list of category names for
   which to model indels, and the multiple alignment to be used
   (required if indel_cats != NULL).  The number of tree_models is
   assumed to be equal to the number of (spooled) categories.  The
   Puzzler will assume responsibility for memory management of the
   tree models and HMM but not the category map. */
PhyloHMM_Puzzler *puz_new(CategoryMap *cm, TreeModel **tree_models, HMM *hmm, 
                          int reflect_hmm, List *pivot_cats, int nratecats,
                          double *scaling_consts, double lambda, 
                          List *indel_cats, int nseqs) {
  int max_nstates, s, cat;

  PhyloHMM_Puzzler *puz = smalloc(sizeof(PhyloHMM_Puzzler));
  GapPatternMap *gpm = NULL;

  puz->cm = cm;
  puz->hmm = hmm;
  puz->mods = tree_models;
  puz->nmods = cm->ncats + 1;
  puz->nratecats = nratecats;
  puz->functional_hmm = puz->autocorr_hmm = NULL;

  if (puz->nratecats < 0) {     /* obtain nratecats from tree models;
                                   all must be the same */
    int mod;
    for (mod = 0; mod < puz->nmods; mod++) {
      if (puz->mods[mod]->tree == NULL) continue; /* ignore weight-matrix 
                                                     tree models */
      if (puz->nratecats < 0) puz->nratecats = puz->mods[mod]->nratecats;
      else if (puz->nratecats != puz->mods[mod]->nratecats) {
        fprintf(stderr, "ERROR: tree models specify different numbers of rate categories.\n");
        exit(1);
      }
    }
  }
  assert(puz->nratecats >= 1);

  max_nstates = hmm->nstates * 2 * puz->nratecats;
  puz->state_to_mod = smalloc(max_nstates * sizeof(int));
  puz->state_to_cat = smalloc(max_nstates * sizeof(int));
  puz->state_to_pattern = smalloc(max_nstates * sizeof(int));
  puz->reverse_compl = smalloc(max_nstates * sizeof(int));
  puz->cat_to_states = smalloc((cm->ncats+1) * sizeof(List*));
  for (cat = 0; cat <= cm->ncats; cat++) 
    puz->cat_to_states[cat] = lst_new_int(5);

  /* we need separate initialization logic for the cases with and
     without indel cats.  When there are no indel cats, the HMM states
     correspond to ordinary unspooled categories, but when there are
     indel cats, they correspond to unspooled *gap* categories.
     Without indel cats, there is a direct correspondence between
     spooled categories and models; with indel cats, multiple
     spooled gap categories map to the same model */

  if (indel_cats != NULL)
    /* start by enlarging catmap to allows for gap cats, and creating
       mappings between cats and gapcats */
    gpm = gp_create_gapcats(puz->cm, indel_cats, nseqs);

  /* now the number of unspooled categories should equal the number of
     states in the HMM */
  if (hmm->nstates != (puz->cm->unspooler == NULL ? puz->cm->ncats + 1 : 
                       puz->cm->unspooler->nstates_unspooled)) {
    fprintf(stderr, "ERROR: number of states in HMM must equal number of site categories (unspooled).\n");
    exit(1);
  }

  for (s = 0; s < max_nstates; s++) {
    puz->reverse_compl[s] = 0;
    puz->state_to_pattern[s] = -1;

    if (s >= hmm->nstates) {
      puz->state_to_cat[s] = puz->state_to_mod[s] = -1;
      continue; 
    }

    /* in both cases below, 'cat' is the spooled ordinary (non-gap)
       category corresponding to state s */
    if (indel_cats == NULL) 
      cat = (cm->unspooler == NULL ? s : 
             cm->unspooler->unspooled_to_spooled[s]);
    else {
      int sp_gapcat = puz->cm->unspooler->unspooled_to_spooled[s];
      cat = gpm->gapcat_to_cat[sp_gapcat];
      puz->state_to_pattern[s] = gpm->gapcat_to_pattern[sp_gapcat];
    }

    puz->state_to_cat[s] = puz->state_to_mod[s] = cat;
    lst_push_int(puz->cat_to_states[cat], s); 
                                /* FIXME: I don't think this gets
                                   propagated through puz_reflect_hmm
                                   and puz_rate_cats yet */
  }

  /* FIXME: tied states -- encapsulate in CM?  Maybe not: it's really
     a separate issue.  CM is concerned with mapping between feature
     types and site labels; tied states are specific to the phyloHMM,
     and have more to do with efficiency of computing emissions. */

  if (reflect_hmm) 
    puz_reflect_hmm(puz, pivot_cats);

  if (puz->nratecats > 1) 
    /* scale models and create HMM cross product */
    puz_rate_cats(puz, scaling_consts, lambda);

  if (gpm != NULL) gp_free_map(gpm);

  return puz;
}

/* Reflect HMM about pivot states corresponding to specified list of
   category names; update all mappings accordingly */
void puz_reflect_hmm(PhyloHMM_Puzzler *puz, List *pivot_cats) {
  int ncats = puz->cm->ncats + 1;
/*   int ncats_unspooled = puz->cm->unspooler != NULL ?  */
/*     puz->cm->unspooler->nstates_unspooled : ncats; */
  int mark[ncats];
  int i;
  int *new_to_old;
  HMM *tmp_hmm;
  List *pivot_states = lst_new_int(50);

  for (i = 0; i < ncats; i++) mark[i] = 0;
  mark[0] = 1;                  /* background is implicit */

  if (pivot_cats != NULL) {
    List *cats = cm_get_category_list(puz->cm, pivot_cats, 0);
    for (i = 0; i < lst_size(cats); i++)
      mark[lst_get_int(cats, i)] = 1;
    lst_free(cats);
  }


    /* find out which spooled cat numbers are implied */
/*     for (i = 0; i < lst_size(pivot_cats); i++) { */
/*       int cat = cm_get_category(puz->cm, lst_get_ptr(pivot_cats, i)); */
/*       for (j = puz->cm->ranges[cat]->start_cat_no;  */
/*            j <= puz->cm->ranges[cat]->start_cat_no; j++) { */
/*         mark[cat] = 1; */
/*         if (puz->cm->unspooler == NULL) */ /* no unspooler; can add immediately */
/*           lst_push_int(pivot_states, cat); */
/*       } */
/*     } */
/*   } */

    /* if using an indel model, we need to consider all gap
       categories derived from the specified categories */
/*     if (gpm != NULL) { */
/*       for (i = 0; i < ncats; i++) { */
/*         if (!mark[i] && mark[gpm->gapcat_to_cat[i]]) { */
/*           mark[i] = 1; */
/*           if (puz->cm->unspooler == NULL)  */
/*             lst_push_int(pivot_states, i); */
/*         } */
/*       } */
/*     } */ /* FIXME: what about when pivot_cats == NULL? */

/*     if (puz->cm->unspooler != NULL) { */
/*       for (i = 0; i < ncats_unspooled; i++) */
/*         if (mark[puz->cm->unspooler->unspooled_to_spooled[i]]) */
/*           lst_push_int(pivot_states, i); */
/*     } */
/*   } */

  /* enumerate states corresponding to pivot categories.  This will
     take care of indel categories as well as ordinary unspooled
     categories */
  for (i = 1; i < puz->hmm->nstates; i++) /* skip 0; already taken care of */
    if (mark[puz->state_to_cat[i]])
      lst_push_int(pivot_states, i);

  /* now reflect HMM */
  new_to_old = smalloc(puz->hmm->nstates * 2 * sizeof(int));
  tmp_hmm = hmm_reverse_compl(puz->hmm, pivot_states, new_to_old); 
  hmm_free(puz->hmm);
  puz->hmm = tmp_hmm;

  /* finally, revise mappings */
  for (i = puz->hmm->nstates-1; i >= 0; i--) {
    /* we use the fact that abs(new_to_old[i]) <= i to do the
       remapping in place */
    puz->state_to_mod[i] = puz->state_to_mod[abs(new_to_old[i])];
    puz->state_to_cat[i] = puz->state_to_cat[abs(new_to_old[i])];
    puz->state_to_pattern[i] = puz->state_to_pattern[abs(new_to_old[i])];
    puz->reverse_compl[i] = (new_to_old[i] < 0);
  }

  free(new_to_old);
  lst_free(pivot_states);
}

/* create an HMM representing the autocorrelation model of Felsenstein
   and Churchill; HMM object must already be allocated with desired
   number of states */
void puz_create_autocorr_hmm(HMM *hmm, double lambda) {
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

/* Scale tree models and replace HMM by its cross-product with an
   autocorr HMM.  Update all mappings appropriately.  Weight-matrix
   tree models are not scaled.  If scaling_consts is non-NULL, it is
   assumed to be an array of length puz->nratecats containing explicit
   scaling constants; otherwise, "auto scaling" is done, using the
   discrete gamma model (see label.c) */
void puz_rate_cats(PhyloHMM_Puzzler *puz, double *scaling_consts,
                   double lambda) {

  TreeModel **new_mods;
  double *sconsts;
  int mod, state, ratecat, thismod_new = 0;
  int *old_mod_to_new_base;     /* mapping from old model index to
                                   base index in new space */

  new_mods = smalloc(sizeof(TreeModel*) * puz->nmods * puz->nratecats);  
  sconsts = (scaling_consts != NULL ? scaling_consts : 
             smalloc(puz->nratecats * sizeof(double)));
  old_mod_to_new_base = smalloc(puz->nmods * sizeof(int));

  for (mod = 0; mod < puz->nmods; mod++) {
    if (puz->mods[mod]->tree == NULL) { /* weight matrix */
      old_mod_to_new_base[mod] = thismod_new;
      new_mods[thismod_new++] = puz->mods[mod]; /* just reuse same object */
    }
    else {      /* full tree model */
      if (scaling_consts == NULL) { 
        /* in this case, the scaling constants must be determined from
           the model (and may be different for each model) */
        double freqs[puz->nratecats];
        assert(puz->mods[mod]->alpha >= 0); /* FIXME: check in label.c */
        DiscreteGamma(freqs, sconsts, puz->mods[mod]->alpha, 
                      puz->mods[mod]->alpha, puz->nratecats, 0);
      }
          
      tm_reinit(puz->mods[mod], puz->mods[mod]->subst_mod, 1, 0, NULL, NULL);
      /* we don't want to treat the scaled versions as
         discrete-gamma model (would be "counting twice") */

      /* now make nratecats scaled copies of model */
      old_mod_to_new_base[mod] = thismod_new;
      for (ratecat = 0; ratecat < puz->nratecats; ratecat++) {
        new_mods[thismod_new] = tm_create_copy(puz->mods[mod]);

        fprintf(stderr, "Scaling tree model %d by factor %f ...\n",
                mod+1, sconsts[ratecat]); /* remove? */

        tm_scale(new_mods[thismod_new], sconsts[ratecat], 1);
        thismod_new++;
      }
    }
  }

  /* revise mappings */
  for (state = puz->hmm->nstates-1; state >= 0; state--) { 
    for (ratecat = puz->nratecats-1; ratecat >= 0; ratecat--) {
      int new_state = state*puz->nratecats + ratecat;
      puz->state_to_mod[new_state] = 
        old_mod_to_new_base[puz->state_to_mod[state]] + 
        (puz->mods[puz->state_to_mod[state]]->tree == NULL ? 0 : ratecat); 
      puz->state_to_cat[new_state] = puz->state_to_cat[state];
      puz->state_to_pattern[new_state] = puz->state_to_pattern[state];
      puz->reverse_compl[new_state] = puz->reverse_compl[state];
      /* we use the fact that new_state > old_state
         except when state == 0 && ratecat == 0 (the last case
         considered) to do the remapping in place  */
    }
  }

  /* replace HMM with cross product of (functional) HMM and autocorr HMM */
  puz->functional_hmm = puz->hmm;
  puz->autocorr_hmm = hmm_new_nstates(puz->nratecats, 1, 0);
  puz_create_autocorr_hmm(puz->autocorr_hmm, lambda);
  puz->hmm = hmm_new_nstates(puz->nratecats * puz->functional_hmm->nstates, 1, 0);
  hmm_cross_product(puz->hmm, puz->functional_hmm, puz->autocorr_hmm);

  for (mod = 0; mod < puz->nmods; mod++)
    if (puz->mods[mod]->tree != NULL) tm_free(puz->mods[mod]);
  free(puz->mods);
  puz->mods = new_mods;
  puz->nmods = thismod_new;
  free(old_mod_to_new_base);
  if (scaling_consts == NULL) free(sconsts);
}

/* update cross product HMM using new value of lambda */
void puz_update_cross_prod(PhyloHMM_Puzzler *puz, double lambda) {
  puz_create_autocorr_hmm(puz->autocorr_hmm, lambda);
  hmm_cross_product(puz->hmm, puz->functional_hmm, puz->autocorr_hmm);
}

void puz_free(PhyloHMM_Puzzler *puz) {
  int i;
  for (i = 0; i < puz->nmods; i++) tm_free(puz->mods[i]);
  free(puz->mods);
  free(puz->state_to_mod);
  free(puz->state_to_cat);
  free(puz->reverse_compl);
  hmm_free(puz->hmm);
  if (puz->functional_hmm != NULL) hmm_free(puz->functional_hmm);
  if (puz->autocorr_hmm != NULL) hmm_free(puz->autocorr_hmm);
  free(puz);
}

