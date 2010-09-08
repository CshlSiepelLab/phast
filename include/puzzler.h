/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#ifndef PUZ_H
#define PUZ_H

#include <category_map.h>
#include <hmm.h>
#include <tree_model.h>
#include <lists.h>

typedef struct {
  CategoryMap *cm;              /* category map */
  HMM *hmm;                     /* final HMM, after reflection and
                                   allowance for rate categories  */
  HMM *functional_hmm, *autocorr_hmm;
                                /* original HMMs used to create cross
                                   product (NULL if no rate
                                   categories) */
  TreeModel **mods;             /* array of tree models, after
                                   allowance for rate categories  */
  int nmods;                    /* number of tree models (length of
                                   array mods) */
  int *state_to_mod;            /* mapping of HMM state number to tree
                                   model number */
  int *state_to_cat;            /* mapping of HMM state number to (spooled)
                                   category number */
  int *reverse_compl;           /* array of length hmm->nstates
                                   with value 1 for each state
                                   that corresponds to the reverse
                                   strand and value 0 otherwise */
  List **cat_to_states;         /* one to many mapping */
  int *state_to_pattern;        /* gap pattern associated with each
                                   state, when modeling indels (-1 for
                                   no gap pattern) */
  int nratecats;
} PhyloHMM_Puzzler;

PhyloHMM_Puzzler *puz_new(CategoryMap *cm, TreeModel **tree_models, HMM *hmm, 
                          int reflect_hmm, List *pivot_cats, int nratecats,
                          double *scaling_consts, double lambda, 
                          List *indel_cats, int nseqs);
void puz_reflect_hmm(PhyloHMM_Puzzler *puz, List *pivot_cats);
void puz_create_autocorr_hmm(HMM *hmm, double lambda);
void puz_rate_cats(PhyloHMM_Puzzler *puz, double *scaling_consts,
                   double lambda);
void puz_update_cross_prod(PhyloHMM_Puzzler *puz, double lambda);
void puz_free(PhyloHMM_Puzzler *puz);

#endif
