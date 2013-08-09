/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_phastCons.c
The RPHAST handles to the phastCons program

Melissa Hubisz
Last updated: 4/21/2010
*****************************************************/
#include "bgc_hmm.h"
#include "list_of_lists.h"
#include <rph_util.h>
#include <Rdefines.h>
#include <R_ext/Random.h>

SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);

//many more options to be added
//TODO: should either copy the gff or document clearly that it will be
// altered/ruined by this function
SEXP rph_bgc_hmm(SEXP msaP, 
		 SEXP modP, 
		 SEXP foregroundP,
		 SEXP doBgcP,
		 SEXP bgcP,
		 SEXP estimateBgcP,
		 SEXP bgcExpectedLengthP,
		 SEXP estimateBgcExpectedLengthP,
		 SEXP bgcTargetCoverageP,
		 SEXP estimateBgcTargetCoverageP,
		 SEXP selP,
		 SEXP consExpectedLengthP,
		 SEXP consTargetCoverageP,
		 SEXP estimateScaleP,
		 SEXP postProbsP) {
  struct bgchmm_struct *p = bgchmm_struct_new(1);
  
  GetRNGstate(); //seed R's random number generator
  p->msa = (MSA*)EXTPTR_PTR(msaP);
  msa_register_protect(p->msa);
  p->mod = (TreeModel*)EXTPTR_PTR(modP);
  tm_register_protect(p->mod);
  if (foregroundP != R_NilValue)
    p->foregd_branch = copy_charstr(CHARACTER_VALUE(foregroundP));
  p->do_bgc = LOGICAL_VALUE(doBgcP);
  p->bgc = NUMERIC_VALUE(bgcP);
  p->estimate_bgc = LOGICAL_VALUE(estimateBgcP);
  p->sel = NUMERIC_VALUE(selP);
  p->bgc_expected_length = NUMERIC_VALUE(bgcExpectedLengthP);
  p->estimate_bgc_expected_length = LOGICAL_VALUE(estimateBgcExpectedLengthP);
  p->bgc_target_coverage = NUMERIC_VALUE(bgcTargetCoverageP);
  p->estimate_bgc_target_coverage = LOGICAL_VALUE(estimateBgcTargetCoverageP);
  p->rho = 1.0;
  p->sel = NUMERIC_VALUE(selP);
  p->cons_expected_length = NUMERIC_VALUE(consExpectedLengthP);
  p->cons_target_coverage = NUMERIC_VALUE(consTargetCoverageP);
  p->estimate_scale = LOGICAL_VALUE(estimateScaleP);
  //no "WIG" option in rphast for now; simple enough to get this from posteriors anyway
  if (LOGICAL_VALUE(postProbsP)) {
    p->post_probs = FULL;
  } else p->post_probs = NONE;

  bgcHmm(p);

  PutRNGstate();
  return rph_listOfLists_to_SEXP(p->results);
}

SEXP rph_bgc_hmm_get_informative(SEXP msaP, 
				 SEXP modP, 
				 SEXP foregroundP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  TreeModel *mod = (TreeModel*)EXTPTR_PTR(modP);
  char *foregd = copy_charstr(CHARACTER_VALUE(foregroundP));
  ListOfLists *results = lol_new(1);
  int *bgc_informative;
  List *names = lst_new_ptr(mod->tree->nnodes/2);
  msa_register_protect(msa);
  tm_register_protect(mod);
  tm_prune(mod, msa, names);
  bgc_informative = bgchmm_get_informative(msa, foregd, mod->tree);
  bgchmm_print_informative(msa, bgc_informative, results, NULL, 0);
  return rph_listOfLists_to_SEXP(results);
}
