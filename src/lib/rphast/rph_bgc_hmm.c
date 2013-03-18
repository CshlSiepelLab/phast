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
		 SEXP gffP,
		 SEXP doBgcP,
		 SEXP nullModelP,
		 SEXP estimateBgcP,
		 SEXP estimateScaleP,
		 SEXP estimateRhoP,
		 SEXP estimateBgcTransP,
		 SEXP bgcExpectedLengthP,
		 SEXP bgcTargetCoverageP,
		 SEXP initBgcP,
		 SEXP initSelPosP,
		 SEXP initSelNegP,
		 SEXP initWeightsP,
		 SEXP initBgcInP,
		 SEXP initBgcOutP,
		 SEXP initScaleP,
		 SEXP initRhoP,
		 SEXP postProbsP,
		 SEXP randomPathP,
		 SEXP getLikelihoodsP,
		 SEXP informativeOnlyP) {
  struct bgchmm_struct *p = bgchmm_struct_new(1);
  double *doublep;
  
  GetRNGstate(); //seed R's random number generator
  p->msa = (MSA*)EXTPTR_PTR(msaP);
  msa_register_protect(p->msa);
  if (modP != R_NilValue) {
    p->mod = (TreeModel*)EXTPTR_PTR(modP);
    tm_register_protect(p->mod);
  } else {
    if (LOGICAL_VALUE(informativeOnlyP) != TRUE)
      die("bgchmm requires tree model\n");
  }
  if (foregroundP != R_NilValue)
    p->foregd_branch = copy_charstr(CHARACTER_VALUE(foregroundP));
  if (gffP != R_NilValue) {
    p->gff = (GFF_Set*)EXTPTR_PTR(gffP);
    gff_register_protect(p->gff);
  }
  p->do_bgc = LOGICAL_VALUE(doBgcP);
  p->null_model = LOGICAL_VALUE(nullModelP);
  p->estimate_bgc = LOGICAL_VALUE(estimateBgcP);
  p->estimate_scale = LOGICAL_VALUE(estimateScaleP);
  p->estimate_rho = LOGICAL_VALUE(estimateRhoP);
  p->estimate_transitions = LOGICAL_VALUE(estimateBgcTransP);
  p->bgc = NUMERIC_VALUE(initBgcP);
  p->selPos = NUMERIC_VALUE(initSelPosP);
  p->selNeg = NUMERIC_VALUE(initSelNegP);
  p->rho = NUMERIC_VALUE(initRhoP);
  doublep = NUMERIC_POINTER(initWeightsP);
  p->weights[0] = doublep[0];
  p->weights[1] = doublep[1];
  p->bgcInRate = NUMERIC_VALUE(initBgcInP);
  p->bgcOutRate = NUMERIC_VALUE(initBgcOutP);
  if (bgcTargetCoverageP != R_NilValue)
    p->bgcTargetCoverage = NUMERIC_VALUE(bgcTargetCoverageP);
  if (bgcExpectedLengthP != R_NilValue)
    p->bgcExpectedLength = NUMERIC_VALUE(bgcExpectedLengthP);
  p->scale = NUMERIC_VALUE(initScaleP);
  if (postProbsP != R_NilValue)
    p->post_probs = LOGICAL_VALUE(postProbsP);
  p->random_path = INTEGER_VALUE(randomPathP);
  p->get_likelihoods = LOGICAL_VALUE(getLikelihoodsP);
  p->informative_only = LOGICAL_VALUE(informativeOnlyP);
  bgcHmm(p);
  PutRNGstate();
  return rph_listOfLists_to_SEXP(p->results);
}
