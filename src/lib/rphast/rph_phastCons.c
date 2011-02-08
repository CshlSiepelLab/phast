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
#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <getopt.h>
#include <ctype.h>
#include <misc.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <trees.h>
#include <phast_cons.h>
#include <hmm.h>
#include <rph_util.h>
#include <Rdefines.h>
#include <R_ext/Random.h>

SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);

SEXP rph_phastCons(SEXP msaP, 
		   SEXP modP, 
		   SEXP rhoP, 
		   SEXP targetCoverageP,
		   SEXP expectedLengthP, 
		   SEXP transitionsP,
		   SEXP estimateRhoP,
		   SEXP estimateExpectedLengthP,
		   SEXP estimateTransitionsP,
		   SEXP estimateTreesP,
		   SEXP viterbiP, 
		   SEXP scoreViterbiP,
		   SEXP gcP, 
		   SEXP nratesP, 
		   SEXP computeLnlP, 
		   SEXP suppressProbsP,
		   SEXP refIdxP,
		   SEXP hmmP, 
		   SEXP statesP, 
		   SEXP reflectStrandP,
		   SEXP quietP,
		   SEXP categoryMapP) {
  struct phastCons_struct *p = phastCons_struct_new(1);
  int i, *intp, numprotect=0;
  double *doublep;
  SEXP rv;
  
  GetRNGstate(); //seed R's random number generator
  p->msa = (MSA*)EXTPTR_PTR(msaP);
  p->nummod = LENGTH(modP);
  p->mod = (TreeModel**)smalloc(p->nummod*sizeof(TreeModel*));
  for (i=0; i<p->nummod; i++) {
    p->mod[i]=(TreeModel*)EXTPTR_PTR(VECTOR_ELT(modP, i));
    p->mod[i]->use_conditionals = 1;
  }
  if (rhoP != R_NilValue) 
    p->rho = NUMERIC_VALUE(rhoP);
  if (estimateTreesP != R_NilValue) {
    p->estim_trees = LOGICAL_VALUE(estimateTreesP);
  }
  if (estimateRhoP != R_NilValue) {
    p->estim_rho = LOGICAL_VALUE(estimateRhoP);
  }
  if (gcP != R_NilValue) {
    p->gc = NUMERIC_VALUE(gcP);
  }
  if (nratesP != R_NilValue) {
    intp = INTEGER_POINTER(nratesP);
    p->nrates = intp[0];
    if (LENGTH(nratesP)==2)
      p->nrates2 = intp[1];
  }
  if (transitionsP != R_NilValue) {
    doublep = NUMERIC_POINTER(transitionsP);
    p->set_transitions = TRUE;
    p->estim_transitions = LOGICAL_VALUE(estimateTransitionsP);
    p->mu = doublep[0];
    p->nu = doublep[1];
  }
  if (targetCoverageP != R_NilValue) 
    p->gamma = NUMERIC_VALUE(targetCoverageP);
  if (expectedLengthP != R_NilValue) {
    p->omega = NUMERIC_VALUE(expectedLengthP);
    p->mu = 1.0/p->omega;
  }
  p->estim_transitions = FALSE;
  if ((estimateExpectedLengthP != R_NilValue && LOGICAL_VALUE(estimateExpectedLengthP)) ||
      (estimateTransitionsP != R_NilValue && LOGICAL_VALUE(estimateTransitionsP)))
    p->estim_transitions = TRUE;

  if (viterbiP != R_NilValue) 
    p->viterbi = LOGICAL_VALUE(viterbiP);
  if (scoreViterbiP != R_NilValue) {
    p->score = LOGICAL_VALUE(scoreViterbiP);
  }
  if (computeLnlP != R_NilValue) {
    if (LOGICAL_VALUE(computeLnlP)) {
      p->compute_likelihood = TRUE;
    }
  }
  if (suppressProbsP != R_NilValue) 
    p->post_probs = !LOGICAL_VALUE(suppressProbsP);
  if (refIdxP != R_NilValue) 
    p->refidx = INTEGER_VALUE(refIdxP);

  p->seqname = p->msa->names[p->refidx-1];

  if (hmmP != R_NilValue) {
    p->hmm = (HMM*)EXTPTR_PTR(hmmP);
    p->two_state = FALSE;
    p->nummod = p->hmm->nstates;
    hmm_register_protect(p->hmm);
  }
  
  if (statesP != R_NilValue) {
    char tempstr[100];
    p->states = lst_new_ptr(LENGTH(statesP));
    intp = INTEGER_POINTER(statesP);
    for (i=0; i < LENGTH(statesP); i++) {
      sprintf(tempstr, "%i", intp[i]-1);
      lst_push_ptr(p->states, str_new_charstr(tempstr));
    }
  }

  if (reflectStrandP != R_NilValue) {
    char tempstr[100];
    p->pivot_states = lst_new_ptr(LENGTH(reflectStrandP));
    intp = INTEGER_POINTER(reflectStrandP);
    for (i=0; i < LENGTH(reflectStrandP); i++) {
      sprintf(tempstr, "%i", intp[i]-1);
      lst_push_ptr(p->pivot_states, str_new_charstr(tempstr));
    }
  }
  if (LOGICAL_VALUE(quietP))
    p->results_f = NULL;

  if (categoryMapP != R_NilValue)
    p->cm = cm_new_string_or_file(CHARACTER_VALUE(categoryMapP));

  msa_register_protect(p->msa);
  for (i=0; i < LENGTH(modP); i++)
    tm_register_protect((TreeModel*)EXTPTR_PTR(VECTOR_ELT(modP, i)));

  phastCons(p);

  if (p->results != NULL) {
    PROTECT(rv = rph_listOfLists_to_SEXP(p->results));
    numprotect++;
  } else rv=R_NilValue;
  
  PutRNGstate();
  if (numprotect > 0) UNPROTECT(numprotect);
  return rv;
}
