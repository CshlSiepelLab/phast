/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
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
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <trees.h>
#include <phast_cons.h>
#include <Rdefines.h>

SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);

SEXP rph_phastCons(SEXP msaP, SEXP modP, SEXP rhoP, SEXP estimateTreesP,
		   SEXP estimateRhoP, SEXP gcP, SEXP nratesP, 
		   SEXP transitionsP, SEXP initTransitionsP,
		   SEXP targetCoverageP, SEXP expectedLengthP,
		   SEXP initExpectedLengthP, SEXP viterbiP,
		   SEXP scoreViterbiP, SEXP computeLnlP, 
		   SEXP suppressProbsP, SEXP refIdxP) {
  struct phastCons_struct *p = phastCons_struct_new(1);
  int i, *intp, numprotect=0;
  double *doublep;
  SEXP rv;

  p->msa = (MSA*)EXTPTR_PTR(msaP);
  p->mod = (TreeModel**)smalloc(LENGTH(modP)*sizeof(TreeModel*));
  for (i=0; i<LENGTH(modP); i++) {
    p->mod[i]=(TreeModel*)EXTPTR_PTR(VECTOR_ELT(modP, i));
  }
  p->nummod = LENGTH(modP);
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
    p->estim_transitions = FALSE;
    p->mu = doublep[0];
    p->nu = doublep[1];
  }
  if (initTransitionsP != R_NilValue) {
    doublep = NUMERIC_POINTER(initTransitionsP);
    p->set_transitions = TRUE;
    p->estim_transitions=TRUE;  //this is the default anyway
    p->mu = doublep[0];
    p->nu = doublep[1];
  }
  if (targetCoverageP != R_NilValue) 
    p->gamma = NUMERIC_VALUE(targetCoverageP);
  if (expectedLengthP != R_NilValue) {
    p->estim_transitions = FALSE;
    p->omega = NUMERIC_VALUE(expectedLengthP);
    p->mu = 1.0/p->omega;
  }
  if (initExpectedLengthP != R_NilValue) {
    p->estim_transitions = TRUE;  //default anyway
    p->omega = NUMERIC_VALUE(initExpectedLengthP);
    p->mu = 1.0/p->omega;
  }
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

  phastCons(p);

  //free anything allocated?
  //close outfiles
  if (p->results != NULL) {
    PROTECT(rv = rph_listOfLists_to_SEXP(p->results));
    numprotect++;
    lol_free(p->results);
  } else rv=R_NilValue;
  
  free(p);
  if (numprotect > 0) UNPROTECT(numprotect);
  return rv;
}
