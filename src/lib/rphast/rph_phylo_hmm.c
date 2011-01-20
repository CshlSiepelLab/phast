/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_phylo_hmm.c
The RPHAST handles to functions dealing with
phylo-hmm funtions from the phast package.

Melissa Hubisz
Last updated: 10/27/2010
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <Rdefines.h>
#undef Matrix
#include <hmm.h>
#include <math.h>
#include <rph_util.h>
#include <misc.h>
#include <markov_matrix.h>
#include <vector.h>
#include <list_of_lists.h>
#include <phylo_hmm.h>


void rph_phmm_free(SEXP phmmP) {
  PhyloHmm *phmm = (PhyloHmm*)EXTPTR_PTR(phmmP);
  phast_new_mem_handler();  //new memory handler needed because phmm_free invokes tm_free which invokes tr_free which allocates memory
  phast_unregister_protected(phmmP);
  phmm_free(phmm);
  phast_free_all();
}


SEXP rph_phmm_new_extptr(PhyloHmm *phmm) {
  SEXP result;
  phmm_register_protect(phmm);
  PROTECT(result=R_MakeExternalPtr((void*)phmm, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, rph_phmm_free, 1);
  UNPROTECT(1);
  return result;
}


SEXP rph_phmm_get_treeModel(SEXP phmmP, SEXP whichP) {
  PhyloHmm *phmm = (PhyloHmm*)EXTPTR_PTR(phmmP);
  int which = INTEGER_VALUE(whichP)-1;
  if (phmm==NULL) die("phmm is NULL");
  if (phmm->mods == NULL) die("phmm->mods is NULL");
  if (which < 0 || which >= phmm->hmm->nstates)
    die("which is out of bounds (%i, nstates=%i)\n", which, phmm->hmm->nstates);
  return rph_tm_new_extptr(tm_create_copy(phmm->mods[phmm->state_to_mod[which]]));
}


SEXP rph_phmm_get_hmm(SEXP phmmP) {
  PhyloHmm* phmm = (PhyloHmm*)EXTPTR_PTR(phmmP);
  return rph_hmm_new_extptr(phmm->hmm);
}


SEXP rph_phmm_get_state_to_mod(SEXP phmmP) {
  PhyloHmm *phmm = (PhyloHmm*)EXTPTR_PTR(phmmP);
  SEXP rv;
  int *intp, i;
  PROTECT(rv = NEW_INTEGER(phmm->hmm->nstates));
  intp = INTEGER_POINTER(rv);
  for (i=0; i < phmm->hmm->nstates; i++) 
    intp[i] = phmm->state_to_mod[i];
  UNPROTECT(1);
  return rv;
}


SEXP rph_phmm_reflect_strand(SEXP hmmP, SEXP pivotStatesP, SEXP modsP) {
  TreeModel **mods;
  List *pivot_states;
  int *intp, i;
  HMM *hmm;
  PhyloHmm *phmm;
  char tempstr[100];

  pivot_states = lst_new_ptr(LENGTH(pivotStatesP));
  intp = INTEGER_POINTER(pivotStatesP);
  for (i=0; i < LENGTH(pivotStatesP); i++) {
    sprintf(tempstr, "%i", intp[i]-1);
    lst_push_ptr(pivot_states, str_new_charstr(tempstr));
  }
  mods = smalloc(LENGTH(modsP)*sizeof(TreeModel*));
  for (i=0; i < LENGTH(modsP); i++) {
    mods[i] = (TreeModel*)EXTPTR_PTR(VECTOR_ELT(modsP, i));
    tm_register_protect(mods[i]);
  }
  hmm = (HMM*)EXTPTR_PTR(hmmP);
  hmm_register_protect(hmm);

  //sending pivot_states to phmm_new automatically reflects
  phmm = phmm_new(hmm, mods, NULL, pivot_states, MISSING_DATA);

  return rph_phmm_new_extptr(phmm);
}
