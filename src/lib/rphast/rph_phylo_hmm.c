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

void rph_phmm_protect(PhyloHmm *p) {
  int i;
  if (p == NULL) return;
  rph_mem_protect(p);
  rph_cm_protect(p->cm);
  rph_hmm_protect(p->hmm);
  rph_hmm_protect(p->functional_hmm);
  rph_hmm_protect(p->autocorr_hmm);
  if (p->mods != NULL) {
    rph_mem_protect(p->mods);
    for (i=0; i < p->nmods; i++)
      rph_tm_protect(p->mods[i]);
  }
  rph_gp_protect(p->gpm);
  if (p->state_to_mod != NULL)
    rph_mem_protect(p->state_to_mod);
  if (p->state_to_cat != NULL)
    rph_mem_protect(p->state_to_cat);
  if (p->reverse_compl != NULL)
    rph_mem_protect(p->reverse_compl);
  if (p->cat_to_states != NULL) {
    rph_mem_protect(p->cat_to_states);
    for (i=0; i <= p->cm->ncats; i++) 
      rph_lst_protect(p->cat_to_states[i]);
  }
  if (p->state_to_pattern != NULL)
    rph_mem_protect(p->state_to_pattern);
  if (p->emissions != NULL) {
    for (i=0; i < p->hmm->nstates; i++)
      if (p->state_pos[p->state_to_mod[i]] == i ||
          p->state_neg[p->state_to_mod[i]] == i || 
          p->state_to_pattern[i] >= 0)
        rph_mem_protect(p->emissions[i]);
    rph_mem_protect(p->emissions);
    rph_mem_protect(p->state_pos);
    rph_mem_protect(p->state_neg);
  }
  if (p->forward != NULL) {
    for (i=0; i < p->hmm->nstates; i++) 
      rph_mem_protect(p->forward[i]);
    rph_mem_protect(p->forward);
  }
  //topology isn't freed by phmm_free so shouldn't be protected?
  //  if (p->topology != NULL)
  //    rph_tree_protect(p->topology);

  if (p->alpha != NULL) 
    rph_mem_protect(p->alpha);
  if (p->beta != NULL)
    rph_mem_protect(p->beta);
  if (p->tau != NULL)
    rph_mem_protect(p->tau);
  if (p->T != NULL) {
    rph_mem_protect(p->T);
    rph_mem_protect(p->t);
    for (i=0; i < p->functional_hmm->nstates; i++) {
      rph_mem_protect(p->T[i]);
      rph_mem_protect(p->t[i]);
    }
  }
  if (p->em_data != NULL) {
    rph_mem_protect(p->em_data);
    rph_mat_protect(p->em_data->H);
    //rph_msa_protect(p->em_data->msa);  //assume this is protected elsewhere/
  }

}


void rph_phmm_register_protect(PhyloHmm *phmm) {
  rph_register_protected_object(phmm, (void (*)(void *))rph_phmm_protect);
}


void rph_phmm_free(SEXP phmmP) {
  PhyloHmm *phmm = (PhyloHmm*)EXTPTR_PTR(phmmP);
  rph_unregister_protected(phmmP);
  phmm_free(phmm);
}


SEXP rph_phmm_new_extptr(PhyloHmm *phmm) {
  SEXP result;
  rph_phmm_register_protect(phmm);
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
    rph_tm_register_protect(mods[i]);
  }
  hmm = (HMM*)EXTPTR_PTR(hmmP);
  rph_hmm_register_protect(hmm);

  //sending pivot_states to phmm_new automatically reflects
  phmm = phmm_new(hmm, mods, NULL, pivot_states, MISSING_DATA);

  return rph_phmm_new_extptr(phmm);
}
