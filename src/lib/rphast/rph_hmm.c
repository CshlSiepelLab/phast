/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_hmm.c
The RPHAST handles to functions dealing with multiple
sequence alignment functions from the phast package.

Melissa Hubisz
Last updated: 12/14/08
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <Rdefines.h>
#undef Matrix
#undef nrows
#include <hmm.h>
#include <math.h>
#include <rph_util.h>
#include <misc.h>
#include <markov_matrix.h>
#include <vector.h>
#include <list_of_lists.h>


void rph_hmm_protect(HMM *hmm) {
  int i;
  if (hmm == NULL) return;
  rph_mem_protect(hmm);
  rph_mm_protect(hmm->transition_matrix);
  rph_mat_protect(hmm->transition_score_matrix);
  rph_vec_protect(hmm->begin_transitions);
  rph_vec_protect(hmm->end_transitions);
  rph_vec_protect(hmm->eq_freqs);
  rph_vec_protect(hmm->begin_transition_scores);
  rph_vec_protect(hmm->end_transition_scores);
  for (i=0; i < hmm->nstates; i++) {
    rph_lst_protect(hmm->predecessors[i]);
    rph_lst_protect(hmm->successors[i]);
  }
  rph_mem_protect(hmm->predecessors);
  rph_mem_protect(hmm->successors);
  rph_lst_protect(hmm->begin_successors);
  rph_lst_protect(hmm->end_predecessors);
}


void rph_hmm_free(SEXP hmmP) {
  HMM *hmm = (HMM*)EXTPTR_PTR(hmmP);
  rph_unregister_protected(hmm);
  hmm_free(hmm);
}

void rph_hmm_register_protect(HMM *hmm) {
  rph_register_protected_object(hmm, (void (*)(void *))rph_hmm_protect);
}


SEXP rph_hmm_new_extptr(HMM *hmm) {
  SEXP result;
  rph_hmm_register_protect(hmm);
  PROTECT(result=R_MakeExternalPtr((void*)hmm, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, rph_hmm_free, 1);
  UNPROTECT(1);
  return result;
}

SEXP rph_hmm_new(SEXP matrixP, SEXP eqFreqP, SEXP beginFreqP,
		 SEXP endFreqP) {
  MarkovMatrix *mm;
  Matrix *m;
  int dim, i, j;
  double sum;
  Vector *eqFreq, *beginFreq, *endFreq=NULL;
  HMM *hmm;
  
  //matrix
  m = rph_get_matrix(matrixP);
  dim = m->nrows;
  for (i=0; i<dim; i++) {
    sum = 0.0;
    for (j=0; j<dim; j++) 
      sum += mat_get(m, i, j);
    for (j=0; j<dim; j++)
      mat_set(m, i, j, mat_get(m,i,j)/sum);
  }
  mm = mm_new_from_matrix(m, NULL, DISCRETE);

  eqFreq = rph_get_vector(eqFreqP);
  if (eqFreq->size != dim)
    die("bad dimension of eqFreqP in rph_hmm_new");
  vec_normalize(eqFreq);

  beginFreq = rph_get_vector(beginFreqP);
  if (beginFreq->size != dim)
    die("bad dimension of beginFreqP in rph_hmm_new");
  vec_normalize(beginFreq);

  if (endFreqP != R_NilValue) {
    endFreq = rph_get_vector(endFreqP);
    if (endFreq->size != dim)
      die("bad dimension of endFreqP in rph_hmm_new");
  }
  hmm = hmm_new(mm, eqFreq, beginFreq, endFreq);
  return rph_hmm_new_extptr(hmm);
}


SEXP rph_hmm_new_from_file(SEXP filenameP) {
  HMM *hmm;
  FILE *f = fopen_fname(CHARACTER_VALUE(filenameP), "r");
  hmm = hmm_new_from_file(f);
  fclose(f);
  return rph_hmm_new_extptr(hmm);
}


SEXP rph_hmm_print(SEXP hmmP, SEXP filenameP, SEXP appendP) {
  FILE *f=stdout;
  HMM *hmm;
  char *mode;
  if (filenameP != R_NilValue) {
    if (LOGICAL_VALUE(appendP))
      mode="a";
    else mode="w";
    f = fopen_fname(CHARACTER_VALUE(filenameP), mode);
  }
  hmm = (HMM*)EXTPTR_PTR(hmmP);
  hmm_print(f, hmm);
  if (f != stdout) fclose(f);
  return R_NilValue;
}
  
  


/* return the matrix */
SEXP rph_hmm_transMat(SEXP hmmP) {
  HMM *hmm = (HMM*)EXTPTR_PTR(hmmP);
  ListOfLists *lol;
  SEXP result;
  lol = lol_new(1);
  lol_push_matrix(lol, hmm->transition_matrix->matrix, "trans.mat");
  PROTECT(result = rph_listOfLists_to_SEXP(lol));
  UNPROTECT(1);
  return result;
}


SEXP rph_hmm_eqFreq(SEXP hmmP) {
  HMM *hmm = (HMM*)EXTPTR_PTR(hmmP);
  ListOfLists *lol;
  SEXP result;
  lol = lol_new(1);
  lol_push_dbl(lol, hmm->eq_freqs->data, hmm->eq_freqs->size, "eqFreq");
  PROTECT(result = rph_listOfLists_to_SEXP(lol));
  UNPROTECT(1);
  return result;
}


SEXP rph_hmm_beginFreq(SEXP hmmP) {
  HMM *hmm = (HMM*)EXTPTR_PTR(hmmP);
  ListOfLists *lol;
  SEXP result;
  lol = lol_new(1);
  lol_push_dbl(lol, hmm->begin_transitions->data, hmm->begin_transitions->size, 
	       "beginFreq");
  PROTECT(result = rph_listOfLists_to_SEXP(lol));
  UNPROTECT(1);
  return result;
}

SEXP rph_hmm_endFreq(SEXP hmmP) {
  HMM *hmm = (HMM*)EXTPTR_PTR(hmmP);
  ListOfLists *lol;
  SEXP result;
  if (hmm->end_transitions == NULL)
    return R_NilValue;
  lol = lol_new(1);
  lol_push_dbl(lol, hmm->end_transitions->data, hmm->end_transitions->size, 
	       "endFreq");
  PROTECT(result = rph_listOfLists_to_SEXP(lol));
  UNPROTECT(1);
  return result;
}


