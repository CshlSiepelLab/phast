/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_tree_model.c
The RPHAST handles to functions dealing with multiple
sequence alignment functions from the phast package.

Melissa Hubisz
Last updated: 1/13/10
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <tree_model.h>
#include <matrix.h>
#include <tree_likelihoods.h>
#include <misc.h>
#include <Rdefines.h>

//these are defined as macros in R and we don't want them overriding
// phast's matrix->nrows and matrix->ncols
#undef nrows
#undef ncols

void rph_tm_free(SEXP tmP) {
  tm_free((TreeModel*)EXTPTR_PTR(tmP));
}


SEXP rph_tm_new_extptr(TreeModel *tm) {
  SEXP result;
  PROTECT(result=R_MakeExternalPtr((void*)tm, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, rph_tm_free, 1);
  UNPROTECT(1);
  return result;
}



SEXP rph_tm_tree(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  char *treeStr;
  SEXP result;
  
  treeStr = tr_to_string(tm->tree, 1);
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(treeStr));
  free(treeStr);
  UNPROTECT(1);
  return result;
}



SEXP rph_tm_alphabet(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  if (tm->rate_matrix == NULL ||
      tm->rate_matrix->states == NULL) return R_NilValue;
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(tm->rate_matrix->states));
  UNPROTECT(1);
  return result;
}


SEXP rph_tm_backgd(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  double *resultP;
  int i;

  if (tm->backgd_freqs == NULL)
    return R_NilValue;

  PROTECT(result = NEW_NUMERIC(tm->backgd_freqs->size));
  resultP = NUMERIC_POINTER(result);
  for (i=0; i<tm->backgd_freqs->size; i++) 
    resultP[i] = vec_get(tm->backgd_freqs, i);
  UNPROTECT(1);
  return result;
}


SEXP rph_tm_rateMatrix(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  double *resultP;
  int i, j, pos=0;
  Matrix *m;

  if (tm->rate_matrix == NULL || tm->rate_matrix->matrix == NULL)
    return R_NilValue;
  m = tm->rate_matrix->matrix;
  if (m->nrows != tm->rate_matrix->size ||
      m->ncols != tm->rate_matrix->size)
    die("invalid rate matrix dimensions");
  PROTECT(result = NEW_NUMERIC(m->nrows * m->ncols));
  resultP = NUMERIC_POINTER(result);
  pos=0;
  for (i=0; i<m->nrows; i++)
    for (j=0;  j<m->ncols; j++)
      resultP[pos++] = mat_get(m, i, j);
  UNPROTECT(1);
  return result;
}


SEXP rph_tm_substMod(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(tm_get_subst_mod_string(tm->subst_mod)));
  UNPROTECT(1);
  return result;
}


SEXP rph_tm_order(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  int *resultP;
  
  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  resultP[0] = tm->order;
  UNPROTECT(1);
  return result;
}

SEXP rph_tm_likelihood(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  double *resultP;
  if (tm->lnL == NULL_LOG_LIKELIHOOD)
    return R_NilValue;
  PROTECT(result = NEW_NUMERIC(1));
  resultP = NUMERIC_POINTER(result);
  resultP[0] = tm->lnL;
  UNPROTECT(1);
  return result;
}

SEXP rph_tm_empirical_rates(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  int *resultP;
  
  PROTECT(result = NEW_LOGICAL(1));
  resultP = LOGICAL_POINTER(result);
  resultP[0] = tm->empirical_rates;
  UNPROTECT(1);
  return result;
}



SEXP rph_tm_alpha(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  double *resultP;
  
  PROTECT(result = NEW_NUMERIC(1));
  resultP = NUMERIC_POINTER(result);
  resultP[0] = tm->alpha;
  UNPROTECT(1);
  return result;
}

SEXP rph_tm_nratecats(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  int *resultP;

  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  resultP[0] = tm->nratecats;
  UNPROTECT(1);
  return result;
}

SEXP rph_tm_rK(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  double *resultP;
  int i;

  if (tm->rK == NULL || tm->empirical_rates==0)
    return R_NilValue;
  PROTECT(result = NEW_NUMERIC(tm->nratecats));
  resultP = NUMERIC_POINTER(result);
  for (i=0; i<tm->nratecats; i++)
    resultP[i] = tm->rK[i];
  UNPROTECT(1);
  return result;
}

SEXP rph_tm_freqK(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  double *resultP;
  int i;
  
  if (tm->freqK == NULL || tm->empirical_rates==0)
    return R_NilValue;
  PROTECT(result = NEW_NUMERIC(tm->nratecats));
  resultP = NUMERIC_POINTER(result);
  for (i=0; i<tm->nratecats; i++)
    resultP[i] = tm->freqK[i];
  UNPROTECT(1);
  return result;
}
  

SEXP rph_tm_rootLeaf(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  int *resultP;

  if (tm->root_leaf_id == -1) 
    return R_NilValue;
  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  resultP[0] = tm->root_leaf_id;
  UNPROTECT(1);
  return result;
}



SEXP rph_tm_new(SEXP treeP, SEXP alphabetP, SEXP backgdP, SEXP matrixP, 
		SEXP substModP, SEXP lnlP, SEXP alphaP, SEXP nratecatsP,
		SEXP rKP, SEXP freqKP, SEXP rootLeafP) {
  TreeModel *tm;
  TreeNode *tree;
  MarkovMatrix *rateMatrix;
  Matrix *m;
  int dim=-1, i, j, pos, numProtect=0, nratecats=1, rootLeaf;
  double *doubleP, alpha=0.0;
  Vector *backgd;
  char *alphabet;
  subst_mod_type subst_mod;
  List *rate_consts=NULL;

  //tree
  tree = tr_new_from_string(CHARACTER_VALUE(treeP));

  //alphabet
  if (alphabetP == R_NilValue)
    die("alphabet cannot be NULL");
  alphabet = smalloc((strlen(CHARACTER_VALUE(alphabetP))+1)*sizeof(char));
  strcpy(alphabet, CHARACTER_VALUE(alphabetP));
  dim = strlen(alphabet);

  //backgd
  if (backgdP == R_NilValue)
    backgd = NULL;
  else {
    PROTECT(backgdP = AS_NUMERIC(backgdP));
    numProtect++;
    doubleP = NUMERIC_POINTER(backgdP);
    if (LENGTH(backgdP) != dim)
      die("backgd frequency size does not match alphabet");
    backgd = vec_new(dim);
    for (i=0; i<dim; i++)
      vec_set(backgd, i, doubleP[i]);
  }

  //matrix
  if (matrixP == R_NilValue)
    rateMatrix = NULL;
  else {
    PROTECT(matrixP = AS_NUMERIC(matrixP));
    numProtect++;
    doubleP = NUMERIC_POINTER(matrixP);
    if (dim*dim != LENGTH(matrixP))
      die("size of rate matrix does not match alphabet size");
    m = mat_new(dim, dim);
    pos = 0;
    for (i=0; i<dim; i++)
      for (j=0; j<dim; j++)
	mat_set(m, j, i, doubleP[pos++]);
    rateMatrix = mm_new_from_matrix(m, alphabet, CONTINUOUS);
  }

  //subst mod
  subst_mod = tm_get_subst_mod_type(CHARACTER_VALUE(substModP));
  if (subst_mod == UNDEF_MOD) 
    die("invalid subst mod %s", CHARACTER_VALUE(substModP));

  //alpha
  if (alphaP != R_NilValue) 
    alpha = NUMERIC_VALUE(alphaP);

  //nratecats
  if (nratecatsP != R_NilValue)
    nratecats = INTEGER_VALUE(nratecatsP);

  //rK (stored in rate_consts)
  if (rKP != R_NilValue) {
    PROTECT(rKP = AS_NUMERIC(rKP));
    numProtect++;
    doubleP = NUMERIC_POINTER(rKP);
    if (LENGTH(rKP) != nratecats)
      die("rK should be NULL or length equal to nratecats");
    rate_consts = lst_new_dbl(nratecats);
    for (i=0; i<nratecats; i++)
      lst_push_dbl(rate_consts, doubleP[i]);
  }

  //rootLeafP
  if (rootLeafP == R_NilValue) rootLeaf = -1;
  else {
    TreeNode *n = tr_get_node(tree, CHARACTER_VALUE(rootLeafP));
    if (n == NULL) 
      die("no node named %s", CHARACTER_VALUE(rootLeafP));
    rootLeaf = n->id;
  }

  tm = tm_new(tree, rateMatrix, backgd, subst_mod,
	      alphabet, nratecats, alpha, rate_consts, 
	      rootLeaf);


  if (freqKP != R_NilValue) {
    PROTECT(freqKP = AS_NUMERIC(freqKP));
    numProtect++;
    doubleP = NUMERIC_POINTER(freqKP);
    if (LENGTH(freqKP) != tm->nratecats)
      die("length of rate.weights should equal nratecats");
    if (tm->freqK == NULL)
      tm->freqK = smalloc(tm->nratecats * sizeof(double));
    for (i=0; i<tm->nratecats; i++)
      tm->freqK[i] = doubleP[i];
    normalize_probs(tm->freqK, tm->nratecats);
  }

  if (lnlP != R_NilValue)
    tm->lnL = NUMERIC_VALUE(lnlP);

  //these get copied by tm_new
  free(alphabet);
  if (rate_consts != NULL) lst_free(rate_consts);

  if (numProtect > 0)
    UNPROTECT(numProtect);
  return rph_tm_new_extptr(tm);
}
 

SEXP rph_tm_print(SEXP tmP, SEXP filenameP, SEXP appendP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  FILE *outfile;
  char *mode = "w";
  if (filenameP == R_NilValue)
    outfile = stdout;
  else {
    if (LOGICAL_VALUE(appendP)) 
      mode = "a";
    outfile = fopen_fname(CHARACTER_VALUE(filenameP), mode);
  }
  tm_print(outfile, tm);
  if (outfile != stdout) fclose(outfile);
  return R_NilValue;
}


SEXP rph_tm_read(SEXP filenameP) {
  FILE *infile;
  TreeModel *tm;
  if (filenameP == R_NilValue)
    die("filename cannot be NULL");
  infile = fopen_fname(CHARACTER_VALUE(filenameP), "r");
  tm = tm_new_from_file(infile);
  fclose(infile);
  return rph_tm_new_extptr(tm);
}


