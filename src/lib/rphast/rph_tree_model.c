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
#include <getopt.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <tree_model.h>
#include <matrix.h>
#include <tree_likelihoods.h>
#include <list_of_lists.h>
#include <misc.h>
#include <Rdefines.h>
#include <rph_util.h>

//these are defined as macros in R and we don't want them overriding
// phast's matrix->nrows and matrix->ncols
#undef nrows
#undef ncols

SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);

void rph_tm_free(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  rph_unregister_protected(tm);
  tm_free(tm);
}


void rph_tm_rmp_protect(TreeModel *tm) {
  int nparams = tm_get_nparams(tm);
  int i;
  for (i=0; i < nparams; i++) {
    if (tm->rate_matrix_param_row[i] != NULL) {
      rph_lst_protect(tm->rate_matrix_param_row[i]);
      rph_lst_protect(tm->rate_matrix_param_col[i]);
    }
  }
  rph_mem_protect(tm->rate_matrix_param_row);
  rph_mem_protect(tm->rate_matrix_param_col);
}

void rph_tm_altmod_protect(AltSubstMod *am) {
  int i;
  rph_mem_protect(am);
  if (am->backgd_freqs != NULL)
    rph_vec_protect(am->backgd_freqs);
  if (am->rate_matrix != NULL)
    rph_mm_protect(am->rate_matrix);
  if (am->param_list != NULL) {
    rph_lst_protect(am->param_list);
    for (i=0; i < lst_size(am->param_list); i++)
      rph_str_protect(lst_get_ptr(am->param_list, i));
  }
  rph_str_protect(am->defString);
  if (am->noopt_arg != NULL)
    rph_str_protect(am->noopt_arg);
}


void rph_tm_protect(TreeModel *tm) {
  int i, j;
  rph_mem_protect(tm);
  rph_tree_protect(tm->tree);
  rph_vec_protect(tm->backgd_freqs);
  rph_mm_protect(tm->rate_matrix);
  if (tm->msa_seq_idx != NULL) rph_mem_protect(tm->msa_seq_idx);
  if (tm->P != NULL) {
    for (i=0; i < tm->tree->nnodes; i++) {
      for (j=0; j < tm->nratecats; j++)
	rph_mm_protect(tm->P[i][j]);
      rph_mem_protect(tm->P[i]);
    }
    rph_mem_protect(tm->P);
  }
  if (tm->rK != NULL) rph_mem_protect(tm->rK);
  if (tm->freqK != NULL) rph_mem_protect(tm->freqK);
  if (tm->rate_matrix_param_row != NULL) 
    rph_tm_rmp_protect(tm);
  if (tm->in_subtree != NULL) rph_mem_protect(tm->in_subtree);
  if (tm->ignore_branch != NULL) rph_mem_protect(tm->ignore_branch);
  if (tm->alt_subst_mods != NULL) {
    rph_lst_protect(tm->alt_subst_mods);
    for (i=0; i < lst_size(tm->alt_subst_mods); i++)
      rph_tm_altmod_protect(lst_get_ptr(tm->alt_subst_mods, i));
  }
  if (tm->alt_subst_mods_node != NULL)
    rph_mem_protect(tm->alt_subst_mods_node);
  if (tm->all_params != NULL) rph_vec_protect(tm->all_params);
  if (tm->param_map != NULL) rph_mem_protect(tm->param_map);
  if (tm->bound_arg != NULL) {
    rph_lst_protect(tm->bound_arg);
    for (i=0; i < lst_size(tm->bound_arg); i++)
      rph_str_protect(lst_get_ptr(tm->bound_arg, i));
  }
  if (tm->noopt_arg != NULL)
    rph_str_protect(tm->noopt_arg);
  if (tm->iupac_inv_map != NULL) {
    for (i=0; i < 256; i++)
      if (tm->iupac_inv_map[i] != NULL) rph_mem_protect(tm->iupac_inv_map[i]);
    rph_mem_protect(tm->iupac_inv_map);
  }
}

void rph_tm_register_protect(TreeModel *tm) {
  rph_register_protected_object(tm, (void (*)(void*))rph_tm_protect);
}


SEXP rph_tm_new_extptr(TreeModel *tm) {
  SEXP result;
  rph_tm_register_protect(tm);
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


SEXP rph_tm_altmodel_backgd(SEXP tmP, SEXP whichmodP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  AltSubstMod *altmod;
  SEXP result;
  double *resultP;
  int i, whichmod = INTEGER_VALUE(whichmodP);

  if (tm->alt_subst_mods == NULL) 
    die("No alt subst mods in this treeModel");
  if (lst_size(tm->alt_subst_mods) < whichmod)
    die("Not enough alt subst mods in this treeModel");
  altmod = lst_get_ptr(tm->alt_subst_mods, whichmod-1);

  if (altmod->backgd_freqs == NULL)
    return R_NilValue;
  
  PROTECT(result = NEW_NUMERIC(altmod->backgd_freqs->size));
  resultP = NUMERIC_POINTER(result);
  for (i=0; i<altmod->backgd_freqs->size; i++) 
    resultP[i] = vec_get(altmod->backgd_freqs, i);
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


SEXP rph_tm_num_altmodel(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  PROTECT(result = NEW_INTEGER(1));
  int *resultP = INTEGER_POINTER(result);
  if (tm->alt_subst_mods == NULL)
    resultP[0] = 0;
  else resultP[0] = lst_size(tm->alt_subst_mods);
  UNPROTECT(1);
  return result;
}


SEXP rph_tm_altmodel_sel(SEXP tmP, SEXP whichmodP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  AltSubstMod *altmod;
  double *resultP;
  SEXP result;
  int whichmod = INTEGER_VALUE(whichmodP);

  if (tm->alt_subst_mods == NULL) 
    die("No alt subst mods in this treeModel");
  if (lst_size(tm->alt_subst_mods) < whichmod)
    die("Not enough alt subst mods in this treeModel");
  altmod = lst_get_ptr(tm->alt_subst_mods, whichmod-1);
  if  (altmod->selection_idx < 0) return R_NilValue;
  PROTECT(result = NEW_NUMERIC(1));
  resultP = NUMERIC_POINTER(result);
  resultP[0] = altmod->selection;
  UNPROTECT(1);
  return result;
}


SEXP rph_tm_altmodel_bgc(SEXP tmP, SEXP whichmodP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  AltSubstMod *altmod;
  double *resultP;
  SEXP result;
  int whichmod = INTEGER_VALUE(whichmodP);

  if (tm->alt_subst_mods == NULL) 
    die("No alt subst mods in this treeModel");
  if (lst_size(tm->alt_subst_mods) < whichmod)
    die("Not enough alt subst mods in this treeModel");
  altmod = lst_get_ptr(tm->alt_subst_mods, whichmod-1);
  if (altmod->bgc_idx < 0) return R_NilValue;

  PROTECT(result = NEW_NUMERIC(1));
  resultP = NUMERIC_POINTER(result);
  resultP[0] = altmod->bgc;
  UNPROTECT(1);
  return result;
}


SEXP rph_tm_altmodel_def(SEXP tmP, SEXP whichmodP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  AltSubstMod *altmod;
  int whichmod = INTEGER_VALUE(whichmodP);
  SEXP result;
  
  if (tm->alt_subst_mods == NULL) 
    die("No alt subst mods in this treeModel");
  if (lst_size(tm->alt_subst_mods) < whichmod)
    die("Not enough alt subst mods in this treeModel");
  altmod = lst_get_ptr(tm->alt_subst_mods, whichmod-1);

  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(altmod->defString->chars));
  UNPROTECT(1);
  return result;
}


SEXP rph_tm_altmodel_rateMatrix(SEXP tmP, SEXP whichmodP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  AltSubstMod *altmod;
  ListOfLists *lol;
  int whichmod = INTEGER_VALUE(whichmodP);
  SEXP result;

  if (tm->alt_subst_mods == NULL) 
    die("No alt subst mods in this treeModel");
  if (lst_size(tm->alt_subst_mods) < whichmod)
    die("Not enough alt subst mods in this treeModel");
  altmod = lst_get_ptr(tm->alt_subst_mods, whichmod-1);

  if (altmod->rate_matrix == NULL || altmod->rate_matrix->matrix == NULL)
    return R_NilValue;
  lol = lol_new(1);
  lol_push_matrix(lol, altmod->rate_matrix->matrix, "rate.matrix");
  PROTECT(result = rph_listOfLists_to_SEXP(lol));
  UNPROTECT(1);
  return result;
}

SEXP rph_tm_rateMatrix(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  ListOfLists *lol;
  SEXP result;

  if (tm->rate_matrix == NULL || tm->rate_matrix->matrix == NULL)
    return R_NilValue;
  lol = lol_new(1);
  lol_push_matrix(lol, tm->rate_matrix->matrix, "rate.matrix");
  PROTECT(result = rph_listOfLists_to_SEXP(lol));
  UNPROTECT(1);
  return result;
}



SEXP rph_tm_altmodel_substMod(SEXP tmP, SEXP whichmodP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  AltSubstMod *altmod;
  int whichmod = INTEGER_VALUE(whichmodP);
  SEXP result;
  if (tm->alt_subst_mods == NULL) 
    die("No alt subst mods in this treeModel");
  if (lst_size(tm->alt_subst_mods) < whichmod)
    die("Not enough alt subst mods in this treeModel");
  altmod = lst_get_ptr(tm->alt_subst_mods, whichmod-1);
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(tm_get_subst_mod_string(altmod->subst_mod)));
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


//returns a numeric vector of length 2.  The first element will
// be zero if there is no selection paramer.  If there is selection,
// the first element will be one and the second will be the parameter.
SEXP rph_tm_selection(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  SEXP result;
  double *resultP;
  PROTECT(result = NEW_NUMERIC(2));
  resultP = NUMERIC_POINTER(result);
  if  (tm->selection_idx >= 0) {
    resultP[0] = 1.0;
    resultP[1] = tm->selection;
  } else resultP[0] = 0.0;
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
		SEXP rKP, SEXP freqKP, SEXP rootLeafP,
		SEXP selectionP) {
  TreeModel *tm;
  TreeNode *tree;
  MarkovMatrix *rateMatrix;
  Matrix *m;
  int dim=-1, i, numProtect=0, nratecats=1, rootLeaf;
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

  //subst mod
  subst_mod = tm_get_subst_mod_type(CHARACTER_VALUE(substModP));
  if (subst_mod == UNDEF_MOD) 
    die("invalid subst mod %s", CHARACTER_VALUE(substModP));

  dim = int_pow(strlen(alphabet), tm_order(subst_mod)+1);

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
    m = rph_get_matrix(matrixP);
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

  if (selectionP != R_NilValue) {
    tm->selection = NUMERIC_VALUE(selectionP);
    tm->selection_idx = 0;
  }


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
  tm = tm_new_from_file(infile, 0);
  fclose(infile);
  return rph_tm_new_extptr(tm);
}



SEXP rph_tm_add_alt_mod(SEXP tmP, SEXP defStrP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  String *temp = str_new_charstr(CHARACTER_VALUE(defStrP));
  rph_tm_register_protect(tm);
  tm_add_alt_mod(tm, temp);
  return R_NilValue;
}


SEXP rph_tm_altmod_set_subst_mod(SEXP tmP, SEXP whichModP, SEXP substModP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  int whichMod = INTEGER_VALUE(whichModP);
  AltSubstMod *altmod;
  if (tm->alt_subst_mods == NULL || lst_size(tm->alt_subst_mods) < whichMod)
    die("ERROR: not enough alt subst  mods (%i %i)\n",
	tm->alt_subst_mods == NULL ? 0 : lst_size(tm->alt_subst_mods),
	whichMod);
  altmod = lst_get_ptr(tm->alt_subst_mods, whichMod-1);
  altmod->subst_mod = tm_get_subst_mod_type(CHARACTER_VALUE(substModP));
  return R_NilValue;
}


SEXP rph_tm_altmod_set_backgd(SEXP tmP, SEXP whichModP, SEXP backgdP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  int whichMod = INTEGER_VALUE(whichModP), i;
  double *doubleP;
  AltSubstMod *altmod;
  if (tm->alt_subst_mods == NULL || lst_size(tm->alt_subst_mods) < whichMod)
    die("ERROR: not enough alt subst  mods (%i %i)\n",
	tm->alt_subst_mods == NULL ? 0 : lst_size(tm->alt_subst_mods),
	whichMod);
  rph_tm_register_protect(tm);
  altmod = lst_get_ptr(tm->alt_subst_mods, whichMod-1);
  if (altmod->backgd_freqs != NULL)
    vec_free(altmod->backgd_freqs);
  if (backgdP == R_NilValue) {
    altmod->backgd_freqs = NULL;
    return R_NilValue;
  }
  if (tm->rate_matrix ==  NULL)
    die("tm->rate_matrix is NULL in rph_tm_altmod_set_backgd\n");
  if (LENGTH(backgdP) != tm->rate_matrix->size)
    die("bad dimensions in rph_tm_altmod_set_backgd");
  altmod->backgd_freqs = vec_new(LENGTH(backgdP));
  PROTECT(backgdP = AS_NUMERIC(backgdP));
  doubleP = NUMERIC_POINTER(backgdP);
  for (i=0; i < LENGTH(backgdP); i++)
    vec_set(altmod->backgd_freqs, i, doubleP[i]);
  UNPROTECT(1);
  return R_NilValue;
}


SEXP rph_tm_altmod_set_ratematrix(SEXP tmP, SEXP whichModP, SEXP matrixP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  int whichMod = INTEGER_VALUE(whichModP), dim;
  AltSubstMod *altmod;
  Matrix *m;
  if (tm->alt_subst_mods == NULL || lst_size(tm->alt_subst_mods) < whichMod)
    die("ERROR: not enough alt subst  mods (%i %i)\n",
	tm->alt_subst_mods == NULL ? 0 : lst_size(tm->alt_subst_mods),
	whichMod);
  altmod = lst_get_ptr(tm->alt_subst_mods, whichMod-1);
  if (altmod->rate_matrix != NULL)
    mm_free(altmod->rate_matrix);
  m = rph_get_matrix(matrixP);
  if (tm->rate_matrix ==  NULL)
    die("ERROR: tm->rate_matrix is NULL in rph_tm_altmodd_set_ratematrix\n");  
  dim = tm->rate_matrix->size;
  if (dim != m->nrows || dim != m->ncols)
    die("Wrong matrix dimensions in rph_tm_altmod_set_ratematrix %i %i %i\n",
	dim, m->nrows, m->ncols);
  rph_tm_register_protect(tm);
  altmod->rate_matrix = mm_new_from_matrix(m, tm->rate_matrix->states, 
					   CONTINUOUS);
  return R_NilValue;
}


SEXP rph_tm_altmod_set_sel_bgc(SEXP tmP, SEXP whichModP, SEXP selP, SEXP bgcP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  int whichMod = INTEGER_VALUE(whichModP);
  AltSubstMod *altmod;

  if (tm->alt_subst_mods == NULL || lst_size(tm->alt_subst_mods) < whichMod)
    die("ERROR: not enough alt subst  mods (%i %i)\n",
	tm->alt_subst_mods == NULL ? 0 : lst_size(tm->alt_subst_mods),
	whichMod);
  altmod = lst_get_ptr(tm->alt_subst_mods, whichMod-1);
  if (selP != R_NilValue) {
    altmod->selection_idx = 0;  //these will get reset later but setting to them something >= 0 
                                //indicates that the parameters are used.
    altmod->selection = NUMERIC_VALUE(selP);
  }
  if (bgcP != R_NilValue) {
    altmod->bgc_idx = 0;
    altmod->bgc = NUMERIC_VALUE(bgcP);
  }
  return R_NilValue;
}


SEXP rph_tree_model_set_matrix(SEXP tmP, SEXP paramsP, SEXP scaleP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  double *params;
  Vector *paramVec=NULL;
  int numparam, paramlen, scale = LOGICAL_VALUE(scaleP);
  if (scale) tm->scale_during_opt = 1;
  if (paramsP == R_NilValue) {
    params = NULL;
    paramlen = 0;
  } else {
    PROTECT(paramsP = AS_NUMERIC(paramsP));
    params = NUMERIC_POINTER(paramsP);
    paramlen = LENGTH(paramsP);
  }
  numparam = tm_get_nratematparams(tm);
  if (numparam != paramlen) {
    if (params != NULL) UNPROTECT(1);
    die("%s requires %i params, got %i\n", tm_get_subst_mod_string(tm->subst_mod),
	numparam, paramlen);
  }
  if (numparam != 0) 
    paramVec = vec_new_from_array(params, numparam);
  tm_set_rate_matrix(tm, paramVec, 0);
  UNPROTECT(1);
  return R_NilValue;  //don't need to return the value since it's the one passed in
}


SEXP rph_tree_model_get_rate_matrix_params(SEXP tmP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  int i, nparam=tm_get_nratematparams(tm);
  double *resultP;
  Vector *v;
  SEXP result;
  if (nparam == 0) return R_NilValue;
  v = vec_new(nparam);
  tm_rate_params_init_from_model(tm, v, 0, tm->selection, 0.0);
  PROTECT(result = NEW_NUMERIC(nparam));
  resultP = NUMERIC_POINTER(result);
  for (i=0; i < nparam; i++)
    resultP[i] = vec_get(v, i);
  UNPROTECT(1);
  return result;
}


SEXP rph_tm_mod_freqs(SEXP tmP, SEXP newBackgdP) {
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  Vector *newfreqs = vec_new(LENGTH(newBackgdP));
  double *doubleP;

  PROTECT(newBackgdP = AS_NUMERIC(newBackgdP));
  doubleP = NUMERIC_POINTER(newBackgdP);
  newfreqs = vec_new_from_array(doubleP, LENGTH(newBackgdP));
  tm_mod_freqs(tm, newfreqs);
  UNPROTECT(1);
  return tmP;
}

SEXP rph_tm_apply_selection_bgc(SEXP matrixP, SEXP alphabetP, SEXP selectionP, 
				SEXP bgcP) {
  double selection=0.0, bgc=0.0;
  MarkovMatrix *mm;
  ListOfLists *lol;
  SEXP result;
  
  mm = mm_new_from_matrix(rph_get_matrix(matrixP), 
			  CHARACTER_VALUE(alphabetP), CONTINUOUS);
  if (selectionP != R_NilValue)
    selection = NUMERIC_VALUE(selectionP);
  if (bgcP != R_NilValue)
    bgc = NUMERIC_VALUE(bgcP);
  tm_apply_selection_bgc(mm, selection, bgc);
  lol = lol_new(1);
  lol_push_matrix(lol, mm->matrix, "rate.matrix");
  PROTECT(result = rph_listOfLists_to_SEXP(lol));
  lol_free(lol);
  mm_free(mm);
  UNPROTECT(1);
  return result;
}


SEXP rph_tm_unapply_selection_bgc(SEXP matrixP, SEXP alphabetP,
				  SEXP selectionP, SEXP bgcP) {
  double selection=0.0, bgc=0.0;
  MarkovMatrix *mm;
  ListOfLists *lol;
  SEXP result;

  mm = mm_new_from_matrix(rph_get_matrix(matrixP), CHARACTER_VALUE(alphabetP), 
			  CONTINUOUS);
  if (selectionP != R_NilValue)
    selection = NUMERIC_VALUE(selectionP);
  if (bgcP != R_NilValue)
    bgc = NUMERIC_VALUE(bgcP);
  tm_unapply_selection_bgc(mm, selection, bgc);
  lol = lol_new(1);
  lol_push_matrix(lol, mm->matrix, "rate.matrix");
  PROTECT(result = rph_listOfLists_to_SEXP(lol));
  mm_free(mm);
  lol_free(lol);
  UNPROTECT(1);
  return result;
}

