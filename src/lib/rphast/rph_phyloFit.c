/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_phyloFit.c
The RPHAST handles to the phyloFit program

Melissa Hubisz
Last updated: 1/5/2010
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <misc.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <trees.h>
#include <phylo_fit.h>

#include <Rdefines.h>
#include <R_ext/Random.h>

TreeNode* rph_tree_new(SEXP treeStr);
subst_mod_type rph_get_subst_mod(SEXP mod);

struct phyloFit_result_struct {
  List *models;
  List *labels;
  //  List *errors;  //should be one per model if not null
  //...  others to be added
};

void rph_phyloFit_result_struct_free(struct phyloFit_result_struct *pf) {
  int i;
//  printf("rph_phyloFit_result_struct_free\n");
  if (pf->models != NULL) {
    for (i=0; i<lst_size(pf->models); i++)
      tm_free((TreeModel*)lst_get_ptr(pf->models, i));
    lst_free(pf->models);
  }
  if (pf->labels != NULL) {
    lst_free_strings(pf->labels);
    lst_free(pf->labels);
  }
  // need to add more as more types of results are implemented
  free(pf);

}


void rph_phyloFit_result_free(SEXP pfP) {
  rph_phyloFit_result_struct_free((struct phyloFit_result_struct*)EXTPTR_PTR(pfP));
}


SEXP rph_phyloFit_result_new_extptr(struct phyloFit_result_struct *pf) {
  SEXP result;
  PROTECT(result=R_MakeExternalPtr((void*)pf, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, rph_phyloFit_result_free, 1);
  UNPROTECT(1);
  return result;
}



//return the number of models in the result
SEXP rph_phyloFit_result_num_models(SEXP pfResultP) {
  struct phyloFit_result_struct *pfResult;
  SEXP result;
  int *integerP;
  
  pfResult = (struct phyloFit_result_struct*)EXTPTR_PTR(pfResultP);
  PROTECT(result = NEW_INTEGER(1));
  integerP = INTEGER_POINTER(result);
  if (pfResult->models == NULL)
    integerP[0] = 0;
  else integerP[0] = lst_size(pfResult->models);
  UNPROTECT(1);
  return result;
}


SEXP rph_phyloFit_result_has_names(SEXP pfResultP) {
  struct phyloFit_result_struct *pfResult;
  SEXP result;
  int *integerP;
  
  pfResult = (struct phyloFit_result_struct*)EXTPTR_PTR(pfResultP);
  PROTECT(result = NEW_LOGICAL(1));
  integerP = LOGICAL_POINTER(result);
  integerP[0] = (pfResult->labels != NULL);
  UNPROTECT(1);
  return result;
}


SEXP rph_phyloFit_result_get_model(SEXP pfResultP, SEXP whichP) {
  struct phyloFit_result_struct *pfResult;
  int which = INTEGER_VALUE(whichP);
  which--;  //indices are 1-based in R but 0-based in C
  pfResult = (struct phyloFit_result_struct*)EXTPTR_PTR(pfResultP);
  if (pfResult->models == NULL ||
      which >= lst_size(pfResult->models))
    die("internal error, phyloFit_result_model has no %ith element",
	which+1);
  //note: do not register this pointer for cleanup as the
  //phyloFit_result_struct will be cleaned
  return R_MakeExternalPtr(lst_get_ptr(pfResult->models, which),
			   R_NilValue, R_NilValue);
}

SEXP rph_phyloFit_result_get_name(SEXP pfResultP, SEXP whichP) {
  struct phyloFit_result_struct *pfResult;
  int which;
  SEXP result;
  String *label;
  which = INTEGER_VALUE(whichP);
  which--;
  pfResult = (struct phyloFit_result_struct*)EXTPTR_PTR(pfResultP);
  if (pfResult->labels == NULL) return R_NilValue;
  if (which >= lst_size(pfResult->labels))
    die("internal error, phyloFit_result_get_name has no %ith element",
	which+1);
  label = (String*)lst_get_ptr(pfResult->labels, which);
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(label->chars));
  UNPROTECT(1);
  return result;
}


SEXP rph_phyloFit(SEXP msaP, 
		  SEXP treeStrP, 
		  SEXP substModP,
		  SEXP scaleOnlyP,
		  SEXP scaleSubtreeP,
		  SEXP nratesP,
		  SEXP alphaP,
		  SEXP rateConstantsP,
		  SEXP initModP,
		  SEXP noFreqsP,
		  SEXP noRatesP,
		  SEXP initRandomP,
		  SEXP initParsimonyP,
		  SEXP clockP,
		  SEXP emP,
		  SEXP precisionP,
		  SEXP gffP,
		  SEXP ninfSitesP,
		  SEXP quietP) {
  struct phyloFit_struct *pf;
  int numProtect=0, i;
  double *doubleP;
  char *die_message=NULL;
  struct phyloFit_result_struct *result=NULL;

  GetRNGstate(); //seed R's random number generator
  pf = phyloFit_struct_new(1);  //sets appropriate defaults for RPHAST mode

  pf->msa = (MSA*)EXTPTR_PTR(msaP);

  if (treeStrP != R_NilValue) 
    pf->tree = rph_tree_new(treeStrP);

  pf->subst_mod = rph_get_subst_mod(substModP);  //not allowed to be null
  
  pf->estimate_scale_only = LOGICAL_VALUE(scaleOnlyP);
  
  if (scaleSubtreeP != R_NilValue) {
    pf->subtree_name = smalloc((1+strlen(CHARACTER_VALUE(scaleSubtreeP)))*sizeof(char));
    strcpy(pf->subtree_name, CHARACTER_VALUE(scaleSubtreeP));
  }
  
  pf->nratecats = INTEGER_VALUE(nratesP);
  
  if (alphaP != R_NilValue)
    pf->alpha = NUMERIC_VALUE(alphaP);

  if (rateConstantsP != R_NilValue) {
    PROTECT(rateConstantsP = AS_NUMERIC(rateConstantsP));
    numProtect++;
    doubleP = NUMERIC_POINTER(rateConstantsP);
    pf->rate_consts = lst_new_dbl(LENGTH(rateConstantsP));
    for (i=0; i<LENGTH(rateConstantsP); i++)
      lst_push_dbl(pf->rate_consts, doubleP[i]);
  }

  if (initModP != R_NilValue)
    pf->input_mod = (TreeModel*)EXTPTR_PTR(initModP);

  pf->random_init = LOGICAL_VALUE(initRandomP);

  pf->no_freqs = LOGICAL_VALUE(noFreqsP);

  pf->no_rates = LOGICAL_VALUE(noRatesP);
  
  pf->init_parsimony = LOGICAL_VALUE(initParsimonyP);
  
  pf->assume_clock = LOGICAL_VALUE(clockP);

  pf->use_em = LOGICAL_VALUE(emP);

  if (strcmp(CHARACTER_VALUE(precisionP), "LOW")==0)
    pf->precision = OPT_LOW_PREC;
  else if (strcmp(CHARACTER_VALUE(precisionP), "MED")==0)
    pf->precision = OPT_MED_PREC;
  else if (strcmp(CHARACTER_VALUE(precisionP), "HIGH")==0)
    pf->precision = OPT_HIGH_PREC;
  else {
    die_message = "invalid precision";
    goto rph_phyloFit_end;
  }

  pf->estimated_models = lst_new_ptr(50);
  pf->model_labels = lst_new_ptr(50);
  if (gffP != R_NilValue)
    pf->gff = (GFF_Set*)EXTPTR_PTR(gffP);

  if (ninfSitesP != R_NilValue)
    pf->nsites_threshold = INTEGER_VALUE(ninfSitesP);
  
  pf->quiet = LOGICAL_VALUE(quietP);

  run_phyloFit(pf);

  if (pf->estimated_models == NULL)
    die("ERROR in rph_phyloFit: estimated_models is NULL\n");
  result = smalloc(sizeof(struct phyloFit_result_struct));
  result->models = pf->estimated_models;
  result->labels = pf->model_labels;

 rph_phyloFit_end:
  if (pf->tree != NULL)
    tr_free(pf->tree);
  if (pf->subtree_name != NULL)
    free(pf->subtree_name);
  if (pf->rate_consts != NULL)
    lst_free(pf->rate_consts);
  free(pf);
  PutRNGstate();
  if (numProtect > 0) 
    UNPROTECT(numProtect);
  if (die_message != NULL) die(die_message);
  return rph_phyloFit_result_new_extptr(result);
}
