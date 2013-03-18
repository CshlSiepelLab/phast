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
#include <getopt.h>
#include <ctype.h>
#include <misc.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <trees.h>
#include <phylo_fit.h>
#include <rph_util.h>

#include <Rdefines.h>
#include <R_ext/Random.h>

TreeNode* rph_tree_new(SEXP treeStr);
subst_mod_type rph_get_subst_mod(SEXP mod);


SEXP rph_phyloFit(SEXP msaP, 
		  SEXP treeStrP, 
		  SEXP substModP,
		  SEXP scaleOnlyP,
		  SEXP scaleSubtreeP,
		  SEXP nratesP,
		  SEXP alphaP,
		  SEXP rateConstantsP,
		  SEXP initModP,
		  SEXP initBackgdFromDataP,
		  SEXP initRandomP,
		  SEXP initParsimonyP,
		  SEXP clockP,
		  SEXP emP,
		  SEXP maxEmItsP,
		  SEXP precisionP,
		  SEXP gffP,
		  SEXP ninfSitesP,
		  SEXP quietP,
		  SEXP noOptP,
		  SEXP boundP,
		  SEXP logFileP,
		  SEXP selectionP) {
  struct phyloFit_struct *pf;
  int numProtect=0, i;
  double *doubleP;
  char *die_message=NULL;
  SEXP rv=R_NilValue;
  List *new_rate_consts = NULL;
  List *new_rate_weights = NULL;

  GetRNGstate(); //seed R's random number generator
  pf = phyloFit_struct_new(1);  //sets appropriate defaults for RPHAST mode

  pf->msa = (MSA*)EXTPTR_PTR(msaP);

  if (treeStrP != R_NilValue) 
    pf->tree = rph_tree_new(treeStrP);

  pf->use_em = LOGICAL_VALUE(emP);

  if (rateConstantsP != R_NilValue) {
    PROTECT(rateConstantsP = AS_NUMERIC(rateConstantsP));
    numProtect++;
    doubleP = NUMERIC_POINTER(rateConstantsP);
    new_rate_consts = lst_new_dbl(LENGTH(rateConstantsP));
    for (i=0; i < LENGTH(rateConstantsP); i++)
      lst_push_dbl(new_rate_consts, doubleP[i]);
//    pf->use_em = 1;
  }

  if (initModP != R_NilValue) {
    pf->input_mod = (TreeModel*)EXTPTR_PTR(initModP);
    pf->subst_mod = pf->input_mod->subst_mod;
    tm_register_protect(pf->input_mod);
    
    if (new_rate_consts == NULL && pf->input_mod->rK != NULL && pf->input_mod->nratecats > 1) {
      new_rate_consts = lst_new_dbl(pf->input_mod->nratecats);
      for (i=0; i < pf->input_mod->nratecats; i++) 
	lst_push_dbl(new_rate_consts, pf->input_mod->rK[i]);
//      pf-> = 1;
    }

    if (pf->input_mod->empirical_rates && pf->input_mod->freqK != NULL && pf->input_mod->nratecats > 1) {
      new_rate_weights = lst_new_dbl(pf->input_mod->nratecats);
      for (i=0; i < pf->input_mod->nratecats; i++)
	lst_push_dbl(new_rate_weights, pf->input_mod->freqK[i]);
    }

    tm_reinit(pf->input_mod, 
	      rph_get_subst_mod(substModP),
	      nratesP == R_NilValue ? pf->input_mod->nratecats : INTEGER_VALUE(nratesP),
	      NUMERIC_VALUE(alphaP),
	      new_rate_consts,
	      new_rate_weights);
  } else {
    if (nratesP != R_NilValue)
      pf->nratecats = INTEGER_VALUE(nratesP);
    if (alphaP != R_NilValue)
      pf->alpha = NUMERIC_VALUE(alphaP);
    if (rateConstantsP != R_NilValue) {
      pf->rate_consts = new_rate_consts;
      if (nratesP == R_NilValue)
	pf->nratecats = lst_size(new_rate_consts);
      else if (lst_size(new_rate_consts) != pf->nratecats) 
	die("length of new_rate_consts does not match nratecats\n");
    }
  }
  pf->subst_mod = rph_get_subst_mod(substModP);
  
  pf->estimate_scale_only = LOGICAL_VALUE(scaleOnlyP);
  
  if (scaleSubtreeP != R_NilValue) {
    pf->subtree_name = smalloc((1+strlen(CHARACTER_VALUE(scaleSubtreeP)))*sizeof(char));
    strcpy(pf->subtree_name, CHARACTER_VALUE(scaleSubtreeP));
  }
  
  pf->random_init = LOGICAL_VALUE(initRandomP);

  pf->init_backgd_from_data = LOGICAL_VALUE(initBackgdFromDataP);
  
  pf->init_parsimony = LOGICAL_VALUE(initParsimonyP);
  
  pf->assume_clock = LOGICAL_VALUE(clockP);

  if (maxEmItsP != R_NilValue)
    pf->max_em_its = INTEGER_VALUE(maxEmItsP);

  pf->precision = get_precision(CHARACTER_VALUE(precisionP));
  if (pf->precision == OPT_UNKNOWN_PREC) {
    die_message = "invalid precision";
    goto rph_phyloFit_end;
  }

  if (gffP != R_NilValue) {
    pf->gff = (GFF_Set*)EXTPTR_PTR(gffP);
    gff_register_protect(pf->gff);
  }

  if (ninfSitesP != R_NilValue)
    pf->nsites_threshold = INTEGER_VALUE(ninfSitesP);
  
  pf->quiet = LOGICAL_VALUE(quietP);

  if (noOptP != R_NilValue) {
    int len=LENGTH(noOptP), pos=0;
    char *temp;
    for (i=0; i < LENGTH(noOptP); i++) 
      len += strlen(CHARACTER_VALUE(STRING_ELT(noOptP, i)));
    temp = smalloc(len*sizeof(char));
    for (i=0; i < LENGTH(noOptP); i++) {
      if (i != 0) temp[pos++] = ',';
      sprintf(&temp[pos], "%s", CHARACTER_VALUE(STRING_ELT(noOptP, i)));
      pos += strlen(CHARACTER_VALUE(STRING_ELT(noOptP, i)));
    }
    if (pos != len-1) die("ERROR parsing noOpt len=%i pos=%i\n", len, pos);
    temp[pos] = '\0';
    pf->nooptstr = str_new_charstr(temp);
  }

  if (boundP != R_NilValue) {
    pf->bound_arg = lst_new_ptr(LENGTH(boundP));
    for (i=0; i < LENGTH(boundP); i++) {
      String *temp = str_new_charstr(CHARACTER_VALUE(STRING_ELT(boundP, i)));
      lst_push_ptr(pf->bound_arg, temp);
    }
  }

  if (logFileP != R_NilValue) {
    if (IS_CHARACTER(logFileP)) 
      pf->logf = phast_fopen(CHARACTER_VALUE(logFileP), "w+");
    else if (IS_LOGICAL(logFileP) &&
	     LOGICAL_VALUE(logFileP)) {
      pf->logf = stdout;
    }
  }

  if (selectionP != R_NilValue) {
    pf->use_selection = TRUE;
    pf->selection = NUMERIC_VALUE(selectionP);
  }

  msa_register_protect(pf->msa);

  run_phyloFit(pf);
  rv = PROTECT(rph_listOfLists_to_SEXP(pf->results));
  numProtect++;

 rph_phyloFit_end:
  if (pf->logf != NULL && pf->logf != stdout && pf->logf != stderr)
    phast_fclose(pf->logf);
  PutRNGstate();
  if (die_message != NULL) die(die_message);
  if (numProtect > 0) 
    UNPROTECT(numProtect);
  return rv;
}
