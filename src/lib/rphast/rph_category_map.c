/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_category_map.c
The RPHAST handles to functions dealing with
category maps in the PHAST package.

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
#include <category_map.h>

#include <Rdefines.h>
#include <rph_util.h>

void rph_cm_protect(CategoryMap *cm) {
  int i;
  if (cm == NULL) return;
  rph_mem_protect(cm);
  for (i=0; i <= cm->ncats; i++) {
    int len=0;
    if (cm->ranges[i] != NULL) {
      len = cm->ranges[i]->end_cat_no - cm->ranges[i]->start_cat_no;
      rph_lst_protect(cm->ranges[i]->feature_types);
      rph_mem_protect(cm->ranges[i]);
      rph_lst_protect(cm->conditioned_on[i]);
      i += len;
    }
  }
  rph_mem_protect(cm->ranges);
  rph_mem_protect(cm->conditioned_on);
  rph_mem_protect(cm->labelling_precedence);
  rph_mem_protect(cm->fill_precedence);
  if (cm->unspooler != NULL) 
    rph_mem_protect(cm->unspooler);
}


void rph_cm_free(SEXP cmP) {
  CategoryMap *cm = (CategoryMap*)EXTPTR_PTR(cmP);
  rph_unregister_protected(cm);
  cm_free(cm);
}


void rph_cm_register_protect(CategoryMap *cm) {
  rph_register_protected_object(cm, (void (*)(void *))rph_cm_protect);
}


SEXP rph_cm_new_extptr(CategoryMap *cm) {
  SEXP result;
  rph_cm_register_protect(cm);
  PROTECT(result=R_MakeExternalPtr((void*)cm, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, rph_cm_free, 1);
  UNPROTECT(1);
  return result;
}


SEXP rph_cm_new_from_gff(SEXP gff) {
  return rph_cm_new_extptr(cm_new_from_features((GFF_Set*)EXTPTR_PTR(gff)));
}


SEXP rph_cm_new_from_str(SEXP str) {
  return rph_cm_new_extptr(cm_new_string_or_file(CHARACTER_VALUE(str)));
}
