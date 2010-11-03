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


void rph_cm_free(SEXP cmP) {
  cm_free((CategoryMap*)EXTPTR_PTR(cmP));
}


SEXP rph_cm_new_extptr(CategoryMap *cm) {
  SEXP result;
  PROTECT(result=R_MakeExternalPtr((void*)cm, R_NilValue, R_NilValue));
  UNPROTECT(1);
  return result;
}


SEXP rph_cm_new_from_gff(SEXP gff) {
  return rph_cm_new_extptr(cm_new_from_features((GFF_Set*)EXTPTR_PTR(gff)));
}


SEXP rph_cm_new_from_str(SEXP str) {
  return rph_cm_new_extptr(cm_new_string_or_file(CHARACTER_VALUE(str)));
}
