/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_subst_mod.c
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
#include <subst_mods.h>
#include <external_libs.h>
#include <misc.h>

#include <Rdefines.h>


subst_mod_type rph_get_subst_mod(SEXP mod) {
  return tm_get_subst_mod_type(CHARACTER_VALUE(mod));
}


SEXP rph_subst_mods_is_valid_string(SEXP mod) {
  subst_mod_type subst_mod = tm_get_subst_mod_type(CHARACTER_VALUE(mod));
  SEXP result;
  int *resultP;
  
  PROTECT(result = NEW_LOGICAL(1));
  resultP = LOGICAL_POINTER(result);
  resultP[0] = (subst_mod != UNDEF_MOD);
  UNPROTECT(1);
  return result;
}


SEXP rph_subst_mods_list_all(SEXP nilvalue) {
  SEXP result;
  int i, total=0;
  if (nilvalue != R_NilValue)
    die("rph_subst_mods_list_all doesn't really take an argument");
      
  for (i=0; ; i++) {
    if ((subst_mod_type)i == UNDEF_MOD) break;
    total++;
  }
  PROTECT(result = NEW_CHARACTER(total));
  for (i=0; i<total; i++) 
    SET_STRING_ELT(result, i, mkChar(tm_get_subst_mod_string((subst_mod_type)i)));
  UNPROTECT(1);
  return result;
}
  

