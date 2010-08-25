/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_lists.c
The RPHAST handles to functions dealing with phast's
lists structure

Melissa Hubisz
Last updated: 7/24/10
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <lists.h>
#include <Rdefines.h>
#include <misc.h>


void rph_list_free(SEXP listP) {
  List *l;
  l = (List*)EXTPTR_PTR(listP);
  lst_free(l);
}


SEXP rph_list_new_extptr(List *l) {
  SEXP result;
  PROTECT(result=R_MakeExternalPtr((void*)l, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, rph_list_free, 1);
  UNPROTECT(1);
  return result;
}
