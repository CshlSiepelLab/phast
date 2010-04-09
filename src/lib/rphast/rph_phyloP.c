/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_phyloP.c
The RPHAST handles to the phyloP program

Melissa Hubisz
Last updated: 4/8/2010
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <trees.h>
#include <phylo_p.h>

#include <Rdefines.h>

SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);

SEXP rph_phyloP(SEXP modP, SEXP msaP, SEXP methodP, SEXP modeP,
		SEXP gffP, SEXP subtreeP, SEXP wigP, SEXP baseByBaseP,
		SEXP outputFileP) {

  struct phyloP_struct *p = phyloP_struct_new(1);
  SEXP rv;
  char tempstr[1000];

  p->mod = (TreeModel*)EXTPTR_PTR(modP);
  if (msaP != R_NilValue)
    p->msa = (MSA*)EXTPTR_PTR(msaP);
  if (methodP != R_NilValue) {
    strcpy(tempstr, CHARACTER_VALUE(methodP));
    if (strcmp("SPH", tempstr)==0) p->method = SPH;
    else if (strcmp("LRT", tempstr)==0) p->method = LRT;
    else if (strcmp("SCORE", tempstr)==0) p->method = SCORE;
    else if (strcmp("GERP", tempstr)==0) p->method = GERP;
    else 
      die("unknown phyloP method %s\n", tempstr);
  }
  if (modeP != R_NilValue) {
    strcpy(tempstr, CHARACTER_VALUE(modeP));
    if (strcmp("CON", tempstr)==0) p->mode = CON;
    else if (strcmp("ACC", tempstr)==0) p->mode = ACC;
    else if (strcmp("NNEUT", tempstr)==0) p->mode = NNEUT;
    else if (strcmp("CONACC", tempstr)==0) p->mode = CONACC;
    else 
      die("unknown phyloP mode %s\n", tempstr);
  }
  if (gffP != R_NilValue)
    p->feats = (GFF_Set*)EXTPTR_PTR(gffP);
  if (subtreeP != R_NilValue) {
    p->subtree_name = smalloc((1+strlen(CHARACTER_VALUE(subtreeP)))*sizeof(char));
    strcpy(p->subtree_name, CHARACTER_VALUE(subtreeP));
  } 
  if (wigP != R_NilValue) 
    p->output_wig = LOGICAL_VALUE(wigP);
  if (baseByBaseP != R_NilValue)
    p->base_by_base = LOGICAL_VALUE(baseByBaseP);
  if (outputFileP != R_NilValue) 
    p->outfile = fopen_fname(CHARACTER_VALUE(outputFileP), "w");
  
  phyloP(p);
  
  if (p->subtree_name != NULL) free(p->subtree_name);
  if (p->outfile != NULL) fclose(p->outfile);
  PROTECT(rv = rph_listOfLists_to_SEXP(p->results));
  free(p);
  UNPROTECT(1);
  return rv;
}
