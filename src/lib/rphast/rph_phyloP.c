/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
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
#include <getopt.h>
#include <ctype.h>
#include <misc.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <trees.h>
#include <phylo_p.h>
#include <rph_util.h>

#include <Rdefines.h>

SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);

SEXP rph_phyloP(SEXP modP, 
		SEXP msaP, 
		SEXP methodP, 
		SEXP modeP,
		SEXP gffP, 
		SEXP basewiseP, 
		SEXP subtreeP, 
		SEXP branchesP, 
		SEXP refidxP, 
		SEXP outfileP, 
		SEXP outfileOnlyP, 
		SEXP outfileFormatP,
		SEXP priorOnlyP, 
		SEXP nsitesP, 
		SEXP postOnlyP,
		SEXP fitModelP, 
		SEXP epsilonP, 
		SEXP confIntP,
		SEXP quantilesP) {
  struct phyloP_struct *p = phyloP_struct_new(1);
  SEXP rv;
  char tempstr[1000];
  int i, numprotect=0;

  p->mod = (TreeModel*)EXTPTR_PTR(modP);
  rph_tm_register_protect(p->mod);
  if (msaP != R_NilValue) {
    p->msa = (MSA*)EXTPTR_PTR(msaP);
    rph_msa_register_protect(p->msa);
  }
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
  if (gffP != R_NilValue) {
    p->feats = (GFF_Set*)EXTPTR_PTR(gffP);
    rph_gff_register_protect(p->feats);
  }

  if (basewiseP != R_NilValue)
    p->base_by_base = LOGICAL_VALUE(basewiseP);

  if (subtreeP != R_NilValue) 
    p->subtree_name = copy_charstr(CHARACTER_VALUE(subtreeP));

  if (branchesP != R_NilValue) {
    p->branch_name = lst_new_ptr(LENGTH(branchesP));
    for (i=0; i<LENGTH(branchesP); i++)
      lst_push_ptr(p->branch_name, str_new_charstr(CHAR(STRING_ELT(branchesP, i))));
  }
  if (gffP != R_NilValue && refidxP != R_NilValue) 
    p->refidx_feat = INTEGER_VALUE(refidxP);
  else if (refidxP != R_NilValue) {
    p->refidx = INTEGER_VALUE(refidxP);
    if (p->refidx == 0) p->chrom = copy_charstr("align");
  }
  else if (p->msa != NULL)
    p->chrom = copy_charstr(p->msa->names[p->refidx-1]);
  if (outfileP != R_NilValue) {
    p->outfile = fopen_fname(CHARACTER_VALUE(outfileP), "w");
  }
  if (outfileOnlyP != R_NilValue && LOGICAL_VALUE(outfileOnlyP)) {
    lol_free(p->results);
    p->results = NULL;
  }
  if (outfileFormatP != R_NilValue) {
    char *format = copy_charstr(CHARACTER_VALUE(outfileFormatP));
    if (strcmp(format, "default")==0);
    else if (strcmp(format, "wig")==0) {
      p->output_wig = TRUE;
      if (! (p->base_by_base)) die("need basewise scores for wig output");
    }
    else if (strcmp(format, "gff")==0) {
      p->output_gff = TRUE;
      if (p->feats == NULL) die("need features for gff output");
    }
    else die("unknown format string %s", format);
  }

  if (p->method == SPH) {
    if (priorOnlyP != R_NilValue)
      p->prior_only = LOGICAL_VALUE(priorOnlyP);
    if (nsitesP != R_NilValue)
      p->nsites = INTEGER_VALUE(nsitesP);
    if (postOnlyP != R_NilValue)
      p->post_only = LOGICAL_VALUE(postOnlyP);
    if (fitModelP != R_NilValue) {
      p->fit_model = LOGICAL_VALUE(fitModelP);
    }
    if (epsilonP != R_NilValue)
      p->epsilon = NUMERIC_VALUE(epsilonP);
    if (confIntP != R_NilValue)
      p->ci = NUMERIC_VALUE(confIntP);
    if (quantilesP != R_NilValue)
      p->quantiles_only = LOGICAL_VALUE(quantilesP);
  }

  GetRNGstate();
  phyloP(p);
  PutRNGstate();
  
  if (p->outfile != NULL && p->outfile != stdout && p->outfile != stderr) 
    fclose(p->outfile);
  if (p->results != NULL) {
    PROTECT(rv = rph_listOfLists_to_SEXP(p->results));
    numprotect++;
  }
  else rv = R_NilValue;
  fflush(stdout);
  if (numprotect > 0) UNPROTECT(numprotect);
  return rv;
}
