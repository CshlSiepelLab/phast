/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_wig.c
The RPHAST handles to functions dealing with wigs from
the phast package.

Melissa Hubisz
Last updated: 4/10/2012
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <misc.h>
#include <gff.h>
#include <list_of_lists.h>
#include <rph_util.h>
#include <wig.h>

#include <Rdefines.h>


SEXP rph_wig_read(SEXP filename) {
  FILE *infile = phast_fopen(CHARACTER_VALUE(filename), "r");
  SEXP rv;
  PROTECT(rv = rph_gff_new_extptr(gff_read_wig(infile)));
  phast_fclose(infile);
  UNPROTECT(1);
  return rv;
}

SEXP rph_wig_print(SEXP gffP, SEXP filename, SEXP append) {
  FILE *outfile;
  char *mode;
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  gff_register_protect(gff);
  if (LOGICAL_VALUE(append)==TRUE) mode="a";
  else mode="w";
  if (filename == R_NilValue)
    outfile = stdout;
  else outfile = phast_fopen(CHARACTER_VALUE(filename), mode);
  
  wig_print(outfile, gff);
  if (outfile != stdout) phast_fclose(outfile);
  return R_NilValue;
}

