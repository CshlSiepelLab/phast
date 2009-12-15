/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_msa.c
The RPHAST handles to functions dealing with multiple
sequence alignment functions from the phast package.

Melissa Hubisz
Last updated: 12/14/08
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
#include <maf.h>
#include <rph_util.h>

#include <Rdefines.h>

SEXP rph_msa_new(SEXP seqsP, SEXP namesP, SEXP nseqsP, SEXP lengthP, SEXP alphabetP) {
  char **seqs=NULL, **names=NULL, *alphabet=NULL;
  int nseqs=0, length=0, i, numProtect=0;
  MSA *msa=NULL;

  nseqs = INTEGER_VALUE(nseqsP);
  length = INTEGER_VALUE(lengthP);

  if (namesP != R_NilValue) {
    PROTECT(namesP = AS_CHARACTER(namesP));
    numProtect++;
    names = malloc(nseqs*sizeof(char*));
    for (i=0; i<nseqs; i++) {
      names[i] = malloc((strlen(CHAR(STRING_ELT(namesP, i)))+1)*sizeof(char));
      strcpy(names[i], CHAR(STRING_ELT(namesP, i)));
    }
  }
  if (seqsP != R_NilValue) {
    PROTECT(seqsP = AS_CHARACTER(seqsP));
    numProtect++;
    seqs = malloc(nseqs*sizeof(char*));
    for (i=0; i<nseqs; i++) {
      seqs[i] = malloc((strlen(CHAR(STRING_ELT(seqsP, i)))+1)*sizeof(char));
      strcpy(seqs[i], CHAR(STRING_ELT(seqsP, i)));
    }
  }
  if (alphabetP != R_NilValue) {
    PROTECT(alphabetP = AS_CHARACTER(alphabetP));
    numProtect++;
    alphabet = malloc((strlen(CHAR(STRING_ELT(alphabetP, 0)))+1)*sizeof(char));
    strcpy(alphabet, CHAR(STRING_ELT(alphabetP, 0)));
  }

  printf("nseqs=%i length=%i\n", nseqs, length);
  if (names != NULL) {
    for (i=0; i<nseqs; i++) 
      printf("names[%i]=%s\n", i, names[i]==NULL ? "NULL" : names[i]);
  }
  else printf("names is NULL\n");
  if (seqs != NULL) {
    for (i=0; i<nseqs; i++)
      printf("seqs[%i]=%s\n", i, seqs[i]==NULL ? "NULL" : seqs[i]);
  } else printf("seqs is NULL\n");
  if (alphabet != NULL) printf("alphabet=%s\n", alphabet);
  else printf("alphabet is NULL\n");
  msa = msa_new(seqs, names, nseqs, length, alphabet);

  //don't free seqs or names because they are used by MSA object
  if (alphabet != NULL) free(alphabet);
  UNPROTECT(numProtect);
  return R_MakeExternalPtr((void*)msa, R_NilValue, R_NilValue);
}


SEXP rph_msa_free(SEXP msaP) {
  MSA *msa;
  printf("freeing msa\n");
  msa = (MSA*)EXTPTR_PTR(msaP);
  printf("msa->nseqs=%i\n", msa->nseqs);
  msa_free(msa);
  printf("done rph_msa_free\n");
  return R_NilValue;
}


SEXP rph_msa_valid_fmt(SEXP formatP) {
  msa_format_type fmt;
  SEXP result;
  int *resultP;
  PROTECT( result = allocVector(LGLSXP, 1));
  resultP = LOGICAL_POINTER(result);
  fmt = msa_str_to_format((char*)CHAR(STRING_ELT(formatP, 0)));
  resultP[0] = (fmt != -1);
  UNPROTECT(1);
  return result;
}


//print sequence in given format.  If format not valid,
//use FASTA, but give no warning.
SEXP rph_msa_printSeq(SEXP msaP, SEXP filenameP, SEXP formatP, SEXP prettyPrintP) {
  MSA *msa;
  msa_format_type fmt;
  msa = (MSA*)EXTPTR_PTR(msaP);
  fmt = msa_str_to_format((char*)CHAR(STRING_ELT(formatP, 0)));
  if ((int)fmt == -1) 
    fmt = FASTA;
  if (filenameP != R_NilValue)
    msa_print_to_file((char*)CHAR(STRING_ELT(filenameP, 0)),
		      msa, fmt, INTEGER_VALUE(prettyPrintP));
  else msa_print(stdout, msa, fmt, INTEGER_VALUE(prettyPrintP));
  return R_NilValue;
}

