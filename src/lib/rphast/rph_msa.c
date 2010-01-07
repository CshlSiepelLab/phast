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
  SEXP result;

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
    alphabet = malloc((strlen(CHARACTER_VALUE(alphabetP))+1)*sizeof(char));
    strcpy(alphabet, CHARACTER_VALUE(alphabetP));
  }

  /*  printf("nseqs=%i length=%i\n", nseqs, length);
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
  else printf("alphabet is NULL\n");*/
  msa = msa_new(seqs, names, nseqs, length, alphabet);

  //don't free seqs or names because they are used by MSA object
  if (alphabet != NULL) free(alphabet);
  
  PROTECT(result = R_MakeExternalPtr((void*)msa, R_NilValue, R_NilValue));
  numProtect++;

  UNPROTECT(numProtect);
  return result;
}


SEXP rph_msa_read(SEXP filenameP, SEXP formatP, SEXP gffP, SEXP do4dP,
		  SEXP alphabetP, SEXP tupleSizeP, SEXP refseqP) {
  int do4d, tupleSize=1;
  GFF_Set *gff=NULL;
  msa_format_type fmt;
  MSA *msa;
  FILE *infile;
  char *alphabet = NULL, *refseq=NULL;

  do4d = INTEGER_VALUE(do4dP);
  if (gffP != R_NilValue) 
    gff = (GFF_Set*)EXTPTR_PTR(gffP);
  fmt = msa_str_to_format(CHARACTER_VALUE(formatP));
  if ((int)fmt == -1) 
    fmt = FASTA;
  if (alphabetP != R_NilValue) {
    alphabet = malloc((strlen(CHARACTER_VALUE(alphabetP))+1)*sizeof(char));
    strcpy(alphabet, CHARACTER_VALUE(alphabetP));
  }
  if (refseqP != R_NilValue) {
    refseq = malloc((strlen(CHARACTER_VALUE(refseqP))+1)*sizeof(char));
    strcpy(refseq, CHARACTER_VALUE(refseqP));
  }
  if (tupleSizeP != R_NilValue)
    tupleSize = INTEGER_VALUE(tupleSizeP);
  infile = fopen_fname(CHARACTER_VALUE(filenameP), "r");
  /*  if (fmt == MAF) {
    msa = maf_read(infile, refseq, tupleSize, alphabet, gff, NULL, 
		   } else {
    msa = msa_new_from_file(infile, fmt, alphabet);
    }*/
  if (alphabet != NULL)
    free(alphabet);
  if (refseq != NULL)
    free(refseq);
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


SEXP rph_msa_valid_fmt_str(SEXP formatP) {
  msa_format_type fmt;
  SEXP result;
  int *resultP;
  PROTECT( result = allocVector(LGLSXP, 1));
  resultP = LOGICAL_POINTER(result);
  fmt = msa_str_to_format(CHARACTER_VALUE(formatP));
  resultP[0] = (fmt != -1);
  UNPROTECT(1);
  return result;
}


//print sequence in given format.  If format not valid,
//use FASTA, but give no warning.
SEXP rph_msa_printSeq(SEXP msaP, SEXP filenameP, SEXP formatP, 
		      SEXP prettyPrintP) {
  MSA *msa;
  msa_format_type fmt;
  msa = (MSA*)EXTPTR_PTR(msaP);
  printf("got msa pointer\n");
  fmt = msa_str_to_format(CHARACTER_VALUE(formatP));
  printf("got format string %i\n", (int)fmt);
  if ((int)fmt == -1) 
    fmt = FASTA;
  if (filenameP != R_NilValue)
    msa_print_to_file(CHARACTER_VALUE(filenameP), 
		      msa, fmt, INTEGER_VALUE(prettyPrintP));
  else msa_print(stdout, msa, fmt, INTEGER_VALUE(prettyPrintP));
  return R_NilValue;
}



/****** Accessor functions**********/

/* Returns a SEXP containing the sequence.
 */
SEXP rph_msa_seqs(SEXP msaP) {
  SEXP result;
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  char *tempseq;
  int seq;
  
  PROTECT(result = NEW_CHARACTER(msa->nseqs));

  if (msa->ss != NULL) {
    /* rather than calling ss_to_msa, which will allocate an entire new
       sequence that we don't need, only convert one species at a time.
     */
    for  (seq = 0; seq < msa->nseqs; seq++) {
      tempseq = ss_get_one_seq(msa, seq);
      SET_STRING_ELT(result, seq, mkChar(tempseq));
      free(tempseq);
    }
  } else {
    for (seq = 0; seq < msa->nseqs; seq++) {
      SET_STRING_ELT(result, seq, mkChar(msa->seqs[seq]));
    }
  }
  UNPROTECT(1);
  return result;
}


SEXP rph_msa_seqlen(SEXP msaP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP result;
  int *resultP;
  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  resultP[0] = msa->length;
  UNPROTECT(1);
  return result;
}


SEXP rph_msa_nseq(SEXP msaP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP result;
  int *resultP;
  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  resultP[0] = msa->nseqs;
  UNPROTECT(1);
  return result;
}


SEXP rph_msa_seqNames(SEXP msaP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP result;
  int i;
  if (msa->names==NULL) return R_NilValue;
  PROTECT(result = NEW_CHARACTER(msa->nseqs));
  //PROTECT(mychar = allocVector(STRSXP, 5));
  for (i=0; i<msa->nseqs; i++) 
    SET_STRING_ELT(result, i, mkChar(msa->names[i]));
  UNPROTECT(1);
  return result;
}

SEXP rph_msa_alphabet(SEXP msaP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP result;
  if (msa->alphabet==NULL) return R_NilValue;
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(msa->alphabet));
  UNPROTECT(1);
  return result;
}

/* TODO (maybe not here but in general): at some point it would
   be nice to have unordered MSA objects without requiring SS format */
SEXP rph_msa_isOrdered(SEXP msaP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP result;
  int *resultP;
  if (msa->ss == NULL) return R_NilValue;

  PROTECT(result = NEW_LOGICAL(1));

  resultP = LOGICAL_POINTER(result);
  resultP[0] = (msa->ss->tuple_idx != NULL);
  UNPROTECT(1);
  return result;
}

SEXP rph_msa_idxOffset(SEXP msaP, SEXP tagP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP result;
  int *resultP;
  
  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  resultP[0] = msa->idx_offset;
  UNPROTECT(1);
  return result;
}
