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


/* Note: this changes the object passed in */
SEXP rph_msa_reduce_to_4d(SEXP msa) {
  CategoryMap  *cm;
  cm = cm_new_string_or_file("NCATS=3; CDS 1-3");
  reduce_to_4d((MSA*)EXTPTR_PTR(msa), cm);
  cm_free(cm);
  return msa;
}


SEXP rph_msa_read(SEXP filenameP, SEXP formatP, SEXP gffP, 
		  SEXP do4dP, SEXP alphabetP, 
		  SEXP tupleSizeP, SEXP refseqP, SEXP orderedP,
		  SEXP docatsP) {
  int do4d, tupleSize=1, ordered, i;
  GFF_Set *gff=NULL;
  msa_format_type fmt;
  CategoryMap *cm=NULL;
  MSA *msa;
  FILE *infile, *refseq=NULL;
  char *alphabet = NULL, *reverse_groups_tag=NULL;
  int numProtect=0;
  List *cats_to_do=NULL, *cats_to_do_str;

  fmt = msa_str_to_format(CHARACTER_VALUE(formatP));
  if ((int)fmt == -1) 
    fmt = FASTA;
  if (alphabetP != R_NilValue) {
    alphabet = malloc((strlen(CHARACTER_VALUE(alphabetP))+1)*sizeof(char));
    strcpy(alphabet, CHARACTER_VALUE(alphabetP));
  }
  if (refseqP != R_NilValue) 
    refseq = fopen_fname(CHARACTER_VALUE(refseqP), "r");
  if (tupleSizeP != R_NilValue)
    tupleSize = INTEGER_VALUE(tupleSizeP);
  if (gffP != R_NilValue) 
    gff = (GFF_Set*)EXTPTR_PTR(gffP);
  ordered = INTEGER_VALUE(orderedP);
  do4d = INTEGER_VALUE(do4dP);
  if (do4d) {
    reverse_groups_tag = "transcript_id";
    cm = cm_new_string_or_file("NCATS=3; CDS 1-3");
    ordered=0;
    tupleSize=3;
  }
  if (gff != NULL && cm==NULL) 
    cm = cm_new_from_features(gff);
  if (docatsP != R_NilValue) {
    int numcats;
    PROTECT(docatsP = AS_CHARACTER(docatsP));
    numProtect++;
    numcats = LENGTH(docatsP);
    cats_to_do_str = lst_new_ptr(numcats);
    for (i=0; i<numcats; i++) {
      lst_push_ptr(cats_to_do_str, 
		   str_new_charstr(CHAR(STRING_ELT(docatsP, i))));
      if (cm_get_category(cm, lst_get_ptr(cats_to_do_str, i))==0)
	Rf_warning("category %s not found in GFF", 
		   (CHAR(STRING_ELT(docatsP, i))));
    }
    cats_to_do = cm_get_category_list(cm, cats_to_do_str, TRUE);
    lst_free_strings(cats_to_do_str);
    lst_free(cats_to_do_str);
  } else if (gff != NULL) {
    cats_to_do = lst_new_int(cm->ncats);
    for (i=1; i<=cm->ncats; i++)
      lst_push_int(cats_to_do, i);
  }

  infile = fopen_fname(CHARACTER_VALUE(filenameP), "r");
  if (fmt == MAF) {  //reads and automatically converts to SS format
    printf("cats_to do is nULL=%i\n", cats_to_do==NULL);
    msa = maf_read_cats(infile, refseq, tupleSize, alphabet, gff, cm, -1, 
			ordered, reverse_groups_tag, NO_STRIP, 0, 
			cats_to_do);
  } else {
    msa = msa_new_from_file(infile, fmt, alphabet);
    if (gff != NULL) {
      if (msa->ss != NULL && msa->ss->tuple_idx == NULL)
	die("ordered representation of alignment required to extract features");
      
      /* convert GFF to coordinate frame of alignment */
      if (msa->idx_offset !=0) {
	for (i=0; i<lst_size(gff->features); i++) {
	  GFF_Feature *f = lst_get_ptr(gff->features, i);
	  f->start -= msa->idx_offset;
	  f->end -= msa->idx_offset;
	}
      }
      msa_map_gff_coords(msa, gff, -1, 0, 0, cm);
      if (reverse_groups_tag != NULL) { /*reverse complement by group */
	if (fmt == SS) 
	  die("ERROR: need an explicit representation of the alignment to reverse complement");
	gff_group(gff, reverse_groups_tag);
	msa_reverse_compl_feats(msa, gff, NULL);
      }
      msa_label_categories(msa, gff, cm);
    }
  }
  if (msa == NULL) die("ERROR reading %s\n", CHARACTER_VALUE(filenameP));
  if (gff != NULL && fmt != MAF) {
    if ((fmt == SS || fmt==MAF) && (msa->ss==NULL || msa->ss->tuple_idx==NULL)) {
      die("ordered representation of alignment required to extract features");
    }
  }
  if (tupleSize != 1 || cats_to_do != NULL) {
    if (msa->ss == NULL)
      ss_from_msas(msa, tupleSize, ordered, cats_to_do, NULL, NULL, -1);
    else {
      if (msa->ss->tuple_size < tupleSize)
	die("ERROR: input tuple size must be at least as large as tupleSize");
      if (msa->ss->tuple_idx != NULL && ordered==0) {
	free(msa->ss->tuple_idx);
	free(msa->ss->tuple_idx);
      }
      if (msa->ss->tuple_size > tupleSize)
	ss_reduce_tuple_size(msa, tupleSize);
    }
  }
  if (do4d) reduce_to_4d(msa, cm);

  if (cats_to_do != NULL)
    lst_free(cats_to_do);
  if (cm != NULL)
    cm_free(cm);
  if (alphabet != NULL)
    free(alphabet);
  if (refseq != NULL)
    fclose(refseq);
  fclose(infile);
  if (numProtect > 0)
    UNPROTECT(numProtect);
  return R_MakeExternalPtr((void*)msa, R_NilValue, R_NilValue);
}


SEXP rph_msa_free(SEXP msaP) {
  MSA *msa;
  msa = (MSA*)EXTPTR_PTR(msaP);
  msa_free(msa);
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
  fmt = msa_str_to_format(CHARACTER_VALUE(formatP));
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
      printf("calling set_string_elt %s\n", msa->seqs[seq]);
      SET_STRING_ELT(result, seq, mkChar(msa->seqs[seq]));
    }
  }
  UNPROTECT(1);
  return result;
}


SEXP rph_msa_seqlen(SEXP msaP, SEXP refseqP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP result;
  int *resultP, seqidx;
  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  if (refseqP == R_NilValue)
    resultP[0] = msa->length;
  else {
    seqidx = msa_get_seq_idx(msa, CHARACTER_VALUE(refseqP));
    if (seqidx == -1) 
      die("sequence %s not found", CHARACTER_VALUE(refseqP));
    resultP[0] = msa_seqlen(msa, seqidx);
  }
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


SEXP rph_msa_sub_alignment(SEXP msaP, SEXP seqsP, SEXP keepP, 
			   SEXP startcolP, SEXP endcolP, 
			   SEXP refseqNameP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP), *subMsa;
  int numseq, i, keep=1, startcol, endcol, numProtect=0, refseq;
  List *seqlist_str, *l=NULL;
  msa_coord_map *map = NULL;

  if (seqsP != R_NilValue) {
    PROTECT(seqsP = AS_CHARACTER(seqsP));
    numProtect++;
    numseq = LENGTH(seqsP);
    seqlist_str = lst_new_ptr(numseq);
    for (i=0; i<numseq; i++)  
      lst_push_ptr(seqlist_str, str_new_charstr(CHAR(STRING_ELT(seqsP, i))));
    l = msa_seq_indices(msa, seqlist_str);
    lst_free_strings(seqlist_str);
    lst_free(seqlist_str);
    keep = LOGICAL_VALUE(keepP);
  }

  if (startcolP == R_NilValue)
    startcol = 1;
  else startcol = INTEGER_VALUE(startcolP);
  
  if (endcolP == R_NilValue)
    endcol = msa->length;
  else endcol = INTEGER_VALUE(endcolP);

  if (refseqNameP != R_NilValue) {
    for (refseq=0; refseq<msa->nseqs; refseq++) 
      if (strcmp(CHARACTER_VALUE(refseqNameP), msa->names[refseq])==0) break;
    if (refseq==msa->nseqs) 
      die("no sequences named %s", CHARACTER_VALUE(refseqNameP));      
    if (msa->ss != NULL && msa->ss->tuple_idx == NULL)
      die("an ordered representation of the alignment is required");
    refseq++;
    map = msa_build_coord_map(msa, refseq);
    startcol = msa_map_seq_to_msa(map, startcol);
    if (endcolP != R_NilValue)
      endcol = msa_map_seq_to_msa(map, endcol);
  }

  if (startcol < 0 || startcol > msa->length)
    die("start column out of range");
  if (endcol < 0 || endcol > msa->length)
    die("end column out of range");
    

  subMsa = msa_sub_alignment(msa, l, keep, startcol-1, endcol);
  assert(subMsa != NULL);
  if (l != NULL) lst_free(l);
  if (map != NULL) msa_map_free(map);
  if (numProtect > 0) UNPROTECT(numProtect);
  return R_MakeExternalPtr(subMsa, R_NilValue, R_NilValue);
}
