/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
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
#include <getopt.h>
#include <ctype.h>
#include <misc.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <maf.h>
#include <tree_likelihoods.h>
#include <hmm.h>
#include <category_map.h>
#include <gff.h>
#include <list_of_lists.h>
#include <matrix.h>
#include <sufficient_stats.h>
#include <phylo_fit.h>
#include <rph_util.h>
#include <memory_handler.h>

#include <Rdefines.h>
#include <R_ext/Random.h>

void rph_msa_free(SEXP msaP) {
  MSA *msa;
  msa = (MSA*)EXTPTR_PTR(msaP);
  phast_unregister_protected(msa);
  msa_free(msa);
}


SEXP rph_msa_new_extptr(MSA *msa) {
  SEXP result;
  msa_register_protect(msa);
  PROTECT(result=R_MakeExternalPtr((void*)msa, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, rph_msa_free, 1);
  UNPROTECT(1);
  return result;
}


SEXP rph_msa_new(SEXP seqsP, SEXP namesP, SEXP nseqsP, SEXP lengthP, 
		 SEXP alphabetP, SEXP orderedP, SEXP idxOffsetP) {
  char **seqs=NULL, **names=NULL, *alphabet=NULL;
  int nseqs=0, length=0, i, numProtect=0, ordered;
  MSA *msa=NULL;
  SEXP result;

  nseqs = INTEGER_VALUE(nseqsP);
  length = INTEGER_VALUE(lengthP);
  ordered = LOGICAL_VALUE(orderedP);

  if (namesP != R_NilValue) {
    PROTECT(namesP = AS_CHARACTER(namesP));
    numProtect++;
    names = smalloc(nseqs*sizeof(char*));
    for (i=0; i<nseqs; i++) {
      names[i] = smalloc((strlen(CHAR(STRING_ELT(namesP, i)))+1)*sizeof(char));
      strcpy(names[i], CHAR(STRING_ELT(namesP, i)));
    }
  }
  if (seqsP != R_NilValue) {
    PROTECT(seqsP = AS_CHARACTER(seqsP));
    numProtect++;
    seqs = smalloc(nseqs*sizeof(char*));
    for (i=0; i<nseqs; i++) {
      seqs[i] = smalloc((strlen(CHAR(STRING_ELT(seqsP, i)))+1)*sizeof(char));
      strcpy(seqs[i], CHAR(STRING_ELT(seqsP, i)));
    }
  }
  if (alphabetP != R_NilValue) {
    PROTECT(alphabetP = AS_CHARACTER(alphabetP));
    numProtect++;
    alphabet = smalloc((strlen(CHARACTER_VALUE(alphabetP))+1)*sizeof(char));
    strcpy(alphabet, CHARACTER_VALUE(alphabetP));
  }

  msa = msa_new(seqs, names, nseqs, length, alphabet);
  if (msa->length > 0 && ! ordered) {
    ss_from_msas(msa, 1, 0, NULL, NULL, NULL, -1, 0);
  } else if (idxOffsetP != R_NilValue)
    msa->idx_offset = INTEGER_VALUE(idxOffsetP);

  //don't free seqs or names because they are used by MSA object
  if (alphabet != NULL) sfree(alphabet);
  
  PROTECT(result = rph_msa_new_extptr(msa));
  numProtect++;

  UNPROTECT(numProtect);
  return result;
}

SEXP rph_msa_copy(SEXP msa) {
  return rph_msa_new_extptr(msa_create_copy((MSA*)EXTPTR_PTR(msa), 0));
}
			   


/* Not really sure if this will work, but seems better than 
   nothing (or is it worse)? How can we test? */
SEXP rph_is_msa(SEXP msaP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  return ScalarLogical(msa != NULL && (msa->seqs != NULL || msa->ss != NULL) &&
		       msa->names != NULL && msa->length > 0);
}

SEXP rph_msa_base_freq(SEXP msaP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  ListOfLists *lol = lol_new(3);
  char **alph;
  int i, alph_size = strlen(msa->alphabet);
  Vector *base_counts;
  double sum=0.0;
  alph = smalloc(alph_size*sizeof(char*));
  for (i=0; i < alph_size; i++) {
    alph[i] = smalloc(2*sizeof(char));
    alph[i][0] = msa->alphabet[i];
    alph[i][1] = '\0';
  }
  lol_push_charvec(lol, alph, alph_size, "states");
  
  base_counts = msa_get_base_counts(msa, -1, -1);
  lol_push_dbl(lol, base_counts->data, base_counts->size, "counts");
  
  sum = 0;
  for (i=0; i < base_counts->size; i++)
    sum += vec_get(base_counts, i);
  if (sum != 0.0)
    vec_scale(base_counts, 1.0/sum);
  lol_push_dbl(lol, base_counts->data, base_counts->size, "freq");
  lol_set_class(lol, "data.frame");
  return rph_listOfLists_to_SEXP(lol);
}


/* Note: this changes the object passed in */
SEXP rph_msa_reduce_to_4d(SEXP msaP, SEXP gffP) {
  CategoryMap  *cm;
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  GFF_Feature *f;
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  String *fourD_refseq = NULL;
  int i, tuple_size;
  if (msa->ss != NULL && msa->ss->tuple_idx == NULL)
    die("cannot extract 4d sites with unordered representation of MSA");
  msa_register_protect(msa);
  gff_register_protect(gff);
  msa_free_categories(msa);
  if (msa->ss != NULL) ss_free_categories(msa->ss);
  for (i=0; i<lst_size(gff->features); i++) {
    checkInterruptN(i, 1000);
    f = lst_get_ptr(gff->features, i);
    if (f->frame == GFF_NULL_FRAME) f->frame = 0;
    if (fourD_refseq == NULL) 
      fourD_refseq = str_new_charstr(f->seqname->chars);
    else if (!str_equals(fourD_refseq, f->seqname))
      die("to obtain 4d sites, all features should have same source");
    if (str_equals_charstr(f->feature, "CDS") && f->strand != '-')
      str_cpy_charstr(f->feature, "CDSplus");
    else if (str_equals_charstr(f->feature, "CDS") && f->strand == '-')
      str_cpy_charstr(f->feature, "CDSminus");
  }
  if (fourD_refseq == NULL)
    die("ERROR rph_msa_reduce_to_4d: fourD_refseq is NULL\n");
  tuple_size = 1;
  if (msa->ss != NULL && msa->ss->tuple_size != 1) 
    ss_reduce_tuple_size(msa, tuple_size);

  cm = cm_new_string_or_file("NCATS=6; CDSplus 1-3; CDSminus 4-6");

  if (msa->idx_offset != 0) {
    for (i=0; i<lst_size(gff->features); i++) {
      f = lst_get_ptr(gff->features, i);
      f->start -= msa->idx_offset;
      f->end -= msa->idx_offset;
    }
    msa_map_gff_coords(msa, gff, -1, 0, 0);
  }
  msa_label_categories(msa, gff, cm);
  msa_strip_gaps(msa, msa_get_seq_idx(msa, fourD_refseq->chars)+1);

  reduce_to_4d(msa, cm);  //this returns SS with tuple size 3
  
  ss_reduce_tuple_size(msa, tuple_size);
  return msaP;
}


SEXP rph_msa_extract_feature(SEXP msaP, SEXP gffP) {
  int i, j, pos=0;
  GFF_Set *gff=NULL;
  MSA *msa;
  CategoryMap *cm;

  msa = (MSA*)EXTPTR_PTR(msaP);
  if (msa->ss != NULL && msa->ss->tuple_idx == NULL)
    die("ordered representation of alignment required to extract features");
  gff=(GFF_Set*)EXTPTR_PTR(gffP);
  msa_register_protect(msa);
  gff_register_protect(gff);
  cm = cm_new_from_features(gff);
  
  /* convert GFF to coordinate frame of alignment */
  if (msa->idx_offset !=0) {
    for (i=0; i<lst_size(gff->features); i++) {
      GFF_Feature *f = lst_get_ptr(gff->features, i);
      checkInterruptN(i, 1000);
      f->start -= msa->idx_offset;
      f->end -= msa->idx_offset;
    }
  }
  msa_map_gff_coords(msa, gff, -1, 0, 0);
  msa_label_categories(msa, gff, cm);

  if (msa->ss != NULL) {
    for (i=0; i<msa->length; i++) {
      checkInterruptN(i, 1000);
      if (msa->categories[i] == 0) {
	msa->ss->counts[msa->ss->tuple_idx[i]]--;
	if (msa->ss->counts[msa->ss->tuple_idx[i]] < 0)
	  die("ERROR msa->ss->counts[msa->ss->tuple_idx[%i]]=%i\n",
	      i, msa->ss->counts[msa->ss->tuple_idx[i]]);
      }
    }
    ss_remove_zero_counts(msa);
    sfree(msa->ss->tuple_idx);
    msa->ss->tuple_idx = NULL;
  }
  if (msa->seqs != NULL) {
    for (i=0; i<msa->length; i++) {
      checkInterruptN(i, 1000);
      if (msa->categories[i] > 0) {
	if (pos != i) {
	  for (j=0; j<msa->nseqs; j++) 
	    msa->seqs[j][pos] = msa->seqs[j][i];
	}
	pos++;
      }
    }
    for (j=0; j<msa->nseqs; j++) 
      msa->seqs[j][pos] = '\0';
  }
  msa_free_categories(msa);
  msa_update_length(msa);
  msa->idx_offset = 0;
  return msaP;
}





SEXP rph_msa_format_for_suffix(SEXP filenameP) {
  msa_format_type t = msa_format_for_suffix(CHARACTER_VALUE(filenameP));
  SEXP result;
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(msa_format_to_str(t)));
  UNPROTECT(1);
  return result;
}


SEXP rph_msa_format_for_content(SEXP filenameP) {
  FILE *infile = phast_fopen(CHARACTER_VALUE(filenameP), "r");
  msa_format_type t = msa_format_for_content(infile, 0);
  SEXP result;
  phast_fclose(infile);
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(msa_format_to_str(t)));
  UNPROTECT(1);
  return result;
}



/* TODO: fix so that if features has only one row it can return them ordered!
 */
SEXP rph_msa_read(SEXP filenameP, SEXP formatP, SEXP gffP, 
		  SEXP do4dP, SEXP alphabetP, 
		  SEXP tupleSizeP, SEXP refseqP, SEXP orderedP,
		  SEXP catsCycleP, SEXP docatsP, SEXP idxOffsetP, 
		  SEXP seqnamesP, SEXP discardSeqnamesP) {
  int do4d, tupleSize=1, ordered, i;
  GFF_Set *gff=NULL;
  msa_format_type fmt;
  CategoryMap *cm=NULL;
  MSA *msa;
  FILE *infile, *refseq=NULL;
  char *alphabet = NULL, *reverse_groups_tag=NULL;
  int numProtect=0, seq_keep=0, cycle_size=-1;
  List *cats_to_do=NULL, *cats_to_do_str, *seqnames;
  String *fourD_refseq = NULL;

  fmt = msa_str_to_format(CHARACTER_VALUE(formatP));
  if (fmt == UNKNOWN_FORMAT) 
    die("Unknown alignment format\n");
  if (alphabetP != R_NilValue) {
    alphabet = smalloc((strlen(CHARACTER_VALUE(alphabetP))+1)*sizeof(char));
    strcpy(alphabet, CHARACTER_VALUE(alphabetP));
  }
  if (refseqP != R_NilValue) 
    refseq = phast_fopen(CHARACTER_VALUE(refseqP), "r");
  if (tupleSizeP != R_NilValue)
    tupleSize = INTEGER_VALUE(tupleSizeP);
  if (gffP != R_NilValue) {
    gff = (GFF_Set*)EXTPTR_PTR(gffP);
    gff_register_protect(gff);
  }
  ordered = INTEGER_VALUE(orderedP);

  if (catsCycleP != R_NilValue)
    cycle_size = INTEGER_VALUE(catsCycleP);

  do4d = INTEGER_VALUE(do4dP);
  if (do4d) {
    cm = cm_new_string_or_file("NCATS=6; CDSplus 1-3; CDSminus 4-6");
    ordered=1;
    tupleSize=1;
    for (i=0; i<lst_size(gff->features); i++) {
      GFF_Feature *f = lst_get_ptr(gff->features, i);
      checkInterruptN(i, 1000);
      if (f->frame == GFF_NULL_FRAME) f->frame = 0;
      if (fourD_refseq == NULL) fourD_refseq = f->seqname;
      else if (!str_equals(fourD_refseq, f->seqname))
	die("get4d requires all features have same source column");
      if (str_equals_charstr(f->feature, "CDS") && f->strand != '-')
	str_cpy_charstr(f->feature, "CDSplus");
      else if (str_equals_charstr(f->feature, "CDS") && f->strand == '-')
	str_cpy_charstr(f->feature, "CDSminus");
    }
    cats_to_do = lst_new_int(6);
    for (i=1; i<=6; i++) lst_push_int(cats_to_do, i);
  } else if (catsCycleP != R_NilValue && docatsP != R_NilValue) {
    int *intp;
    PROTECT(docatsP = AS_INTEGER(docatsP));
    numProtect++;
    intp = INTEGER_POINTER(docatsP);
    cats_to_do = lst_new_int(LENGTH(docatsP));
    for (i=0; i < LENGTH(docatsP); i++)
      lst_push_int(cats_to_do, intp[i]);
  } else if (docatsP != R_NilValue) {
    int numcats, cmStrLen=100;
    char *cmStr;
    PROTECT(docatsP = AS_CHARACTER(docatsP));
    numProtect++;
    numcats = LENGTH(docatsP);
    for (i=0; i<LENGTH(docatsP); i++)
      cmStrLen += strlen(CHAR(STRING_ELT(docatsP, i))) + 10;
    cmStr = smalloc(cmStrLen*sizeof(char));
    sprintf(cmStr, "NCATS = %i", LENGTH(docatsP));
    for (i=0; i<LENGTH(docatsP); i++) {
      cmStrLen = (int)strlen(cmStr);
      sprintf(&cmStr[cmStrLen], "; %s %i", CHAR(STRING_ELT(docatsP, i)), i+1);
    }
    cm = cm_new_string_or_file(cmStr);
    cats_to_do_str = lst_new_ptr(numcats);
    for (i=0; i<numcats; i++) {
      lst_push_ptr(cats_to_do_str, 
		   str_new_charstr(CHAR(STRING_ELT(docatsP, i))));
      if (cm_get_category(cm, lst_get_ptr(cats_to_do_str, i))==0)
	Rf_warning("category %s not found in GFF", 
		   (CHAR(STRING_ELT(docatsP, i))));
    }
    cats_to_do = cm_get_category_list(cm, cats_to_do_str, TRUE);
  } else if (gff != NULL) {
    if (cm ==  NULL)
      cm = cm_new_from_features(gff);
    cats_to_do = lst_new_int(cm->ncats);
    for (i=1; i<=cm->ncats; i++)
      lst_push_int(cats_to_do, i);
  }
  if (seqnamesP == R_NilValue && discardSeqnamesP == R_NilValue) {
    seqnames = NULL;
  } else if (seqnamesP != R_NilValue) {
    if (discardSeqnamesP != R_NilValue)
      die("either seqnames or discard.seqnames must be NULL");
    seq_keep = 1;
    seqnames = lst_new_ptr(LENGTH(seqnamesP));
    for (i=0; i < LENGTH(seqnamesP); i++) {
      lst_push_ptr(seqnames, str_new_charstr(CHAR(STRING_ELT(seqnamesP, i))));
    }
  } else {
    if (discardSeqnamesP == R_NilValue)
      die("ERROR rph_msa_read discardSeqnames = NULL\n");
    seq_keep = 0;
    seqnames = lst_new_ptr(LENGTH(discardSeqnamesP));
    for (i = 0 ; i < LENGTH(discardSeqnamesP); i++) {
      lst_push_ptr(seqnames, str_new_charstr(CHAR(STRING_ELT(discardSeqnamesP, i))));
    }
  }

  infile = phast_fopen(CHARACTER_VALUE(filenameP), "r");
  if (fmt == MAF) {  //reads and automatically converts to SS format
    msa = maf_read_cats_subset(infile, refseq, tupleSize, alphabet, gff, cm, 
			       cycle_size, ordered, reverse_groups_tag,
			       NO_STRIP, 0, cats_to_do, seqnames, seq_keep);
  } else {
    msa = msa_new_from_file_define_format(infile, fmt, alphabet);
    if (fmt != SS && idxOffsetP != R_NilValue) 
      msa->idx_offset = INTEGER_VALUE(idxOffsetP);
  }
  if ((fmt != MAF || do4d) && gff != NULL) {
    if (msa->ss != NULL && msa->ss->tuple_idx == NULL)
      die("ordered representation of alignment required to extract features");
    if (msa->ss != NULL) ss_free_categories(msa->ss);
    
    /* convert GFF to coordinate frame of alignment */
    if (msa->idx_offset !=0) {
      for (i=0; i<lst_size(gff->features); i++) {
	GFF_Feature *f = lst_get_ptr(gff->features, i);
	f->start -= msa->idx_offset;
	f->end -= msa->idx_offset;
      }
    }
    msa_map_gff_coords(msa, gff, -1, 0, 0);
    if (reverse_groups_tag != NULL) { /*reverse complement by group */
      if (fmt == SS) {
	ss_to_msa(msa);
	ss_free(msa->ss);
	msa->ss = NULL;
      }
      gff_group(gff, reverse_groups_tag);
      msa_reverse_compl_feats(msa, gff, NULL);
    }
    msa_label_categories(msa, gff, cm);
  } else if (fmt != MAF && cycle_size > 0) {
    msa->categories = (int*)smalloc(msa->length*sizeof(int));
    msa->ncats = cycle_size;
    for (i=0; i < msa->length; i++)
      msa->categories[i] = (i%cycle_size)+1;
  }
  if (msa == NULL) die("ERROR reading %s\n", CHARACTER_VALUE(filenameP));
  if (tupleSize != 1 || cats_to_do != NULL) {
    if (msa->ss == NULL) {
      ss_from_msas(msa, tupleSize, ordered, cats_to_do, NULL, NULL, -1, 0);
    }
    else {
      if (msa->ss->tuple_size < tupleSize)
	die("ERROR: input tuple size must be at least as large as tupleSize");
      if (msa->ss->tuple_idx != NULL && ordered==0) {
	free(msa->ss->tuple_idx);
	msa->ss->tuple_idx = NULL;
      }
      if (msa->ss->tuple_size > tupleSize)
	ss_reduce_tuple_size(msa, tupleSize);
    }
  }
  if (do4d) {
    //might want to move this to gap strip option at some point
    msa_strip_gaps(msa, msa_get_seq_idx(msa, fourD_refseq->chars)+1);
    reduce_to_4d(msa, cm);
  }

  if (ordered==0 && msa->ss != NULL)
    if (msa->ss->tuple_idx != NULL)
      die("ERROR: rph_msa_read: msa->ss->tuple_idx is not NULL\n");
  if (ordered==0) msa->idx_offset = 0;

  if (refseq != NULL)
    phast_fclose(refseq);
  phast_fclose(infile);
  
  /* If we have caluclated sufficient statistics, they are more up-to-date
     than the sequences.  The sequences may contain the full alignment
     whereas the sufficient statistics only the part of the alignment
     that was requested.  Remove the sequences to prevent downstream
     phast functions from working with them */
  if (msa->ss != NULL) 
    msa_free_seqs(msa);
  
  /* This should probably be somewhere in the main code, but I'm not 
     sure where.
     It seems like, if unordered sufficient statistics are calculated, 
     then the length of the MSA should be reduced to the number of 
     tuples which have been kept after filtering.  But they never are.
   */
  msa_update_length(msa);

  // also, do we need category-specific counts anymore?  We just wanted
  // the subset of the alignment.  
  msa_free_categories(msa);

  if (numProtect > 0)
    UNPROTECT(numProtect);
  return rph_msa_new_extptr((void*)msa);
}


SEXP rph_msa_valid_fmt_str(SEXP formatP) {
  msa_format_type fmt;
  SEXP result;
  int *resultP;
  PROTECT( result = allocVector(LGLSXP, 1));
  resultP = LOGICAL_POINTER(result);
  fmt = msa_str_to_format(CHARACTER_VALUE(formatP));
  resultP[0] = (fmt != UNKNOWN_FORMAT);
  UNPROTECT(1);
  return result;
}


//print sequence in given format.  If format not valid,
//use FASTA, but give no warning.
SEXP rph_msa_printSeq(SEXP msaP, SEXP fileP, SEXP formatP, 
		      SEXP prettyPrintP) {
  MSA *msa;
  msa_format_type fmt;
  msa = (MSA*)EXTPTR_PTR(msaP);
  msa_register_protect(msa);
  fmt = msa_str_to_format(CHARACTER_VALUE(formatP));
  if (fmt == UNKNOWN_FORMAT) 
    fmt = FASTA;
  if (fileP != R_NilValue)
    msa_print_to_file(CHARACTER_VALUE(fileP), 
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
  msa_register_protect(msa);

  if (msa->ss != NULL) {
    /* rather than calling ss_to_msa, which will allocate an entire new
       sequence that we don't need, only convert one species at a time.
     */
    for  (seq = 0; seq < msa->nseqs; seq++) {
      tempseq = ss_get_one_seq(msa, seq);
      SET_STRING_ELT(result, seq, mkChar(tempseq));
      sfree(tempseq);
    }
  } else {
    for (seq = 0; seq < msa->nseqs; seq++) { 
      if (msa->seqs[seq][msa->length] != '\0')
	die("ERROR rph_msa_seqs: bad sequence terminator\n");
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


SEXP rph_msa_ninformative_sites(SEXP msaP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP rv;
  PROTECT(rv = allocVector(INTSXP, 1));
  INTEGER(rv)[0] = (int)msa_ninformative_sites(msa, -1);
  UNPROTECT(1);
  return rv;
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

SEXP rph_msa_isOrdered(SEXP msaP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP result;
  int *resultP;
  PROTECT(result = NEW_LOGICAL(1));
  resultP = LOGICAL_POINTER(result);
  if (msa->ss == NULL) resultP[0] = 1;
  else resultP[0] = (msa->ss->tuple_idx != NULL);
  UNPROTECT(1);
  return result;
}

SEXP rph_msa_idxOffset(SEXP msaP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP result;
  int *resultP;
  
  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  resultP[0] = msa->idx_offset;
  UNPROTECT(1);
  return result;
}


SEXP rph_msa_square_brackets(SEXP msaP, SEXP rowsP, SEXP colsP) {
  MSA *msa, *newMsa;
  char **names, **seqs;
  int *rows=NULL, *cols=NULL, i, j, spec, nrow, ncol, numprotect=0;

  msa = (MSA*)EXTPTR_PTR(msaP);
  msa_register_protect(msa);
  if (msa->ss != NULL && msa->ss->tuple_idx == NULL)
    ss_make_ordered(msa);

  if (rowsP != R_NilValue) {
    nrow = LENGTH(rowsP);
    PROTECT(rowsP = AS_INTEGER(rowsP));
    rows = INTEGER_POINTER(rowsP);
    numprotect++;
  } else nrow = msa->nseqs;


  if (colsP != R_NilValue) {
    ncol = LENGTH(colsP);
    PROTECT(colsP = AS_INTEGER(colsP));
    cols = INTEGER_POINTER(colsP);
    numprotect++;
  } else ncol = msa->length;

  names = smalloc(nrow*sizeof(char*));
  seqs = smalloc(nrow*sizeof(char*));
  for (i=0; i < nrow; i++) {
    checkInterrupt();
    if (rows == NULL) spec = i;
    else {
      spec = rows[i]-1; //convert to 0-based numbers from R indices
      if (spec < 0 || spec >= msa->nseqs) 
	die("invalid row in rph_msa_square_brackets");
    }
    names[i] = copy_charstr(msa->names[spec]);

    if (cols == NULL) {
      if (msa->ss != NULL) seqs[i] = ss_get_one_seq(msa, spec);
      else seqs[i] = copy_charstr(msa->seqs[spec]);
    } else {
      seqs[i] = smalloc((ncol+1)*sizeof(char));
      for (j=0; j<ncol; j++) {
	checkInterrupt();
	if (cols[j] > msa->length) 
	  die("invalid column in rph_msa_square_brackets");
	seqs[i][j] = msa_get_char(msa, spec, cols[j]-1);
      }
      seqs[i][j] = '\0';
    }
  }
  newMsa = msa_new(seqs, names, nrow, ncol, msa->alphabet);
  if (numprotect > 0) UNPROTECT(numprotect);
  return rph_msa_new_extptr(newMsa);
}

// function for [<-.msa
// R code already checks that valueP can be recycled to expected length
// and that rowsP==colsP or rowsP==NULL or colsP==NULL.  If rowsP==NULL
// or colsP == NULL this indicates all rows/columns respectively
SEXP rph_msa_square_bracket_equals(SEXP msaP, SEXP rowsP, SEXP colsP,
				   SEXP valueP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  int i, j, *rows=NULL, *cols=NULL, nrow, ncol, numprotect=0, pos, row, col;
  msa_register_protect(msa);

  if (msa->ss != NULL) {
    //cannot have ss because changing the sequence will change the ss
    ss_to_msa(msa);
    ss_free(msa->ss);
    msa->ss = NULL;
  }
  if (rowsP != R_NilValue) {
    nrow = LENGTH(rowsP);
    PROTECT(rowsP = AS_INTEGER(rowsP));
    rows = INTEGER_POINTER(rowsP);
    numprotect++;
  } else nrow = msa->nseqs;
  
  if (colsP != R_NilValue) {
    ncol = LENGTH(colsP);
    PROTECT(colsP = AS_INTEGER(colsP));
    cols = INTEGER_POINTER(colsP);
    numprotect++;
  } else ncol = msa->length;

  pos = 0;
  for (i=0; i < ncol; i++) {
    if (colsP == R_NilValue) col=i;
    else col = cols[i]-1;
    for (j=0; j < nrow; j++) {
      checkInterrupt();
      if (rowsP == R_NilValue) row = j;
      else row = rows[j]-1;
      msa->seqs[row][col] = CHARACTER_VALUE(STRING_ELT(valueP, pos++))[0];
      if (pos == LENGTH(valueP)) pos = 0;
    }
  }
  if (numprotect > 0)
    UNPROTECT(numprotect);
  return msaP;
}

SEXP rph_msa_sub_alignment(SEXP msaP, SEXP seqsP, SEXP keepP, 
			   SEXP startcolP, SEXP endcolP, 
			   SEXP refseqNameP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP), *subMsa;
  int numseq, i, keep=1, startcol, endcol, numProtect=0, refseq, idx_offset=0;
  List *seqlist_str, *l=NULL;
  msa_coord_map *map = NULL;

  msa_register_protect(msa);
  if (msa->ss != NULL && msa->ss->tuple_idx == NULL)
    ss_make_ordered(msa);

  if (seqsP != R_NilValue) {
    PROTECT(seqsP = AS_CHARACTER(seqsP));
    numProtect++;
    numseq = LENGTH(seqsP);
    seqlist_str = lst_new_ptr(numseq);
    for (i=0; i<numseq; i++)  
      lst_push_ptr(seqlist_str, str_new_charstr(CHAR(STRING_ELT(seqsP, i))));
    l = msa_seq_indices(msa, seqlist_str);
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
    if (refseq == 1 && msa->idx_offset != 0) {   //assume idx_offset refers to first seq?
      if (startcolP != R_NilValue)
	startcol -= msa->idx_offset;
      if (endcolP != R_NilValue)
	endcol -= msa->idx_offset;
    }
    if (refseq == 1) idx_offset = startcol - 1 + msa->idx_offset;
    map = msa_build_coord_map(msa, refseq);
    startcol = msa_map_seq_to_msa(map, startcol);
    if (endcolP != R_NilValue)
      endcol = msa_map_seq_to_msa(map, endcol);
  }

  if (startcol < 1 || startcol > msa->length)
    die("start column (%i) out of range", startcol);
  if (endcol < 1 || endcol > msa->length)
    die("end column (%i) out of range", endcol);

  subMsa = msa_sub_alignment(msa, l, keep, startcol-1, endcol);
  subMsa->idx_offset = idx_offset;
  if (subMsa == NULL)
    die("ERROR rph_msa_sub_alignment got NULL subMsa\n");
  if (numProtect > 0) UNPROTECT(numProtect);
  return rph_msa_new_extptr(subMsa);
}


/*  Note msa_strip_gaps modifies original alignment.  This returns
    same pointer as was passed in.
 */
SEXP rph_msa_strip_gaps(SEXP msaP, SEXP stripModeP, SEXP allOrAnyGaps) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  int stripMode=-1, unordered;

  unordered=(msa->ss!=NULL && msa->ss->tuple_idx==NULL);
  msa_register_protect(msa);
  
  if (allOrAnyGaps != R_NilValue) {
    if (!strcmp(CHARACTER_VALUE(allOrAnyGaps), "all.gaps"))
      stripMode = STRIP_ALL_GAPS;
    else if (!strcmp(CHARACTER_VALUE(allOrAnyGaps), "any.gaps"))
      stripMode = STRIP_ANY_GAPS;
    else die("invalid strip.mode %s", CHARACTER_VALUE(allOrAnyGaps));
  } else stripMode = INTEGER_VALUE(stripModeP);

  msa_strip_gaps(msa, stripMode);

  /* sometimes stripping gaps restores order to the tuples where it
     shouldn't */
  if (unordered) {
    if (msa->ss == NULL) 
      ss_from_msas(msa, 1, 0, NULL, NULL, NULL, -1, 0);
    else if (msa->ss->tuple_idx != NULL) {
      sfree(msa->ss->tuple_idx);
      msa->ss->tuple_idx = NULL;
    }
  }
  return msaP;
}


SEXP rph_msa_postprob(SEXP msaP, SEXP tmP, SEXP doEverySiteP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  ListOfLists *result = lol_new(1);
  int doEverySite = LOGICAL_VALUE(doEverySiteP);
  if (msa->ss == NULL) {
    msa_register_protect(msa);
    ss_from_msas(msa, tm->order+1, doEverySite, NULL, NULL, NULL, -1, 0);
  }
  tm_register_protect(tm);
  print_post_prob_stats(tm, msa, NULL, 1, 0, 0, 0, doEverySite,
			-1, 1, result);
  return rph_listOfLists_to_SEXP(result);
}


SEXP rph_msa_exp_subs(SEXP msaP, SEXP tmP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  ListOfLists *result = lol_new(1);
  if (msa->ss == NULL) {
    msa_register_protect(msa);
    ss_from_msas(msa, tm->order+1, 0, NULL, NULL, NULL, -1, 0);
  }
  tm_register_protect(tm);
  print_post_prob_stats(tm, msa, NULL, 0, 1, 0, 0, 0, -1, 1, result);
  return rph_listOfLists_to_SEXP(result);
}


SEXP rph_msa_exp_tot_subs(SEXP msaP, SEXP tmP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  ListOfLists *result = lol_new(1);
  if (msa->ss == NULL) {
    ss_from_msas(msa, tm->order+1, 0, NULL, NULL, NULL, -1, 0);
    msa_protect(msa);
  }
  tm_register_protect(tm);
  print_post_prob_stats(tm, msa, NULL, 0, 0, 1, 0, 0, -1, 1, result);
  return rph_listOfLists_to_SEXP(result);
}

SEXP rph_msa_exp_col_subs(SEXP msaP, SEXP tmP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  ListOfLists *result = lol_new(1);
  if (msa->ss == NULL) {
    ss_from_msas(msa, tm->order+1, 0, NULL, NULL, NULL, -1, 0);
    msa_protect(msa);
  }
  tm_register_protect(tm);
  print_post_prob_stats(tm, msa, NULL, 0, 0, 0, 1, 0, -1, 1, result);
  return rph_listOfLists_to_SEXP(result);
}


SEXP rph_msa_likelihood(SEXP msaP, SEXP tmP, SEXP gffP, SEXP byColumnP) {
  int by_column, force_order=0, i, j, start, end;
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  TreeModel *tm = (TreeModel*)EXTPTR_PTR(tmP);
  double *col_probs=NULL, likelihood, *resultP, log2=log(2);
  SEXP result;
  GFF_Set *gff=NULL;
  GFF_Feature *feat;

  by_column = LOGICAL_VALUE(byColumnP);
  if (gffP != R_NilValue) {
    gff = (GFF_Set*)EXTPTR_PTR(gffP);
    if (by_column) die("cannot use by.column with features");
    gff_register_protect(gff);
  }
  if (msa->ss == NULL)
    msa_register_protect(msa);
  tm_register_protect(tm);

  tm_set_subst_matrices(tm);

  if (by_column) {
    if (msa->ss != NULL && msa->ss->tuple_idx == NULL) {
      force_order = 1;
      msa->ss->tuple_idx = smalloc(msa->length*sizeof(int));
      for (i=0; i<msa->length; i++)
	msa->ss->tuple_idx[i] = i;
    }
    col_probs = smalloc(msa->length*sizeof(double));
    PROTECT(result = NEW_NUMERIC(msa->length));
  } else if (gff != NULL) {
    if (msa->ss != NULL && msa->ss->tuple_idx == NULL)
      die("cannot get likelihood for features of an un-ordered alignment");
    PROTECT(result = NEW_NUMERIC(lst_size(gff->features)));
    col_probs = smalloc(msa->length*sizeof(double));
  }
  else PROTECT(result = NEW_NUMERIC(1));
  resultP = NUMERIC_POINTER(result);


  likelihood = log2*tl_compute_log_likelihood(tm, msa, col_probs, NULL, -1, NULL);
  if (by_column) {
    for (i=0; i<msa->length; i++)
      resultP[i] = log2*col_probs[i];
  } else if (gff != NULL) {    
    if (msa->idx_offset != 0) {
      for (i=0; i < lst_size(gff->features); i++) {
	feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
	feat->start -= msa->idx_offset;
	feat->end -= msa->idx_offset;
      }
    }
    msa_map_gff_coords(msa, gff, -1, 0, 0);
    for (i=0; i < lst_size(gff->features); i++) {
      checkInterruptN(i, 1000);
      feat=(GFF_Feature*)lst_get_ptr(gff->features, i);
      start = max(feat->start-1, 0);
      end = min(feat->end, msa->length);
      resultP[i] = 0.0;
      for (j=start; j < end; j++)
	resultP[i] += col_probs[j];
      resultP[i] *= log2;
    }
  }
  else resultP[0] = likelihood;
  
  if (force_order) {
    sfree(msa->ss->tuple_idx);
    msa->ss->tuple_idx = NULL;
  }
  UNPROTECT(1);
  return result;
}


SEXP rph_msa_base_evolve(SEXP modP, SEXP nsitesP, SEXP hmmP, 
			 SEXP getFeaturesP) {
  TreeModel **mods;
  int nsites, nstate=1, i, *labels=NULL;
  MSA *msa;
  HMM *hmm=NULL;
  SEXP result, names_sexp;
  char **names;
  GFF_Set *feats;
  GFF_Feature *newfeat;
  int currstart, currstate;
  char *seqname, *src="base.evolve", temp[1000];

  GetRNGstate(); //seed R's random number generator
  nsites = INTEGER_VALUE(nsitesP);
  if (hmmP != R_NilValue) {
    hmm = (HMM*)EXTPTR_PTR(hmmP);
    hmm_register_protect(hmm);
    nstate = hmm->nstates;
    if (LOGICAL_VALUE(getFeaturesP))
      labels = smalloc(nsites*sizeof(int));
  }
  mods = smalloc(nstate*sizeof(TreeModel*));
  for (i=0; i<nstate; i++) {
    mods[i] = (TreeModel*)EXTPTR_PTR(VECTOR_ELT(modP, i));
    tm_register_protect(mods[i]);
  }
  msa = tm_generate_msa(nsites, hmm, mods, labels);
  seqname = msa->names[0];
  PutRNGstate();
  if (labels != NULL) {
    ListOfLists *result_lol = lol_new(2);
    feats = gff_new_set();
    names = smalloc(nstate*sizeof(char*));
    PROTECT(names_sexp = GET_NAMES(modP));
    if (names_sexp == R_NilValue) {
      for (i=0; i < nstate; i++) {
	sprintf(temp, "state%i", i+1);
	names[i] = copy_charstr(temp);
      }
    } else {
      for (i = 0 ; i < nstate; i++)
	names[i] = copy_charstr(CHAR(STRING_ELT(names_sexp, i)));
    }
    currstart = 0;
    currstate = labels[0];
    for (i=1; i < nsites; i++) {
      checkInterruptN(i, 1000);
      if (labels[i] != currstate) {
	sprintf(temp, "id \"%s\"", names[currstate]);
	newfeat = gff_new_feature_copy_chars(seqname, src, 
					     names[currstate], 
					     currstart+1, i, 0, '+', 
					     GFF_NULL_FRAME, temp, TRUE);
	lst_push_ptr(feats->features, newfeat);
	currstate = labels[i];
	currstart = i;
      }
    }
    sprintf(temp, "id \"%s\"", names[currstate]);
    newfeat = gff_new_feature_copy_chars(seqname, src, names[currstate],
					 currstart+1, i, 0, '+', 
					 GFF_NULL_FRAME,
					 temp, TRUE);
    lst_push_ptr(feats->features, newfeat);
    lol_push_msa_ptr(result_lol, msa, "msa");
    lol_push_gff_ptr(result_lol, feats, "feats");
    PROTECT(result = rph_listOfLists_to_SEXP(result_lol));
    UNPROTECT(2);
    return result;
  }
  return rph_msa_new_extptr((void*)msa);
}


SEXP rph_msa_concat(SEXP aggregate_msaP, SEXP source_msaP) {
  MSA *aggregate_msa = (MSA*)EXTPTR_PTR(aggregate_msaP);
  MSA *source_msa = (MSA*)EXTPTR_PTR(source_msaP);
  msa_register_protect(aggregate_msa);
  msa_concatenate(aggregate_msa, source_msa);
  return aggregate_msaP;
}

SEXP rph_lst_new_extptr(List *l);

SEXP rph_msa_split_by_gff(SEXP msaP, SEXP gffP) {
  MSA *msa, **msas;
  ListOfLists *rv;
  GFF_Set *gff;
  int i, num_msa;
  
  msa = (MSA*)EXTPTR_PTR(msaP);
  gff = (GFF_Set*)EXTPTR_PTR(gffP);
  msa_register_protect(msa);
  gff_register_protect(gff);

  msas = msa_split_by_gff(msa, gff);
  if (msas == NULL) return R_NilValue;
  num_msa = lst_size(gff->features);
  rv = lol_new(num_msa);
  for (i=0; i < num_msa; i++) 
    lol_push_msa_ptr(rv, msas[i], NULL);
  return rph_listOfLists_to_SEXP(rv);
}


SEXP rph_msaList_get(SEXP listP, SEXP idxP) {
  List* l= (List*)EXTPTR_PTR(listP);
  int idx = INTEGER_VALUE(idxP);
  MSA *msa;

  msa = lst_get_ptr(l, idx-1);
  return rph_msa_new_extptr(msa);
}



SEXP rph_msa_reverse_complement(SEXP msaP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  msa_register_protect(msa);
  msa_reverse_compl(msa);
  return msaP;
}


SEXP rph_msa_informative_feats(SEXP msaP,
			       SEXP minInformativeP,
			       SEXP specP,
			       SEXP refseqP,
			       SEXP gapsAreInformativeP) {
  GFF_Set *feats;
  MSA *msa;
  int min_informative, *spec, numprotect=0, refseq, i, gaps_inf;
  List *speclist=NULL;
  SEXP rph_gff_new_extptr(GFF_Set *gff);

  msa = (MSA*)EXTPTR_PTR(msaP);
  msa_register_protect(msa);
  if (msa->ss != NULL && msa->ss->tuple_idx == NULL)
    ss_make_ordered(msa);
  min_informative = INTEGER_VALUE(minInformativeP);
  gaps_inf = LOGICAL_VALUE(gapsAreInformativeP);
  if (specP != R_NilValue) {
    PROTECT(specP = AS_INTEGER(specP));
    numprotect++;
    spec = INTEGER_POINTER(specP);
    speclist = lst_new_int(LENGTH(specP));
    for (i = 0 ; i < LENGTH(specP); i++) {
      lst_push_int(speclist, spec[i]-1);  //convert to zero-based species indices
    }
  }
  refseq = INTEGER_VALUE(refseqP);

  feats = msa_get_informative_feats(msa, min_informative, speclist,
				    refseq, gaps_inf);
  if (numprotect > 0) UNPROTECT(numprotect);
  return rph_gff_new_extptr(feats);
}


SEXP rph_msa_codon_clean(SEXP msaP, SEXP refseqP, SEXP strandP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  char strand='+';
  if (strandP != R_NilValue) {
    if (strcmp("+", CHARACTER_VALUE(strandP))==0);
    else if (strcmp("-", CHARACTER_VALUE(strandP))==0)
      strand='-';
    else die("Unknown strand %s\n", CHARACTER_VALUE(strandP));
  }
  msa_register_protect(msa);
  msa_codon_clean(msa, CHARACTER_VALUE(refseqP), strand);
  return msaP;
}


SEXP rph_msa_get_base_freqs_tuples(SEXP msaP, SEXP modP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  TreeModel *mod = (TreeModel*)EXTPTR_PTR(modP);
  int nstate, i;
  SEXP rv;
  double *doubleP;

  if (msa->ss == NULL) {
    msa_register_protect(msa);
    ss_from_msas(msa, mod->order+1, 0, NULL, NULL, NULL, -1,
                 subst_mod_is_codon_model(mod->subst_mod));
  }
  nstate = int_pow((int)strlen(msa->alphabet), mod->order+1);
  if (mod->backgd_freqs == NULL) {
    tm_register_protect(mod);
    mod->backgd_freqs = vec_new(nstate);
  }
  msa_get_base_freqs_tuples(msa, mod->backgd_freqs, mod->order+1, -1);
  PROTECT(rv = NEW_NUMERIC(nstate));
  doubleP = NUMERIC_POINTER(rv);
  for (i=0; i < nstate; i++) 
    doubleP[i] = vec_get(mod->backgd_freqs, i);
  UNPROTECT(1);
  return rv;
}


SEXP rph_msa_freq3x4(SEXP msaP) {
  Vector *backgd = vec_new(64);
  double *doubleP;
  int i;
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP rv;

  msa_get_backgd_3x4(backgd, msa);
  PROTECT(rv = NEW_NUMERIC(64));
  doubleP = NUMERIC_POINTER(rv);
  for (i=0; i < 64; i++)
    doubleP[i] = vec_get(backgd, i);
  UNPROTECT(1);
  return rv;
}


SEXP rph_msa_fraction_pairwise_diff(SEXP msaP, SEXP seq1P, SEXP seq2P, 
				    SEXP ignoreMissingP, SEXP ignoreGapsP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  SEXP rv;
  double *dp, val;
  ListOfLists *lol;
  Matrix *m;
  int ignore_missing = LOGICAL_VALUE(ignoreMissingP),
    ignore_gaps = LOGICAL_VALUE(ignoreGapsP), seq1=-1, seq2=-1;
  
  
  if (seq1P != R_NilValue) seq1 = INTEGER_VALUE(seq1P)-1;
  if (seq2P != R_NilValue) seq2 = INTEGER_VALUE(seq2P)-1;

  msa_register_protect(msa);
  if (seq1 != -1 && seq2 != -1) {
    PROTECT(rv = NEW_NUMERIC(1));
    dp = NUMERIC_POINTER(rv);
    dp[0] = msa_fraction_pairwise_diff(msa, seq1, seq2,
				       ignore_missing, ignore_gaps);
  } else if (seq1 != -1) {
    PROTECT(rv = NEW_NUMERIC(msa->nseqs));
    dp = NUMERIC_POINTER(rv);
    for (seq2=0; seq2 < msa->nseqs; seq2++)
      dp[seq2] = msa_fraction_pairwise_diff(msa, seq1, seq2,
					    ignore_missing, ignore_gaps);
  } else {
    m = mat_new(msa->nseqs, msa->nseqs);
    for (seq1=0; seq1 < msa->nseqs; seq1++) {
      mat_set(m, seq1, seq1, 0.0);
      for (seq2=seq1+1; seq2 < msa->nseqs; seq2++) {
	val = msa_fraction_pairwise_diff(msa, seq1, seq2,
					 ignore_missing, ignore_gaps);
	mat_set(m, seq1, seq2, val);
	mat_set(m, seq2, seq1, val);
      }
    }
    lol = lol_new(1);
    lol_push_matrix(lol, m, "pairwise.diff");
    PROTECT(rv = rph_listOfLists_to_SEXP(lol));
  }
  UNPROTECT(1);
  return rv;
}


SEXP rph_msa_translate(SEXP msaP, SEXP oneFrameP, SEXP frameP) {
  MSA *msa = (MSA*)EXTPTR_PTR(msaP);
  int i, *frame, oneframe;
  char **trans;
  SEXP result;
  oneframe = LOGICAL_VALUE(oneFrameP);
  frame = INTEGER_POINTER(frameP);

  trans = msa_translate(msa, oneframe, frame);

  PROTECT(result = NEW_CHARACTER(msa->nseqs));
  for (i=0; i < msa->nseqs; i++)
    SET_STRING_ELT(result, i, mkChar(trans[i]));
  UNPROTECT(1);
  return result;
  
}


