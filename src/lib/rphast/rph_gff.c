/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_gff.c
The RPHAST handles to functions dealing with GFFs from
the phast package.

Melissa Hubisz
Last updated: 1/5/2010
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
#include <gff.h>

#include <Rdefines.h>

SEXP rph_gff_free(SEXP gffPtr) {
  gff_free_set((GFF_Set*)EXTPTR_PTR(gffPtr));
  return R_NilValue;
}


SEXP rph_gff_read(SEXP filename) {
  return R_MakeExternalPtr((void*)gff_read_set(fopen_fname(CHARACTER_VALUE(filename), "r")), R_NilValue, R_NilValue);
}


SEXP rph_gff_dataframe(SEXP gffPtr) {
  GFF_Set *gff;
  GFF_Feature *feat;
  SEXP result, names, src, feature, start, end, score, strand, frame, attribute, header;
  int i, len, listlen, *intp;
  double *doublep;
  char strandStr[2];
  //first five columns are required; others may not be defined
  char gffCols[9][20] = {"seqname", "src", "feature", "start", "end", "score", "strand", "frame", "attribute"};
  int have[9] = {1, 1, 1, 1, 1, 0, 0, 0, 0};
  int scorePos = 5, strandPos = 6, framePos = 7, attributePos = 8;
  SEXP vec[9];


  gff = (GFF_Set*)EXTPTR_PTR(gffPtr);

  len = lst_size(gff->features);

  //first five columns are required: name, src, feature, start, end
  PROTECT(names = allocVector(STRSXP, len));
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    SET_STRING_ELT(names, i, mkChar(feat->seqname->chars));
  }
  vec[0] = names;

  PROTECT(src = allocVector(STRSXP, len));
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    SET_STRING_ELT(src, i, mkChar(feat->source->chars));
  }
  vec[1] = src;
  
  PROTECT(feature=allocVector(STRSXP, len));
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    SET_STRING_ELT(feature, i, mkChar(feat->feature->chars));
  }
  vec[2] = feature;

  PROTECT(start=NEW_INTEGER(len));
  intp = INTEGER_POINTER(start);
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    intp[i] = feat->start;
  }
  vec[3] = start;

  PROTECT(end = NEW_INTEGER(len));
  intp = INTEGER_POINTER(end);
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    intp[i] = feat->end;
  }
  vec[4] = end;
    
  PROTECT(score = NEW_NUMERIC(len));
  doublep = NUMERIC_POINTER(score);
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    if (feat->score_is_null)
      doublep[i] = NA_REAL; //may have to include R_ext/Arith.h
    else {
      doublep[i] = feat->score;
      have[scorePos] = 1;
    }
  }
  vec[5] = score;

  PROTECT(strand = allocVector(STRSXP, len));
  strandStr[1] = '\0';
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    strandStr[0] = feat->strand;
    SET_STRING_ELT(strand, i, mkChar(strandStr));
    if (feat->strand != '.')
      have[strandPos] = 1;
  }
  vec[6] = strand;
      
  PROTECT(frame = NEW_INTEGER(len));
  intp = INTEGER_POINTER(frame);
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    if (feat->frame == GFF_NULL_FRAME)
      intp[i] = NA_INTEGER;
    else {
      have[framePos] = 1;
      intp[i] = feat->frame;
      if (feat->frame == 0)
	intp[i] = 0;
      else if (feat->frame==1)
	intp[i] = 2;
      else if (feat->frame==2)
	intp[i] = 1;
      else die("invalid frame %i in GFF", feat->frame);
    }
  }
  vec[7] = frame;

  PROTECT(attribute = allocVector(STRSXP, len));
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    //suspect mkChar is not dealing well with empty string?
    //    SET_STRING_ELT(attribute, i, mkChar(feat->attribute->chars));
    if (feat->attribute->length != 0) {
      have[attributePos] = 1;
      SET_STRING_ELT(attribute, i, mkChar(feat->attribute->chars));
    } else 
      SET_STRING_ELT(attribute, i, mkChar("."));
  }
  vec[8] = attribute;
  
  listlen = 0;
  for (i=0; i<9; i++) listlen += have[i];
  
  PROTECT(header = allocVector(STRSXP, listlen));
  PROTECT(result = allocVector(VECSXP, listlen));
  listlen = 0;
  for (i=0; i<9; i++) {
    if (have[i]) {
      SET_STRING_ELT(header, listlen, mkChar(gffCols[i]));
      SET_VECTOR_ELT(result, listlen++, vec[i]);
    }
  }
  SET_NAMES(result, header);

  UNPROTECT(11);
  return result;
}



SEXP rph_gff_print(SEXP filename, SEXP gff) {
  FILE *outfile;
  if (filename == R_NilValue)
    outfile = stdout;
  else outfile = fopen_fname(CHARACTER_VALUE(filename), "w");
  
  gff_print_set(outfile, (GFF_Set*)EXTPTR_PTR(gff));
  if (outfile != stdout) fclose(outfile);
  return R_NilValue;
}


SEXP rph_gff_numrow(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  SEXP result;
  int *resultP;
  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  resultP[0] = lst_size(gff->features);
  UNPROTECT(1);
  return result;
}


//take the elements of a GFF in R and make a GFF object in C; return pointer
//Assume length of vectors are all equal (except optional elements can be NULL)
SEXP rph_gff_new(SEXP seqnameP, SEXP srcP, SEXP featureP, SEXP startP, SEXP endP,
		 SEXP scoreP, SEXP strandP, SEXP frameP, SEXP attributeP) {
  GFF_Set *gff;
  GFF_Feature *feat;
  int gfflen, i;
  int haveScore=0, haveStrand=0, haveFrame=0, haveAttribute=0, numProtect=5;
  String *seqname, *source, *feature, *attribute;
  int *start, *end, frame, *frameVec;
  double *scoreVec, score;
  char strand;
  
  PROTECT(seqnameP = AS_CHARACTER(seqnameP));
  PROTECT(srcP = AS_CHARACTER(srcP));  
  PROTECT(featureP = AS_CHARACTER(featureP));  
  PROTECT(startP = AS_INTEGER(startP));
  start = INTEGER_POINTER(startP);
  PROTECT(endP = AS_INTEGER(endP));
  end = INTEGER_POINTER(endP);
  if (scoreP != R_NilValue) {
    PROTECT(scoreP = AS_NUMERIC(scoreP));
    haveScore = 1;
    scoreVec = NUMERIC_POINTER(scoreP);
  } else score=0;
  if (strandP != R_NilValue) {
    PROTECT(strandP = AS_CHARACTER(strandP));
    haveStrand=1;
  } else strand='.';
  if (frameP != R_NilValue) {
    PROTECT(frameP = AS_INTEGER(frameP));
    haveFrame=1;
    frameVec = INTEGER_POINTER(frameP);
  } else frame=GFF_NULL_FRAME;
  if (attributeP != R_NilValue) {
    PROTECT(attributeP = AS_CHARACTER(attributeP));
    haveAttribute=1;
  }

  numProtect += (haveScore + haveStrand + haveFrame + haveAttribute);

  gfflen = LENGTH(seqnameP);
  gff = gff_new_set_len(gfflen);

  for (i=0; i<gfflen; i++) {
    seqname = str_new_charstr(CHAR(STRING_ELT(seqnameP, i)));
    source = str_new_charstr(CHAR(STRING_ELT(srcP, i)));
    feature = str_new_charstr(CHAR(STRING_ELT(featureP, i)));
    if (haveScore) score = scoreVec[i];
    if (haveStrand) strand = (CHAR(STRING_ELT(strandP, i)))[0];
    if (haveFrame) {
      if (frameVec[i] == 0) frame = 0;
      else if (frameVec[i] == 1) frame = 2;
      else if (frameVec[i] == 2) frame = 1;
    }
    if (haveAttribute) attribute = str_new_charstr(CHAR(STRING_ELT(attributeP, i)));
    else attribute = str_new_charstr("");

    if (seqname == NULL) die("seqname is NULL\n");
    if (source == NULL) die ("source is NULL\n");
    if (feature ==  NULL) die("feature is NULL\n");
    if (attribute == NULL) die("attribute is NULL\n");
    if (strand != '+' && strand != '-' && strand!='.') die("strand is %c\n", strand);
    if (frame != GFF_NULL_FRAME && (frame<0 || frame>2)) die("frame is %i\n", frame);

    feat = gff_new_feature(seqname, source, feature, start[i], end[i], score, strand,
			   frame, attribute, haveScore==0);
    lst_push_ptr(gff->features, feat);
  }

  UNPROTECT(numProtect);
  return R_MakeExternalPtr((void*)gff, R_NilValue, R_NilValue);
}

