/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
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
#include <getopt.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <gff.h>
#include <misc.h>
#include <list_of_lists.h>
#include <rph_util.h>
#include <limits.h>

#include <Rdefines.h>


void rph_gff_free(SEXP gffPtr) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffPtr);
  phast_unregister_protected(gff);
  gff_free_set(gff);
}


SEXP rph_gff_new_extptr(GFF_Set *gff) {
  SEXP result;
  PROTECT(result=R_MakeExternalPtr((void*)gff, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, rph_gff_free, 1);
  gff_register_protect(gff);
  UNPROTECT(1);
  return result;
}


SEXP rph_gff_copy(SEXP gffP) {
  return rph_gff_new_extptr(gff_copy_set_no_groups((GFF_Set*)EXTPTR_PTR(gffP)));
}


SEXP rph_gff_read(SEXP filename) {
  FILE *infile = phast_fopen(CHARACTER_VALUE(filename), "r");
  SEXP rv;
  PROTECT(rv = rph_gff_new_extptr(gff_read_set(infile)));
  phast_fclose(infile);
  UNPROTECT(1);
  return rv;
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
  gff_register_protect(gff);

  len = lst_size(gff->features);

  //first five columns are required: name, src, feature, start, end
  PROTECT(names = allocVector(STRSXP, len));
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    SET_STRING_ELT(names, i, mkChar(feat->seqname->chars));
  }
  vec[0] = names;
  checkInterrupt();

  PROTECT(src = allocVector(STRSXP, len));
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    SET_STRING_ELT(src, i, mkChar(feat->source->chars));
  }
  vec[1] = src;
  checkInterrupt();

  PROTECT(feature=allocVector(STRSXP, len));
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    SET_STRING_ELT(feature, i, mkChar(feat->feature->chars));
  }
  vec[2] = feature;
  checkInterrupt();

  PROTECT(start=NEW_INTEGER(len));
  intp = INTEGER_POINTER(start);
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    intp[i] = feat->start;
  }
  vec[3] = start;
  checkInterrupt();

  PROTECT(end = NEW_INTEGER(len));
  intp = INTEGER_POINTER(end);
  for (i=0; i<len; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    intp[i] = feat->end;
  }
  vec[4] = end;
  checkInterrupt();


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
  checkInterrupt();


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
  checkInterrupt();
      
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
  checkInterrupt();

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
  checkInterrupt();
  
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
  else outfile = phast_fopen(CHARACTER_VALUE(filename), "w");
  
  gff_print_set(outfile, (GFF_Set*)EXTPTR_PTR(gff));
  if (outfile != stdout) phast_fclose(outfile);
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
  int *start, *end, frame=GFF_NULL_FRAME, *frameVec=NULL;
  double *scoreVec=NULL, score;
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
  }
  if (attributeP != R_NilValue) {
    PROTECT(attributeP = AS_CHARACTER(attributeP));
    haveAttribute=1;
  }

  numProtect += (haveScore + haveStrand + haveFrame + haveAttribute);

  gfflen = LENGTH(seqnameP);
  gff = gff_new_set_len(gfflen);

  for (i=0; i<gfflen; i++) {
    checkInterruptN(i, 1000);
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
  return rph_gff_new_extptr(gff);
}

SEXP rph_gff_minCoord(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  GFF_Feature *f;
  int i, mincoord=-1;
  SEXP rv;
  for (i=0; i<lst_size(gff->features); i++) {
    checkInterruptN(i, 1000);
    f = (GFF_Feature*)lst_get_ptr(gff->features, i);
    if (i==0 || f->start < mincoord)
      mincoord = f->start;
  }
  PROTECT(rv = allocVector(INTSXP, 1));
  INTEGER(rv)[0] = mincoord;
  UNPROTECT(1);
  return rv;
}



SEXP rph_gff_maxCoord(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  GFF_Feature *f;
  int i, maxcoord=-1;
  SEXP rv;
  for (i=0; i<lst_size(gff->features); i++) {
    checkInterruptN(i, 1000);
    f = (GFF_Feature*)lst_get_ptr(gff->features, i);
    if (i==0 || f->end > maxcoord)
      maxcoord = f->end;
  }
  PROTECT(rv = allocVector(INTSXP, 1));
  INTEGER(rv)[0] = maxcoord;
  UNPROTECT(1);
  return rv;
}


SEXP rph_gff_starts(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  GFF_Feature *f;
  int i;
  SEXP rv;
  PROTECT(rv = allocVector(INTSXP, lst_size(gff->features)));
  for (i=0; i<lst_size(gff->features); i++) {
    checkInterruptN(i, 1000);
    f = (GFF_Feature*)lst_get_ptr(gff->features, i);
    INTEGER(rv)[i] = f->start;
  }
  UNPROTECT(1);
  return rv;
}


SEXP rph_gff_ends(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  GFF_Feature *f;
  int i;
  SEXP rv;
  PROTECT(rv = allocVector(INTSXP, lst_size(gff->features)));
  for (i=0; i<lst_size(gff->features); i++) {
    checkInterruptN(i, 1000);
    f = (GFF_Feature*)lst_get_ptr(gff->features, i);
    INTEGER(rv)[i] = f->end;
  }
  UNPROTECT(1);
  return rv;
}


SEXP rph_gff_scores(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  GFF_Feature *f;
  int i;
  SEXP rv;
  PROTECT(rv = allocVector(REALSXP, lst_size(gff->features)));
  for (i=0; i<lst_size(gff->features); i++) {
    checkInterruptN(i, 1000);
    f = (GFF_Feature*)lst_get_ptr(gff->features, i);
    REAL(rv)[i] = f->score;
  }
  UNPROTECT(1);
  return rv;
}


SEXP rph_gff_strand(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  GFF_Feature *f;
  int i;
  SEXP rv;
  PROTECT(rv = allocVector(STRSXP, lst_size(gff->features)));
  for (i=0; i < lst_size(gff->features); i++) {
    checkInterruptN(i, 1000);
    f = (GFF_Feature*)lst_get_ptr(gff->features, i);
    if (f->strand == '+')
      SET_STRING_ELT(rv, i, mkChar("+"));
    else if (f->strand == '-')
      SET_STRING_ELT(rv, i, mkChar("-"));
    else SET_STRING_ELT(rv, i, mkChar("."));
  }
  UNPROTECT(1);
  return rv;
}



SEXP rph_gff_seqnames(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  GFF_Feature *f;
  int i;
  SEXP rv;
  PROTECT(rv = allocVector(STRSXP, lst_size(gff->features)));
  for (i=0; i < lst_size(gff->features); i++) {
    checkInterruptN(i, 1000);
    f = (GFF_Feature*)lst_get_ptr(gff->features, i);
    SET_STRING_ELT(rv, i, mkChar(f->seqname->chars));
  }
  UNPROTECT(1);
  return rv;
}


SEXP rph_gff_features(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  GFF_Feature *f;
  int i;
  SEXP rv;
  PROTECT(rv = allocVector(STRSXP, lst_size(gff->features)));
  for (i=0; i < lst_size(gff->features); i++) {
    checkInterruptN(i, 1000);
    f = (GFF_Feature*)lst_get_ptr(gff->features, i);
    SET_STRING_ELT(rv, i, mkChar(f->feature->chars));
  }
  UNPROTECT(1);
  return rv;
}


SEXP rph_gff_one_attribute(SEXP gffP, SEXP tagP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  GFF_Feature *f;
  ListOfLists *lol;
  List *l1, *l2;
  int numtag, numval, i, j, k, resultLen, maxResultLen=10;
  String *currStr, *tag, *currTag;
  char **result;
  SEXP rv;
  SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);

  
  if (lst_size(gff->features) == 0) return R_NilValue;
  gff_register_protect(gff);
  result = smalloc(maxResultLen*sizeof(char*));    
  tag = str_new_charstr(CHARACTER_VALUE(tagP));
  str_double_trim(tag);
  lol = lol_new(lst_size(gff->features));
  l1 = lst_new_ptr(10);
  l2 = lst_new_ptr(10);
  for (i=0; i < lst_size(gff->features); i++) {
    checkInterruptN(i, 1000);
    resultLen=0;
    f = (GFF_Feature*) lst_get_ptr(gff->features, i);
    numtag = str_split_with_quotes(f->attribute, ";", l1);  //split tags
    for (j=0; j < numtag; j++) {
      currStr = (String*)lst_get_ptr(l1, j);
      str_double_trim(currStr);
      
      //first try gff version 3, see if we have tag=val format
      numval = str_split_with_quotes(currStr, "=", l2);
      if (numval == 2) {
	currTag = (String*)lst_get_ptr(l2, 0);
	str_double_trim(currTag);
	if (str_equals(tag, currTag)) {  // tag matches target, add all values to list
	  currStr = str_new_charstr(((String*)lst_get_ptr(l2, 1))->chars);
	  lst_free_strings(l2);
	  numval = str_split_with_quotes(currStr, ",", l2);
	  str_free(currStr);
	  for (k=0; k < numval; k++) {
	    currStr = lst_get_ptr(l2, k);
	    str_double_trim(currStr);
	    str_remove_quotes(currStr);
	    if (resultLen > maxResultLen) {
	      maxResultLen += 100;
	      result = srealloc(result, maxResultLen*sizeof(char*));
	    }
	    result[resultLen++] = copy_charstr(currStr->chars);
	  }
	}
      } else {
	lst_free_strings(l2);
	
	//gff version 2
	//split into tag val val ... by whitespace unless enclosed in quotes
	numval =  str_split_with_quotes(currStr, NULL, l2);
	if (numval > 1) {
	  currStr = (String*)lst_get_ptr(l2, 0);
	  str_double_trim(currStr);
	  if (str_equals(tag, currStr)) {  //tag matches target, add all values to list
	    for (k=1; k < numval; k++) {
	      currStr = (String*)lst_get_ptr(l2, k);
	      str_double_trim(currStr);
	      str_remove_quotes(currStr);
	      if (resultLen > maxResultLen) {
		maxResultLen += 100;
		result = srealloc(result, maxResultLen*sizeof(char*));
	      }
	      result[resultLen++] = copy_charstr(currStr->chars);
	    }
	  }
	}
	lst_free_strings(l2);
      }
    }
    if (resultLen == 0) 
      result[resultLen++] = copy_charstr("");  //empty string will be converted to NA later
    lol_push_charvec(lol, result, resultLen, NULL);
    for (j=0; j < resultLen; j++) sfree(result[j]);
  }
  PROTECT(rv = rph_listOfLists_to_SEXP(lol));
  UNPROTECT(1);
  return rv;
}


SEXP rph_gff_overlapSelect(SEXP gffP, SEXP filter_gffP, 
			   SEXP numbaseOverlapP,
			   SEXP percentOverlapP,
			   SEXP nonOverlappingP,
			   SEXP overlappingFragmentsP) {
  GFF_Set *gff, *filter_gff, *overlapping_gff=NULL;
  int numbaseOverlap, nonOverlapping;
  double percentOverlap, overlappingFragments;

  gff = (GFF_Set*)EXTPTR_PTR(gffP);
  gff_register_protect(gff);
  filter_gff = (GFF_Set*)EXTPTR_PTR(filter_gffP);
  if (percentOverlapP == R_NilValue)
    percentOverlap = -1.0;
  else percentOverlap = NUMERIC_VALUE(percentOverlapP);
  if (nonOverlappingP == R_NilValue)
    nonOverlapping = FALSE;
  else nonOverlapping = LOGICAL_VALUE(nonOverlappingP);
  if (numbaseOverlapP == R_NilValue)
    numbaseOverlap = -1;
  else numbaseOverlap = INTEGER_VALUE(numbaseOverlapP);
  if (overlappingFragmentsP == R_NilValue)
    overlappingFragments = FALSE;
  else overlappingFragments = LOGICAL_VALUE(overlappingFragmentsP);
  
  if (overlappingFragments) overlapping_gff = gff_new_set();

  filter_gff = gff_overlap_gff(gff, filter_gff,
			       numbaseOverlap, percentOverlap, nonOverlapping,
			       overlappingFragments, overlapping_gff);
  if (overlappingFragments) {
    ListOfLists *rv = lol_new(2);
    lol_push_gff_ptr(rv, filter_gff, "frags");
    lol_push_gff_ptr(rv, overlapping_gff, "filter.frags");
    return rph_listOfLists_to_SEXP(rv);
  }
  return rph_gff_new_extptr(filter_gff);
}


SEXP rph_gff_add_UTRs(SEXP gffP) {
  GFF_Set *gff;
  gff = (GFF_Set*)EXTPTR_PTR(gffP);
  gff_register_protect(gff);
  gff_group(gff, "transcript_id");
  gff_create_utrs(gff);
  return gffP;
}


SEXP rph_gff_add_introns(SEXP gffP) {
  GFF_Set *gff;
  gff = (GFF_Set*)EXTPTR_PTR(gffP);
  gff_register_protect(gff);
  gff_group(gff, "transcript_id");
  gff_create_introns(gff);
  return gffP;
}


SEXP rph_gff_add_signals(SEXP gffP) {
  GFF_Set *gff;
  gff = (GFF_Set*)EXTPTR_PTR(gffP);
  gff_register_protect(gff);
  gff_group(gff, "transcript_id");
  gff_create_signals(gff);
  return gffP;
}


SEXP rph_gff_fix_start_stop(SEXP gffP) {
  GFF_Set *gff;
  gff = (GFF_Set*)EXTPTR_PTR(gffP);
  gff_register_protect(gff);
  gff_group(gff, "transcript_id");
  gff_fix_start_stop(gff);
  return gffP;
}




SEXP rph_gff_inverse(SEXP gffP, SEXP regionP) {
  GFF_Set *gff, *region, *notgff;
  gff = (GFF_Set*)EXTPTR_PTR(gffP);
  gff_register_protect(gff);
  region = (GFF_Set*)EXTPTR_PTR(regionP);  
  gff_register_protect(region);
  notgff = gff_inverse(gff, region);
  return rph_gff_new_extptr(notgff);
}


SEXP rph_gff_featureBits(SEXP gffListP, SEXP orP, SEXP returnGffP) {
  int numGff, i, j, or, returnGff;
  long numbit = 0;
  List *gfflist;
  GFF_Set *gff, *newgff=NULL;
  GFF_Feature *feat, *newfeat;
  SEXP rv;

  numGff = length(gffListP);
  gfflist = lst_new_ptr(numGff);
  //  Rf_PrintValue(gffListP);
  for (i = 0; i < numGff; i++) {
    gff = (GFF_Set*)EXTPTR_PTR(VECTOR_ELT(gffListP, i));
    lst_push_ptr(gfflist, gff);
    gff_register_protect(gff);
  }
  or = LOGICAL_VALUE(orP);
  returnGff = LOGICAL_VALUE(returnGffP);
  if (!or && numGff >= 2) {
    newgff = gff_overlap_gff(lst_get_ptr(gfflist, 0),
			     lst_get_ptr(gfflist, 1),
			     1, -1.0, FALSE, TRUE, NULL);
    numbit = gff_flatten_mergeAll(newgff);
    for (i=2; i < numGff; i++) {
      checkInterrupt();
      gff = gff_overlap_gff(newgff,
			    lst_get_ptr(gfflist, i),
			    1, -1.0, FALSE, TRUE, NULL);
      numbit = gff_flatten_mergeAll(gff);
      gff_free_set(newgff);
      newgff = gff;
    }
  } else {
    newgff = gff_new_set();
    for (i=0; i< numGff; i++) {
      gff = (GFF_Set*)lst_get_ptr(gfflist, i);
      for (j=0; j < lst_size(gff->features); j++) {
	checkInterruptN(j, 1000);
	feat = lst_get_ptr(gff->features, j);
	newfeat = gff_new_feature_copy(feat);
	lst_push_ptr(newgff->features, newfeat);
      }
    }
    numbit = gff_flatten_mergeAll(newgff);
  }
  if (returnGff) 
    return rph_gff_new_extptr(newgff);
  
  if (numbit > INT_MAX) {
    PROTECT(rv = allocVector(REALSXP, 1));
    REAL(rv)[0] = numbit;
  } else {
    PROTECT(rv = allocVector(INTSXP, 1));
    INTEGER(rv)[0] = numbit;
  }
  UNPROTECT(1);
  return rv;
}


SEXP rph_gff_append(SEXP gffListP) {
  GFF_Set *newgff = gff_new_set(), *gff;
  int i, j;
  for (i=0 ; i<length(gffListP); i++) {
    gff = (GFF_Set*)EXTPTR_PTR(VECTOR_ELT(gffListP, i));
    gff_register_protect(gff);
    for (j=0; j < lst_size(gff->features); j++) {
      checkInterruptN(j, 1000);
      lst_push_ptr(newgff->features, 
		   gff_new_feature_copy(lst_get_ptr(gff->features, j)));
    }
  }
  return rph_gff_new_extptr(newgff);
}


SEXP rph_gff_split(SEXP gffP, SEXP maxLengthP, SEXP dropP, SEXP splitFromRightP) {
  GFF_Set *gff, *newgff;
  int *maxlen, maxlen_size, *splitFromRight, splitFrom_size, drop;
  gff = (GFF_Set*)EXTPTR_PTR(gffP);
  gff_register_protect(gff);
  drop = LOGICAL_VALUE(dropP);
  PROTECT(maxLengthP = AS_INTEGER(maxLengthP));
  maxlen = INTEGER_POINTER(maxLengthP);
  maxlen_size = LENGTH(maxLengthP);
  PROTECT(splitFromRightP = AS_INTEGER(splitFromRightP));
  splitFromRight = INTEGER_POINTER(splitFromRightP);
  splitFrom_size = LENGTH(splitFromRightP);
  newgff = gff_split(gff, maxlen, maxlen_size, drop, 
		     splitFromRight, splitFrom_size);
  UNPROTECT(2);
  return rph_gff_new_extptr(newgff);
}


SEXP rph_gff_sort(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  gff_register_protect(gff);
  gff_group_by_seqname(gff);
  gff_sort(gff);
  return gffP;
}


SEXP rph_gff_nonOverlapping_genes(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  if (lst_size(gff->features) == 0) return gffP;
  gff_register_protect(gff);
  gff_group(gff, "transcript_id");
  gff_remove_overlaps(gff, NULL);
  return gffP;
}


SEXP rph_gff_flatten(SEXP gffP) {
  GFF_Set *gff = (GFF_Set*)EXTPTR_PTR(gffP);
  gff_register_protect(gff);
  gff_group_by_seqname(gff);
  gff_flatten_within_groups(gff);
  return gffP;
}
