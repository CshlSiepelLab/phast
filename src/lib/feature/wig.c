/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: wig.c,v 1.37 2008-11-12 02:07:59 acs Exp $ */

#include <misc.h>
#include <wig.h>
#include <gff.h>
#include <stringsplus.h>


/* returns 1 if line is a wig header, 0 otherwise
   sets fixed=1 if a fixedStep, 0 for variableStep
   if fixed, sets chrom, start, step, and span (span=1 if not given, all else should be there or error).
   if variable, sets chrom and span (span=1 if not given.  start and step will not be set).
   Will return 0 and not set anything if not valid wig header.
   Can send fixed, chrom, start, step, span as NULL pointers if just want to
   test a line to see if it is a wig header.
 */
int wig_parse_header(String *line, int *fixed, char *chrom, 
		     int *start, int *step, int *span) {
  List *substrs=NULL, *fieldEquals=NULL;
  String *field,  *val;
  int haveChrom=0, haveStart=0, haveStep=0, isValidWig=1, isFixed, 
    spanVal=1, startVal, stepVal, numfield, i;
  char chromVal[STR_LONG_LEN];

  if (str_starts_with_charstr(line, "fixedStep"))
    isFixed = 1;
  else if (str_starts_with_charstr(line, "variableStep"))
    isFixed = 0;
  else return 0;

  substrs = lst_new_ptr(5);
  numfield = str_split(line, NULL, substrs);
  if (isFixed && numfield != 4 && numfield != 5) {
    isValidWig = 0; goto parseHeader_return;
  }
  if ((!isFixed) && numfield != 2 && numfield != 3) {
    isValidWig = 0; goto parseHeader_return;
  }

  fieldEquals = lst_new_ptr(2);
  for (i=1; i < lst_size(substrs); i++) {
    if (2 != str_split(lst_get_ptr(substrs, i), "=", fieldEquals)) {
      isValidWig = 0; goto parseHeader_return;
    }
    field = lst_get_ptr(fieldEquals, 0);
    val = lst_get_ptr(fieldEquals, 1);
    if (str_equals_charstr(field, "chrom")) {
      haveChrom = 1;
      strcpy(chromVal, val->chars);
    } else if (str_equals_charstr(field, "span")) {
      if (0 != str_as_int(val, &spanVal)) {
	isValidWig = 0; goto parseHeader_return;
      }
    } else if (str_equals_charstr(field, "start") && isFixed) {
      haveStart = 1;
      if (0 != str_as_int(val, &startVal)) {
	isValidWig=0; goto parseHeader_return;
      }
    } else if (str_equals_charstr(field, "step") && isFixed) {
      haveStep = 1;
      if (0 != str_as_int(val, &stepVal)) {
	isValidWig=0; goto parseHeader_return;
      }
    } else {
      isValidWig=0; goto parseHeader_return;
    }
    lst_free_strings(fieldEquals);
  }
  parseHeader_return:
  if (fieldEquals != NULL) {
    lst_free_strings(fieldEquals);
    lst_free(fieldEquals);
  }
  lst_free_strings(substrs);
  lst_free(substrs);

  if (isValidWig) {
    if (!haveChrom) isValidWig=0;
    else if (isFixed && (haveStart==0 || haveStep==0)) isValidWig=0;
    else {
      if (fixed != NULL) *fixed = isFixed;
      if (chrom != NULL) strcpy(chrom, chromVal);
      if (span != NULL) *span = spanVal;
      if (isFixed) {
	if (start != NULL) *start = startVal;
	if (step != NULL) *step = stepVal;
      }
    }
  }
  return isValidWig;
}
  

GFF_Set *gff_read_wig(FILE *F) {
  String *line = str_new(STR_LONG_LEN);
  char chrom[STR_LONG_LEN];
  int span, fixed, start, step, numfield;
  double score;
  GFF_Set *gff = gff_new_set();
  GFF_Feature *newfeat;
  List *substrs=NULL;
  
  chrom[0]='\0';
  while (EOF != str_readline(line, F)) {
    str_double_trim(line);
    if (!wig_parse_header(line, &fixed, chrom, &start, &step, &span)) {
      if (chrom[0] == '\0')
	die("Error reading wig file; no initial header line\n");
      if (fixed) {
	if (0 != str_as_dbl(line, &score))
	  die("Error parsing score from wig file, line=%s\n", line->chars);
      } else {
	if (substrs == NULL) substrs = lst_new_ptr(2);
	numfield = str_split(line, NULL, substrs);
	if (numfield != 2) 
	  die("Error parsing variableStep wig file, expected 2 fields, got line=%s\n", line->chars);
	if (0 != str_as_int(lst_get_ptr(substrs, 0), &start))
	  die("Error parsing variableStep wig file; first column should be integer in line %s\n", line->chars);
	if (0 != str_as_dbl(lst_get_ptr(substrs, 1), &score))
	  die("Error parsing variableStep wig file; second column should be score in line %s\n", line->chars);
	lst_free_strings(substrs);
      }
      newfeat = gff_new_feature_copy_chars(chrom, 
					   fixed ? "fixedWig" : "variableWig",
					   "wig_feature",
					   start, 
					   start + span - 1,
					   score, 
					   '+',
					   GFF_NULL_FRAME, ".", FALSE);
      lst_push_ptr(gff->features, newfeat);
      if (fixed) start += step;
    }
  }
  str_free(line);
  if (substrs != NULL) lst_free(substrs);
  return gff;
}


/* prints wig in fixedStep format.  
   Side-effect: groups features by seqname and sorts */
void wig_print(FILE *outfile, GFF_Set *set) {
  GFF_FeatureGroup *group;
  GFF_Feature *feat;
  int i, j, span=-1, step=-1, lastStart;
  gff_group_by_seqname(set);
  gff_sort(set);
  
  /* first make sure we can print this as wig.  Feature length needs
     to be the same across file */
  for (i=0; i < lst_size(set->groups); i++) {
    group = lst_get_ptr(set->groups, i);
    lastStart = -1;
    for (j=0; j < lst_size(group->features); j++) {
      feat = lst_get_ptr(group->features, j);
      if (feat->score_is_null) 
	die("Cannot print this feature as wig; some elements have no score");
      if (span == -1) span = feat->end - feat->start + 1;
      else if (feat->end - feat->start + 1 != span) 
	die("Cannot print this feature as a fixedStep wig due to variable span\n (all elements should have same length)");
      if (j != 0) {
	if (step == -1 || feat->start - lastStart < step)
	  step = feat->start - lastStart;
      }
      lastStart = feat->start;
    }
  }
  
  for (i=0; i < lst_size(set->groups); i++) {
    group = lst_get_ptr(set->groups, i);
    for (j=0; j < lst_size(group->features); j++) {
      feat = lst_get_ptr(group->features, j);
      if (j==0 || lastStart + step != feat->start) {
	fprintf(outfile, "fixedStep chrom=%s start=%i step=%i",
		feat->seqname->chars, feat->start, step);
	if (span != 1) fprintf(outfile, " span=%i\n", span);
	else fprintf(outfile, "\n");
	
      }
      fprintf(outfile, "%g\n", feat->score);
      lastStart = feat->start;
    }
  }
}

	
