/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* \file bed.c
    Reading and writing of BED files.  See http://genome.ucsc.edu/goldenPath/help/customTrack.html
    \ingroup feature
*/

#include <gff.h>
#include <ctype.h>
#include <hashtable.h>
#include <misc.h>

/** Fill out a GFF_Set from a BED file. */
void gff_read_from_bed(GFF_Set *gff, FILE *F) {
  String *line = str_new(STR_MED_LEN);
  List *l = lst_new_ptr(12), *block_sizes = lst_new_ptr(10), 
    *block_starts = lst_new_ptr(10);
  int i, is_error = 0, lineno = 0, id = 1;
  char group[STR_MED_LEN];
  Hashtable *hash = hsh_new(10000);

  /* source and feature type are constant and unknown */
  String *source = str_new_charstr("bed");
  String *feature = str_new_charstr("bed_feature");

  while (str_readline(line, F) != EOF) {
    lineno++;
    checkInterruptN(lineno, 1000);

    str_trim(line);
    if (line->length == 0) continue;

    str_split(line, "\t", l);
    if (str_equals_nocase_charstr(lst_get_ptr(l, 0), "Track")) {
      /* for now do nothing with Track info */
    }
    else {
      int start = 0, end = 0, score = 0, score_is_null = 1;
      String *chrom = NULL;
      char strand = '.';

      if (lst_size(l) < 3 ||
          str_as_int(lst_get_ptr(l, 1), &start) != 0 ||
          str_as_int(lst_get_ptr(l, 2), &end) != 0)  
        is_error = 1;
      else {
        chrom = lst_get_ptr(l, 0);
        start++;                /* switch to GFF coord convention */
      }

      if (lst_size(l) >= 4) {
        String *bed_name = lst_get_ptr(l, 3);
        int num;
        /* make sure unique name is assigned */
        if ((num = hsh_get_int(hash, bed_name->chars)) > 0) {
          num++;
          sprintf(group, "id \"%s.%d\"", bed_name->chars, num);
          hsh_reset_int(hash, bed_name->chars, num);
        }
        else {
          sprintf(group, "id \"%s\"", bed_name->chars);
          hsh_put_int(hash, bed_name->chars, 1);
        }
      }
      else
        sprintf(group, "id \"bed.%d\"", id++);
        
      if (lst_size(l) >= 5) {
        if ((str_as_int(lst_get_ptr(l, 4), &score) != 0)) {
	  phast_warning("score columns should contain integer\n");
          is_error = 1;
	}
        score_is_null = 0;
      }
      if (lst_size(l) >= 6 && 
          (((String*)lst_get_ptr(l, 5))->length != 1 ||
           ((strand = ((String*)lst_get_ptr(l, 5))->chars[0]) != '+' &&
            strand != '-' && strand != '.')))
        is_error = 1;
      if (lst_size(l) >= 10 && !is_error) { /* multiple features for line */
        if (lst_size(l) < 12) is_error = 1;
        else {
          int bl_size = 0, bl_start = 0;
          /* just ignore block count */
          str_split(lst_get_ptr(l, 10), ",", block_sizes);
          str_split(lst_get_ptr(l, 11), ",", block_starts);
          if (lst_size (block_sizes) != lst_size(block_starts)) is_error = 1;
          for (i = 0; !is_error && i < lst_size(block_sizes); i++) {
            if (str_as_int(lst_get_ptr(block_sizes, i), &bl_size) != 0 ||
                str_as_int(lst_get_ptr(block_starts, i), &bl_start) != 0)
              is_error = 1;
            else {
              lst_push_ptr(gff->features, 
                           gff_new_feature(str_dup(chrom), str_dup(source), 
                                           str_dup(feature), bl_start + start, 
                                           bl_start + start + bl_size - 1, score, 
                                           strand, GFF_NULL_FRAME, 
                                           str_new_charstr(group), 0));
            }
          }          
          lst_free_strings(block_sizes);
          lst_free_strings(block_starts);
        }
      }      
      else if (!is_error) {        /* single feature for line */
        lst_push_ptr(gff->features, 
                     gff_new_feature(str_dup(chrom), str_dup(source), 
                                     str_dup(feature), start, end, score, 
                                     strand, GFF_NULL_FRAME, 
                                     str_new_charstr(group), score_is_null));
      }

      if (is_error) 
        die("ERROR in line %d of BED file.\n", lineno);
    }

    lst_free_strings(l);
  }

  str_free(line);
  str_free(source);
  lst_free(l);
  lst_free(block_sizes);
  lst_free(block_starts);
  hsh_free(hash);
}

/** Write a GFF_Set in BED format. */
void gff_print_bed(FILE *OUTF,  GFF_Set *gff, int use_groups) {
  GFF_Feature *feat;
  int i, j;

  if (lst_size(gff->features) == 0) return; /* now can assume at least one feature */

  if (!use_groups) {
    Regex *tag_val_re = str_re_new("[[:alnum:]_.]+[[:space:]]+(\"[^\"]*\"|[^[:space:]]+)");
    List *l = lst_new_ptr(2);
    int ncols = 4;

    /* figure out how many columns to print.  Assume all features
       are consistent with first one */
    feat = lst_get_ptr(gff->features, 0);
    if (feat->strand != '.') ncols = 6;
    else if (!feat->score_is_null) ncols = 5;

    /* now print one line per feature */
    for (i = 0; i < lst_size(gff->features); i++) {
      String *name = NULL;
      checkInterruptN(i, 100);
      feat = lst_get_ptr(gff->features, i);

      /* try to extract name from attribute field (first value in
         tag-value pairs) */
      lst_clear(l);
      if (feat->attribute->length > 0 && 
          str_re_match(feat->attribute, tag_val_re, l, 1) >= 0) {
        name = lst_get_ptr(l, 1);
        str_remove_quotes(name);
      }
        
      /* required four columns */
      fprintf(OUTF, "%s\t%d\t%d\t%s", feat->seqname->chars, feat->start - 1, 
              feat->end, name == NULL ? "" : name->chars);

      /* optional additional columns */
      if (ncols >= 5)
        fprintf(OUTF, "\t%.0f", feat->score_is_null ? 0 : feat->score);
      
      if (ncols == 6)
        fprintf(OUTF, "\t%c", feat->strand);

      fprintf(OUTF, "\n");
      lst_free_strings(l);
    }

    str_re_free(tag_val_re);
    lst_free(l);
  }

  else {                        /* using groups */
    if (gff->groups == NULL)
      die("ERROR: groups required in gff_print_bed if use_groups == TRUE.\n");

    for (i = 0; i < lst_size(gff->groups); i++) {
      GFF_FeatureGroup *group = lst_get_ptr(gff->groups, i);
      double score = 0;
      checkInterrupt();

      if (lst_size(group->features) == 0) continue;

      /* if 'null' group, output individually */
      if (group->name == NULL || group->name->length == 0) {
        fprintf(stderr, "WARNING: skipping %d ungrouped features.\n",
                lst_size(group->features));
        continue;               /* a group may exist of features
                                   without tags; we'll ignore these */
      }

      /* use sum of non-null scores */
      for (j = 0; j < lst_size(group->features); j++) {
        feat = lst_get_ptr(group->features, j);
        if (!feat->score_is_null) score += feat->score;
      }

      feat = lst_get_ptr(group->features, 0);
      fprintf(OUTF, "%s\t%d\t%d\t%s\t%.0f\t%c\t%d\t%d\t0\t%d\t", 
              feat->seqname->chars, group->start - 1, group->end, 
              group->name->chars, score, feat->strand, 
              group->start - 1, group->end, lst_size(group->features));

      for (j = 0; j < lst_size(group->features); j++) {
        feat = lst_get_ptr(group->features, j);
        fprintf(OUTF, "%d,", feat->end - feat->start + 1);
      }    

      fprintf(OUTF, "\t");

      for (j = 0; j < lst_size(group->features); j++) 
        fprintf(OUTF, "%d,", 
                ((GFF_Feature*)lst_get_ptr(group->features, j))->start - 
                group->start);

      fprintf(OUTF, "\n");
    }
  }
}
