/* $Id: bed.c,v 1.4 2004-06-14 22:52:16 acs Exp $
   Written by Adam Siepel, 2004
   Copyright 2004, Adam Siepel, University of California */

/** \file bed.c
    Reading and writing of BED files.  See http://genome.ucsc.edu/goldenPath/help/customTrack.html
    \ingroup feature
*/

#include <gff.h>
#include <assert.h>
#include <ctype.h>
#include <misc.h>

/** Fill out a GFF_Set from a BED file. */
void gff_read_from_bed(GFF_Set *gff, FILE *F) {
  String *line = str_new(STR_MED_LEN);
  List *l = lst_new_ptr(12), *block_sizes = lst_new_ptr(10), 
    *block_starts = lst_new_ptr(10);
  int i, error = 0, lineno = 0;

  /* source and feature type are constant and unknown */
  String *source = str_new_charstr("bed");
  String *feature = str_new_charstr("bed_feature");

  while (str_readline(line, F) != EOF) {
    lineno++;

    str_trim(line);
    if (line->length == 0) continue;

    str_split(line, "\t", l);
    if (str_equals_nocase_charstr(lst_get_ptr(l, 0), "Track")) {
      /* for now do nothing with Track info */
    }
    else {
      int start = 0, end = 0, score = 0, score_is_null = 1;
      String *bed_name = NULL, *chrom = NULL;
      char strand = '.';

      if (lst_size(l) < 3 ||
          str_as_int(lst_get_ptr(l, 1), &start) != 0 ||
          str_as_int(lst_get_ptr(l, 2), &end) != 0)  
        error = 1;
      else {
        chrom = lst_get_ptr(l, 0);
        start++;                /* switch to GFF coord convention */
      }
      if (lst_size(l) >= 4) bed_name = lst_get_ptr(l, 3);
      if (lst_size(l) >= 5) {
        if ((str_as_int(lst_get_ptr(l, 4), &score) != 0 || 
           score < 0 || score > 1000)) 
          error = 1;
        score_is_null = 0;
      }
      if (lst_size(l) >= 6 && 
          (((String*)lst_get_ptr(l, 5))->length != 1 ||
           ((strand = ((String*)lst_get_ptr(l, 5))->chars[0]) != '+' &&
            strand != '-' && strand != '.')))
        error = 1;
      if (lst_size(l) >= 10 && !error) { /* multiple features for line */
        if (lst_size(l) < 12) error = 1;
        else {
          int bl_size = 0, bl_start = 0;
          /* just ignore block count */
          str_split(lst_get_ptr(l, 10), ",", block_sizes);
          str_split(lst_get_ptr(l, 11), ",", block_starts);
          if (lst_size (block_sizes) != lst_size(block_starts)) error = 1;
          for (i = 0; !error && i < lst_size(block_sizes); i++) {
            if (str_as_int(lst_get_ptr(block_sizes, i), &bl_size) != 0 ||
                str_as_int(lst_get_ptr(block_starts, i), &bl_start) != 0)
              error = 1;
            else {
              lst_push_ptr(gff->features, 
                           gff_new_feature(str_dup(chrom), str_dup(source), 
                                           str_dup(feature), bl_start + start, 
                                           bl_start + start + bl_size - 1, score, 
                                           strand, GFF_NULL_FRAME, 
                                           str_dup(bed_name), 0));
            }
          }          
          lst_free_strings(block_sizes);
          lst_free_strings(block_starts);
        }
      }      
      else if (!error) {        /* single feature for line */
        lst_push_ptr(gff->features, 
                     gff_new_feature(str_dup(chrom), str_dup(source), 
                                     str_dup(feature), start, end, score, 
                                     strand, GFF_NULL_FRAME, 
                                     bed_name != NULL ? str_dup(bed_name) : 
                                     str_new(STR_SHORT_LEN), score_is_null));
      }

      if (error) {
        fprintf(stderr, "ERROR in line %d of BED file.\n", lineno);
        exit(1);
      }
    }

    lst_free_strings(l);
  }

  str_free(line);
  str_free(source);
  lst_free(l);
  lst_free(block_sizes);
  lst_free(block_starts);
}

/* used in gff_print_bed, below */
void gff_print_bed_line(FILE *OUTF, List *features, 
                        int bed_start, int bed_end, char *name) {
  GFF_Feature *feat;
  int j;
  double score;

  assert(lst_size(features) > 0);
  feat = lst_get_ptr(features, 0);

  /* use sum of non-null scores */
  score = 0;
  for (j = 0; j < lst_size(features); j++) {
    feat = lst_get_ptr(features, j);
    if (!feat->score_is_null) 
      score += feat->score;
  }

  if (name == NULL) name = ".";

  fprintf(OUTF, "%s\t%d\t%d\t%s\t%.0f\t%c\t%d\t%d\t0\t%d\t", 
          feat->seqname->chars, bed_start - 1, bed_end, 
          name, score, feat->strand, 
          bed_start - 1, bed_end, lst_size(features));

  for (j = 0; j < lst_size(features); j++) {
    feat = lst_get_ptr(features, j);
    fprintf(OUTF, "%d,", feat->end - feat->start + 1);
  }    

  fprintf(OUTF, "\t");

  for (j = 0; j < lst_size(features); j++) 
    fprintf(OUTF, "%d,", ((GFF_Feature*)lst_get_ptr(features, j))->start - 
            bed_start);

  fprintf(OUTF, "\n");
}

/** Write a GFF_Set in BED format.  */
void gff_print_bed(FILE *OUTF,  /**< output stream  */
                   GFF_Set *gff, 
                                /**< Set to write */
                   char *groupby, 
                                /**< Group features according to this
                                   tag.  All members of a group will
                                   be reported on a single line.  Use
                                   NULL to supress grouping.  */
                   List *include
                                /**< Only write features of the
                                   specified types (list of String
                                   objects).  If NULL, all features
                                   will be written.  */
                   ) {
  GFF_Feature *feat;
  String *chrom;
  char strand;
  int i, j, bed_start, bed_end;
  List *keepers = lst_new_ptr(10);

  if (groupby == NULL) {
    for (i = 0; i < lst_size(gff->features); i++) {
      GFF_Feature *feat = lst_get_ptr(gff->features, i);
      if (include == NULL || str_in_list(feat->feature, include)) {
        lst_clear(keepers);
        lst_push_ptr(keepers, feat);
        gff_print_bed_line(OUTF, keepers, feat->start, feat->end, 
                           feat->attribute->chars);
      }
    }
  }

  else {
    gff_group(gff, groupby);

    for (i = 0; i < lst_size(gff->groups); i++) {
      GFF_FeatureGroup *group = lst_get_ptr(gff->groups, i);

      assert(lst_size(group->features) > 0);

      /* if 'null' group, output individually */
      if (group->name == NULL || group->name->length == 0) {
        for (j = 0; j < lst_size(group->features); j++) {
          feat = lst_get_ptr(group->features, j);
          if (include == NULL || str_in_list(feat->feature, include)) {
            lst_clear(keepers);
            lst_push_ptr(keepers, feat);
            gff_print_bed_line(OUTF, keepers, feat->start, feat->end,
                               feat->attribute->chars);
          }
        }
      }
      else {                      /* true (non-null) group */
        /* scan set for shared data */
        bed_start = -1;
        bed_end = 0;
        chrom = NULL;
        strand = '\0';
        lst_clear(keepers);
        for (j = 0; j < lst_size(group->features); j++) {
          feat = lst_get_ptr(group->features, j);
          if (include == NULL || str_in_list(feat->feature, include)) {

            /* check that shared information is the same for all features
               in group */
            if (chrom == NULL) {    /* first time */
              chrom = feat->seqname;
              strand = feat->strand;
            }
            else if (!str_equals(chrom, feat->seqname) || 
                     strand != feat->strand)
              die("ERROR: features of group '%s' have inconsistent chromosome or strand.\n", group->name->chars);

            /* take beg and end of interval to be min and max of indiv
               features */
            if (bed_start == -1 || feat->start < bed_start) 
              bed_start = feat->start;
            if (feat->end > bed_end) bed_end = feat->end;
            lst_push_ptr(keepers, feat);
          }
        }

        if (lst_size(keepers) > 0) 
          gff_print_bed_line(OUTF, keepers, bed_start, bed_end, 
                             group->name->chars);
      }
    }
  }
  lst_free(keepers);
}
