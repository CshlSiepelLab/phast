/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: local_alignment.c,v 1.2 2008-11-12 02:07:59 acs Exp $

   Functions for dealing with pairwise local alignments, primarily as
   produced by BLASTZ.  This code is somewhat experimental.  '

   Throughout, "query" is assumed to be the "reference" sequence, and
   "target" the sequence that is aligned to it.  Global-like alignment
   representations will have the entire query sequence but only
   aligned portions of the target sequence (gaps will appear at
   unaligned positions.

   Loose ends:

    - what about strand info?  Automatically reverse complement?

    - should be able to handle formats such as AXT which have the
    actual sequence as well (in these cases, the secondary sequence is
    never read in).  (version of PSL as well?)
*/



#include "local_alignment.h"
#include "misc.h"
#include <math.h>

LocalPwAlignment *la_new() {
  LocalPwAlignment *lpwa = (LocalPwAlignment*)smalloc(sizeof(LocalPwAlignment));
  lpwa->query_name = lpwa->target_name = lpwa->query_seq = 
    lpwa->target_seq = NULL;
  lpwa->query_len = lpwa->target_len = -1;
  lpwa->alignment_blocks = lst_new_ptr(10);
  return lpwa;
}

LocalPwAlignment *la_read_lav(FILE *F, int read_seqs) {
  String *line = str_new(STR_MED_LEN);
  int line_no=0;
  LocalPwAlignment *lpwa = la_new();
  List *fields = lst_new_ptr(6);
  Regex *stanza_start_re = str_re_new("^([dshaxm])[[:space:]]*{");
  AlignmentBlock *aln_block = NULL;
  char stanza_type = '\0';
  int i;
  int done_with[256];
  done_with[(int)'d'] = done_with[(int)'s'] = done_with[(int)'h'] = 
    done_with[(int)'x'] = done_with[(int)'m'] = 0;

  while (str_readline(line, F) != EOF) {
    str_trim(line);
    if (line->length == 0) continue;

    checkInterruptN(line_no, 1000);
    line_no++;
    if (line_no == 1) {
      if (!str_equals_charstr(line, "#:lav")) {
        die("ERROR: lav file missing header.\n");
      }
    } 
    else if (str_re_match(line, stanza_start_re, fields, 1) >= 0) {
      String *tmpstr = lst_get_ptr(fields, 1);
      stanza_type = tmpstr->chars[0];
      str_free(tmpstr);
      str_free(lst_get_ptr(fields, 0));

      if (stanza_type != 'a' && done_with[(int)stanza_type]) {
        die("ERROR: multiple '%c' stanzas in lav file.\n", 
                stanza_type);
      }

      if (stanza_type == 'a') {
        aln_block = la_new_alignment_block(-1, -1, -1, -1, -1, NULL);
        lst_push_ptr(lpwa->alignment_blocks, aln_block);
      }
    }

    /* end current stanza */
    else if (str_equals_charstr(line, "}")) {
      if (stanza_type == '\0') {
        die("ERROR: end stanza without matching begin.\n");
      }
      done_with[(int)stanza_type] = 1;
      stanza_type = '\0';
    }

    else if (stanza_type == 'd') {
      ; /* do nothing for now */
    }
    else if (stanza_type == 's') {
      int beg, end;
      String *tmpstr, *fname, *seq=NULL;
      FILE *F2;

      str_double_trim(line);
      str_split(line, NULL, fields);
      if (lst_size(fields) != 3 || 
          str_as_int(lst_get_ptr(fields, 1), &beg) != 0 || 
          str_as_int(lst_get_ptr(fields, 2), &end) != 0) {
        die("ERROR: bad line in 's' stanza in lav file.\n");
      }
      tmpstr = lst_get_ptr(fields, 0);
      fname = str_new(tmpstr->length-2); /* remove quotes */
      str_substring(fname, tmpstr, 1, tmpstr->length-2);
      if (read_seqs) {
        if ((F2 = fopen(fname->chars, "r")) == NULL) {
          die("ERROR: cannot read sequence from %s.\n", fname->chars);
        }
        seq = msa_read_seq_fasta(F2);
      }

      for (i = 0; i < lst_size(fields); i++) 
        str_free(lst_get_ptr(fields, i));
      
      if (beg != 1) {
        die("ERROR: unexpected begin index in 's' stanza of lav file (begin index currently must be 1).\n");
      }
      if (lpwa->query_len == -1) {
        lpwa->query_len = end;
        if (read_seqs) lpwa->query_seq = seq;
      }
      else if (lpwa->target_len == -1) {
        lpwa->target_len = end;
        if (read_seqs) lpwa->target_seq = seq;
      }
      else {
        die("ERROR: too many sequences listed in 's' stanza of lav file.\n");
      }
      str_free(fname);
    }
    else if (stanza_type == 'h') {
      String *name;

      str_double_trim(line);
      name = str_new(line->length-3); /* get rid of quotes and
                                         leading '>' */
      str_substring(name, line, 2, line->length-3);
      if (lpwa->query_name == NULL) lpwa->query_name = name;
      else if (lpwa->target_name == NULL) lpwa->target_name = name;
      else {
        die("ERROR: too many entries in 'h' stanza of lav file.\n");
      }
    }
    else if (stanza_type == 'a') {
      String *type;
      int val[6];
      if (!done_with[(int)'s'] || !done_with[(int)'d'] || 
          !done_with[(int)'h']) {
        die("ERROR: 'a' stanza appears in lav file before 'd', 's', or 'h' stanza.\n");
      }

      str_double_trim(line);
      str_split(line, NULL, fields);
      type = lst_get_ptr(fields, 0);
      if (lst_size(fields) > 6) {
        die("ERROR: illegal line in 'a' stanza.\n");
      }
      for (i = 1; i < lst_size(fields); i++) {
        str_as_int(lst_get_ptr(fields, i), &val[i]);
        str_free(lst_get_ptr(fields, i));
      }

      if (type->chars[0] == 's') 
        aln_block->score = val[1];
      else if (type->chars[0] == 'b') {
        aln_block->query_beg = val[1];
        aln_block->target_beg = val[2];      
      }
      else if (type->chars[0] == 'e') {
        aln_block->query_end = val[1];
        aln_block->target_end = val[2];      
      }
      else if (type->chars[0] == 'l') 
        lst_push_ptr(aln_block->gapless_alns, 
                     la_new_gapless_aln(val[1], val[3], val[2], val[4]));

      str_free(type);
    }     
  }
  str_free(line);
  lst_free(fields);
  str_re_free(stanza_start_re);

  return lpwa;
}

LocalPwAlignment* la_read_psl() {
  return NULL;
}

LocalPwAlignment* la_read_axt() {
  return NULL;
}

GaplessAlignment *la_new_gapless_aln(int query_beg, int query_end, 
                                     int target_beg, int target_end) {
  GaplessAlignment *retval = 
    (GaplessAlignment*)smalloc(sizeof(GaplessAlignment));
  retval->query_beg = query_beg;
  retval->query_end = query_end;
  retval->target_beg = target_beg;
  retval->target_end = target_end;
  return retval;  
}

/* note: uses same target_seq object (does not copy) */
AlignmentBlock *la_new_alignment_block(int query_beg, int query_end, 
                                       int target_beg, int target_end, 
                                       double score, String *target_seq) {
  AlignmentBlock *retval = 
    (AlignmentBlock*)smalloc(sizeof(AlignmentBlock));
  retval->query_beg = query_beg;
  retval->query_end = query_end;
  retval->target_beg = target_beg;
  retval->target_end = target_end;
  retval->score = score;
  retval->target_seq = target_seq;
  retval->gapless_alns = lst_new_ptr(10);
  return retval;  
}

/* String *la_make_target_seq_string(AlignmentBlock *block,  */
/*                                   String *whole_target_seq) { */
/*   String *retval = str_new(??); */
/*   GaplessAlignment *ga, *last_ga = NULL; */
/*   int i, j; */
/*   for (i = 0; i < lst_size(block->gapless_alns); i++) { */
/*     ga = lst_get_ptr(block->gapless_alns); */
/*     if (last_ga != NULL) */
/*       for (j = last_ga->target_end + 1; j < ga->target_beg; j++) */
/*         str_append_char(retval, GAP_CHAR); */
/*     str_append_substring(retval, ga->target_beg,  */
/*                          ga->target_end - ga->target_beg + 1); */
/*     last_ga = ga; */
/*   } */
/*   return retval; */
/* } */

/* Also want to be able to represent *only* as sufficient stats */
/* Want easy transition to tiled multiple alignment; can take a set of
   local alignments linked to same reference seq and make a mult
   alignment ... but can't read all of the local alignments into
   memory at once */

/* represent a local pairwise alignment as a global alignment, with
   gaps in the target at all unaligned positions */
/* assumes gapless alignments are in order wrt coords of query sequence */
/* currently pretty inefficient */
/* should directly create a sufficient statistics representation */
/* force_global option causes all of target seq to be represented
   also.  Only makes sense for an alignment in which order is
   maintained in both seqs */ 
/* if query_seq and target_seq are NULL, they will attempt to be read
   from the filenames contained in the LocalPwAlignment object (FASTA
   format is assumed) */
MSA* la_to_msa(LocalPwAlignment *lpwa, int force_global) {
  int i, j, k, len;
  char **names = (char**)smalloc(2 * sizeof(char*));
  char **seqs = (char**)smalloc(2 * sizeof(char*));
  String *query_seq = lpwa->query_seq, *target_seq = lpwa->target_seq;
  String *qseq = str_new(query_seq->length);
  String *tseq = str_new(target_seq->length);
  GaplessAlignment *lga = NULL;

  names[0] = copy_charstr(lpwa->query_name->chars);
  names[1] = copy_charstr(lpwa->target_name->chars);

  for (i = 0; i < lst_size(lpwa->alignment_blocks); i++) {
    AlignmentBlock* b = lst_get_ptr(lpwa->alignment_blocks, i);
    checkInterrupt();
    for (j = 0; j < lst_size(b->gapless_alns); j++) {
      GaplessAlignment *ga = lst_get_ptr(b->gapless_alns, j);

      if (lga == NULL) {
        for (k = 0; k < ga->query_beg-1; k++) {
          str_append_char(qseq, query_seq->chars[k]);
          str_append_char(tseq, GAP_CHAR);
        }
        if (force_global) {
          for (k = 0; k < ga->target_beg-1; k++) {
            str_append_char(qseq, GAP_CHAR);
            str_append_char(tseq, target_seq->chars[k]);
          }
        }
      }
      else {
        if (lga->query_end >= ga->query_beg ||
            (force_global && lga->target_end >= ga->target_beg)) {
          die("ERROR: overlapping alignment segments.\n");
        }
        if (j > 0 && lga->query_end == ga->query_beg-1) {
                                /* gap in query seq */
          for (k = lga->target_end; k < ga->target_beg-1; k++) {
            str_append_char(qseq, GAP_CHAR);
            str_append_char(tseq, target_seq->chars[k]);
          }
        }
        else {                  /* gap in target seq */
          for (k = lga->query_end; k < ga->query_beg-1; k++) {
            str_append_char(qseq, query_seq->chars[k]);
            str_append_char(tseq, GAP_CHAR);
          }
          if (force_global) {
            for (k = lga->target_end; k < ga->target_beg-1; k++) {
              str_append_char(qseq, GAP_CHAR);
              str_append_char(tseq, target_seq->chars[k]);
            }
          }
        }
      }

      for (k = 0; k < ga->query_end - ga->query_beg + 1; k++) {
        str_append_char(qseq, query_seq->chars[ga->query_beg + k - 1]);
        str_append_char(tseq, target_seq->chars[ga->target_beg + k - 1]);
      }

      lga = ga;
    }
  }

  for (k = lga->query_end; k < query_seq->length; k++) {
    str_append_char(qseq, query_seq->chars[k]);
    str_append_char(tseq, GAP_CHAR);
  }

  if (force_global) {
    for (k = lga->target_end; k < target_seq->length; k++) {
      str_append_char(qseq, GAP_CHAR);
      str_append_char(tseq, target_seq->chars[k]);
    }
  }

  seqs[0] = qseq->chars;
  seqs[1] = tseq->chars;
  qseq->chars = NULL; tseq->chars = NULL;
  len = qseq->length;
  str_free(qseq); str_free(tseq);
  return msa_new(seqs, names, 2, len, NULL);  
}


/* build a multiple alignment from a reference sequence and a set of
   local pairwise alignments, all of which have the reference sequence
   as their query sequence.  The multiple alignment will be a
   projection on the reference sequence -- that is, that sequence will
   contain no gaps.  It will be represented in terms of its sufficient
   statistics, with ordering information optional */
MSA* la_create_multiple_pw_alignment(int keep_order) {
  return NULL;
}


/* based on a local pairwise alignment object, estimate the coordinate
   in the target sequence corresponding to the specified coordinate in
   the query sequence.  Currently assumes alignment blocks, and
   gapless alignments within them, are in sorted order (wrt query
   sequence).  Also currently does linear search (should use binary
   search).  A value of -1 is returned if no reasonable estimate is
   possible.  This routine should only be used with relatively close,
   generally orthologous sequences, having good synteny. */
int la_get_target_coord(LocalPwAlignment *lpwa, int query_coord, 
                        adjust_dir adjust) {
  int i, j;
  int q1 = -1, q2 = -2, t1 = -1, t2 = -1;
  AlignmentBlock *last_ab = NULL;
  GaplessAlignment *last_ga = NULL;
  /* find query coords q1 and q2 bracketing the position in question,
     and corresponding target coords t1 and t2 */
  for (i = 0; i < lst_size(lpwa->alignment_blocks); i++) {
    AlignmentBlock *ab = lst_get_ptr(lpwa->alignment_blocks, i);
    checkInterrupt();
    if (!(last_ab == NULL || last_ab->query_end < query_coord))
      die("ERROR la_get_target_coord: bad value for last_ab\n");
    if (ab->query_beg > query_coord) {
      /* coord falls between alignment blocks */
      if (last_ab == NULL) {    /* occurs at beginning */
        q1 = t1 = 0;
        q2 = ab->query_beg;
        t2 = ab->query_end;
      }
      else {
        q1 = last_ab->query_end;  
        q2 = ab->query_beg;
        t1 = last_ab->target_end;
        t2 = ab->target_beg;
      }
      break;
    }
    else if (ab->query_end >= query_coord) {
                                /* coord falls within an alignment
                                   block; need to look at the gapless
                                   alignments  */

      for (j = 0; j < lst_size(ab->gapless_alns); j++) {
        GaplessAlignment *ga = lst_get_ptr(ab->gapless_alns, j);
        if (!(last_ga == NULL || last_ga->query_end < query_coord))
	  die("ERROR la_get_target_coord: bad value for last_ga\n");
        if (ga->query_beg > query_coord) {
          q1 = last_ga->query_end;  
          q2 = ga->query_beg;
          t1 = last_ga->target_end;
          t2 = ga->target_beg;          
          break;
        }
        else if (ga->query_end >= query_coord) 
          /* coord falls within gapless alignment -- this case is easy */
          return (query_coord - ga->query_beg + ga->target_beg);

        last_ga = ga;           /* keep looking */
      }
      if (q1 == -1)
	die("ERROR la_get_target_coord: bad coords\n");
                                /* coords must be assigned above;
                                   otherwise the coords for the block
                                   must have been wrong */
      break;
    }

    /* keep looking */
    last_ab = ab;
  }

  if (q1 == -1) {               /* coord must occur *beyond* all
                                   alignment blocks */
    q1 = last_ab->query_end;
    q2 = lpwa->query_len-1;
    t1 = last_ab->target_end;
    t2 = lpwa->target_len-1;
  }

  if (t2 < t1) return -1;

  return (adjust == ADJUSTRIGHT ? t2 : t1);
}


/* Transform the coordinates of all features in a GFF according to a
   local alignment.  Each feature in the original GFF will be replaced
   by zero or more features with transformed begin and end
   coordinates.  The original features are "projected" onto the
   aligned (target) sequence vis the alignment, in such a way that if
   a feature contains no aligned bases, then it will not be
   represented, and if a feature contains bases that align to multiple
   "blocks", then it will be split into several features, one for each
   block.  

   The general idea is that the new features should cover only those
   bases in the target sequence that align to bases in the query
   sequence.  Currently, however, insertions in the target sequence
   between gapless alignments of the same block are ignored, so that a
   transformed feature may contain some bases that do not directly
   align to the query sequence.  The rationale is that these
   insertions should generally be small, and should reflect
   small-scale events that do not radically disrupt the local
   properties of the sequence. */
void la_gff_transform(LocalPwAlignment *lpwa, GFF_Set *gff) {
  int i, j, k;
  int new_beg, new_end;
  List *new_features = lst_new_ptr(lst_size(gff->features));
  GFF_Feature *feat, *new_feat;

  for (i = 0; i < lst_size(gff->features); i++) {
    checkInterruptN(i, 1000);
    feat = lst_get_ptr(gff->features, i);
    for (j = 0; j < lst_size(lpwa->alignment_blocks); j++) {
                                /* this is a somewhat inefficient way
                                   to proceed, but the number of
                                   features and the number of
                                   alignment blocks is usually pretty
                                   small; will adjust strategy as
                                   needed */
      AlignmentBlock *ab = lst_get_ptr(lpwa->alignment_blocks, j);
      new_beg = new_end = -1;
      if ((ab->query_beg >= feat->start && ab->query_beg <= feat->end) ||
          (ab->query_end >= feat->start && ab->query_end <= feat->end) ||
          (feat->start >= ab->query_beg && feat->end <= ab->query_end)) {
                                /* block and feature overlap */
        if (feat->start <= ab->query_beg)
                                /* feature extends to the left of the
                                   alignment block; use beg of
                                   block */
          new_beg = ab->target_beg;
        else {                  /* ab->query_beg < feat->start */
          /* find first corresponding base within a gapless alignment */
          for (k = 0; k < lst_size(ab->gapless_alns); k++) {
            GaplessAlignment *ga = lst_get_ptr(ab->gapless_alns, k);
            if (ga->query_beg >= feat->start) {
                                /* gapless alignment overlaps the
                                   feature and the feature extends to
                                   the left (equal to or) beyond the
                                   ga; use the start of the ga */
              new_beg = ga->target_beg;
              break;
            }
            else if (ga->query_end >= feat->start) {
                                /* gapless alignment overlaps the
                                   feature and the ga extends to the
                                   left beyond the feature; use the
                                   aligned base within the ga */
              new_beg = ga->target_beg + (feat->start - ga->query_beg);
              break;
            }
          }
        }
        if (feat->end >= ab->query_end) 
                                /* feature extends to the right of the
                                   alignment block; use end of
                                   block */
          new_end = ab->target_end;
        else {
          /* find last corresponding base within a gapless alignment */
          for (k = lst_size(ab->gapless_alns)-1; k >= 0; k--) {
            GaplessAlignment *ga = lst_get_ptr(ab->gapless_alns, k);
            if (ga->query_end <= feat->end) {
                                /* gapless alignment overlaps the
                                   feature and the feature extends to
                                   the right (equal to or) beyond the
                                   ga; use the end of the ga */
              new_end = ga->target_end;
              break;
            }
            else if (ga->query_beg <= feat->end) {
                                /* gapless alignment overlaps the
                                   feature and the ga extends to the
                                   right beyond the feature; use the
                                   aligned base within the ga */
              new_end = ga->target_beg + (feat->end - ga->query_beg);
              break;
            }
          }
        }
        
        if (!(new_beg != -1 && new_end != -1))
	  die("ERROR: la_gff_transform: new_beg=%i new_end=%i\n",
	      new_beg, new_end);
/*         fprintf(stderr, "(%d, %d) -> (%d, %d)\n", feat->start, feat->end, new_beg, new_end); */
        new_feat = gff_new_feature_copy(feat);
        new_feat->start = new_beg;
        new_feat->end = new_end;
        lst_push_ptr(new_features, new_feat);
      }
    }
  }

  for (i = 0; i < lst_size(gff->features); i++)
    gff_free_feature(lst_get_ptr(gff->features, i));
  lst_free(gff->features);
  gff->features = new_features;
  gff_sort(gff);
}


