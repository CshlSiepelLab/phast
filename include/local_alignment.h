/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#ifndef LOC_ALN
#define LOC_ALN

#include "lists.h"
#include "stringsplus.h"
#include "msa.h"
#include "hashtable.h"

typedef struct {
  int query_beg, query_end, target_beg, target_end;
  /* eventually might keep track of percent id, number of mismatches,
     number of matches, as allowed (variously) in PSL, LAV formats */
 } GaplessAlignment;

typedef struct {
  int query_beg, query_end, target_beg, target_end;
  double score;                 /* optional */
  String *target_seq;
  List *gapless_alns;
} AlignmentBlock;

typedef struct{
  String *query_name, *target_name;
  String *query_seq, *target_seq;
  int query_len, target_len;    /* total length, not just aligned portions */
  List *alignment_blocks;
} LocalPwAlignment;

typedef enum {ADJUSTLEFT, ADJUSTRIGHT} adjust_dir;
                                /* defines direction in which to
                                   adjust coordinates when
                                   transforming via a pairwise
                                   alignment; see
                                   la_get_target_coord */

LocalPwAlignment *la_new();
LocalPwAlignment *la_read_lav(FILE *F, int read_seqs);
GaplessAlignment *la_new_gapless_aln(int query_beg, int query_end, 
                                     int target_beg, int target_end);
AlignmentBlock *la_new_alignment_block(int query_beg, int query_end, 
                                       int target_beg, int target_end, 
                                       double score, String *target_seq);
MSA* la_to_msa(LocalPwAlignment *lpwa, int force_global);
int la_get_target_coord(LocalPwAlignment *lpwa, int query_coord, 
                        adjust_dir adjust);
void la_gff_transform(LocalPwAlignment *lpwa, GFF_Set *gff);

#endif
