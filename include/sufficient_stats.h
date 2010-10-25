
/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* sufficient_stats.h - Representation of multiple alignments in terms of their sufficient statistics */

/* $Id: sufficient_stats.h,v 1.13 2009-02-19 19:41:03 agd27 Exp $ */

#ifndef MSA_SS_H
#define MSA_SS_H

#include "hashtable.h"
#include "lists.h"
#include "msa.h"
#include "external_libs.h"

/* sufficient statistics for an alignment */
/* for now, allow only one tuple_size per object */
struct msa_ss_struct {
  int tuple_size;               /* number of adjacent columns to
                                   consider as a 'column tuple' */
  int ntuples;                  /* number of distinct tuples */
  char **col_tuples;            /* the actual column tuples;
                                   col_tuples[i] is string of length
                                   ss->nseqs * tuple_size */
  int *tuple_idx;               /* defines order of column tuples in
                                   alignment; tuple_idx[i] is the
                                   index in col_tuples of the tuple
                                   that appears at position i of the
                                   MSA */
  double *counts;
  double **cat_counts;
  MSA *msa;                     /* parent alignment */
  int alloc_len, alloc_ntuples; /* for ss_realloc */
};

typedef struct msa_ss_struct MSA_SS; /* note: declared as an
                                        incomplete type in "msa.h" */

/* a set of multiple alignments, whose statistics are "pooled" (e.g.,
   to serve as a training set for a combined phylogenetic and hidden
   Markov model).  The (unordered) sufficient stats for the pool will
   be in pooled_msa->ss, and the (possibly ordered) sufficient stats
   for the individual alignments will remain in the 'ss' attributes of
   the source MSAs.  All MSAs (source and pooled) will use the same
   'col_tuples' array. */
typedef struct {
  MSA *pooled_msa;              /* contains pooled statistics */
  List *source_msas;            /* elements expected to be pointers to MSAs */
  int *lens;                    /* length of each alignment */
  int **tuple_idx_map;          /* mapping from tuple indices of
                                   source msas to those of pooled msa */
} PooledMSA;

void ss_from_msas(MSA *msa, int tuple_size, int store_order, 
                  List *cats_to_do, MSA *source_msa, 
                  Hashtable *existing_hash, int idx_offset,
		  int non_overlapping);
PooledMSA *ss_pooled_from_msas(List *source_msas, int tuple_size, 
                               int ncats, List *cats_to_do);
void ss_free_pooled_msa(PooledMSA *pmsa);
MSA *ss_aggregate_from_files(List *fnames, msa_format_type format,
                             List *seqnames, int tuple_size, 
                             List *cats_to_do, int cds_mode);
void ss_to_msa(MSA *msa);
char *ss_get_one_seq(MSA *msa, int spec);
void msa_read_AXT(MSA *msa, List *axt_fnames);
void ss_write(MSA *msa, FILE *F, int show_order);
MSA* ss_read(FILE *F, char *alphabet);
void ss_free_categories(MSA_SS *ss);
void ss_free(MSA_SS *ss);
void ss_update_categories(MSA *msa);
void ss_new(MSA *msa, int tuple_size, int max_ntuples, int do_cats, 
            int store_order);
void ss_realloc(MSA *msa, int tuple_size, int max_ntuples, int do_cats, 
                int store_order);
void ss_compact(MSA_SS *ss);
MSA* ss_alt_msa(MSA *orig_msa, int new_tuple_size, int store_order, 
                int col_offset);
MSA *ss_sub_alignment(MSA *msa, char **new_names, List *include_list, 
                      int start_col, int end_col);
void ss_reverse_compl(MSA *msa);
void ss_reorder_rows(MSA *msa, int *new_to_old, int new_nseqs);
void ss_remove_zero_counts(MSA *msa);
void ss_unique(MSA *msa);
void ss_collapse_missing(MSA *msa, int do_gaps);
int ss_seqlen(MSA *msa, int seqidx);
void ss_strip_gaps(MSA *msa, int gap_strip_mode);
void ss_strip_missing(MSA *msa, int refseq);
int ss_is_4d(MSA *msa, int tuple);
void ss_reduce_tuple_size(MSA *msa, int new_tuple_size);
void ss_add_seq(MSA *msa, int new_nseq);
int ss_lookup_coltuple(char *coltuple_str, Hashtable *tuple_hash, MSA *msa);
void ss_add_coltuple(char *coltuple_str, void *val, Hashtable *tuple_hash, MSA *msa);
void ss_make_ordered(MSA *msa);

/* Produce a string representation of an alignment column tuple, given
   the model order; str must be allocated externally to size
   msa->nseqs * (tuple_size) + 1 */
/* Note: this has been redefined so that any sequence which contains
   only values of msa->missing[0] are removed if there are not any
   sequences with greater indices with non-missing data.  This enables
   the string to be the same regardless of whether sequences containing
   only missing data have been added to the msa yet or not */
static PHAST_INLINE 
void col_to_string(char *str, MSA *msa, int col, int tuple_size) {
  int col_offset, j, pos = 0;
  for (j = 0; j < msa->nseqs; j++) 
    for (col_offset = -1 * (tuple_size-1); col_offset <= 0; col_offset++)
      str[pos++] = 
	(col + col_offset >=0 ? msa->seqs[j][col+col_offset] : GAP_CHAR);
  str[pos] = '\0';
  if (pos != msa->nseqs * tuple_size)
    die("ERROR col_to_string: pos (%i) != msa->nseqs*tuple_size (%i*%i=%i)\n",
	pos, msa->nseqs, tuple_size, msa->nseqs*tuple_size);
}

/* Given a string representation of a column tuple, return the
   character corresponding to the specified sequence and column.  Here
   sequence indexing begins with 0.  The index 'col_offset' is defined
   relative to the last column in the tuple.  Specifically, if
   col_offset == 0 then the last column is considered, if col_offset
   == -1 then the preceding one is considered, and so on. */
static PHAST_INLINE 
char col_string_to_char(MSA *msa, char *str, int seqidx, int tuple_size, int col_offset) {
  return str[tuple_size*seqidx + tuple_size - 1 + col_offset];
                                /* FIXME: WRONG */
}

static PHAST_INLINE
void set_col_char_in_string(MSA *msa, char *str, int seqidx, int tuple_size, int col_offset, char c) {
  str[tuple_size*seqidx + tuple_size - 1 + col_offset] = c;
}

/* return character for specified sequence given index of tuple */
static PHAST_INLINE
char ss_get_char_tuple(MSA *msa, int tupleidx, int seqidx, 
                       int col_offset) {
  return col_string_to_char(msa, msa->ss->col_tuples[tupleidx], seqidx, 
                            msa->ss->tuple_size, col_offset);
}

/* return character for specified sequence at specified alignment
   position; requires representation of column order */
static PHAST_INLINE
char ss_get_char_pos(MSA *msa, int position, int seqidx,
                     int col_offset) {
  if (msa->ss->tuple_idx == NULL)
    die("ERROR ss_get_char_pos: msa->ss->tuple_idx is NULL\n");
  return col_string_to_char(msa, 
                            msa->ss->col_tuples[msa->ss->tuple_idx[position]], 
                            seqidx, msa->ss->tuple_size, col_offset);
}

/* Produce a printable representation of the specified tuple.  Strings
   representing each column will be separated by spaces */
static PHAST_INLINE
void tuple_to_string_pretty(char *str, MSA *msa, int tupleidx) {
  int stridx = 0, offset, j;
  for (offset = -1 * (msa->ss->tuple_size-1); offset <= 0; offset++) {
    for (j = 0; j < msa->nseqs; j++) {
      str[stridx++] = col_string_to_char(msa, msa->ss->col_tuples[tupleidx], 
                                         j, msa->ss->tuple_size, offset);
    }
    if (offset < 0) str[stridx++] = ' ';
  }
  str[stridx] = '\0';
}

/* fill out 'tuplestr' with tuple of characters present in specified
   sequence and column tuple; tuplestr must be allocated to at least
   tuple_size (null terminator will not be added, but if present will
   be left unchanged). */
static PHAST_INLINE
void ss_get_tuple_of_chars(MSA *msa, int tupleidx, int seqidx,
                           char *tuplestr) {
  int offset;
  for (offset = -1 * (msa->ss->tuple_size-1); offset <= 0; offset++) {
    tuplestr[msa->ss->tuple_size + offset - 1] =
      col_string_to_char(msa, msa->ss->col_tuples[tupleidx],
                         seqidx, msa->ss->tuple_size, offset);
  }
}


#endif
