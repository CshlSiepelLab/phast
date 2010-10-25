/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: msa.h,v 1.23 2009-01-09 22:01:00 mt269 Exp $ */

/** \file msa.h
   Multiple sequence alignments.
   \ingroup msa
*/

#ifndef MSA_H
#define MSA_H

#include <stdio.h>
#include <gff.h>
#include <category_map.h>
#include <markov_matrix.h>
#include <vector.h>
#include <hashtable.h>

/* maximum "order" allowed when indexing columns (see
   msa_index_cols) */
#define MSA_MAX_ORDER 5         

struct msa_ss_struct;           /* using an incomplete type in lieu of
                                   including "sufficient_stats.h"
                                   avoids problems caused by
                                   reciprocal dependencies in header
                                   files */

/** Multiple sequence alignment object */
typedef struct {
  int nseqs;                    /**< Number of sequences */
  unsigned int length;          /**< Number of columns */
  char *alphabet;               /**< Alphabet (see #DEFAULT_ALPHABET) */
  int inv_alphabet[NCHARS];
  char *missing;                /**< Recognized missing data characters */
  int is_missing[NCHARS];       /**< Fast lookup of whether character
                                   is missing data char */
  char **names;
  char **seqs;
  int *categories;
  struct msa_ss_struct *ss;
  int ncats;
  int alloc_len;
  int idx_offset;
  int *is_informative;          /**< If non-NULL, indicates which
                                   sequences are to be considered
                                   "informative", e.g., for
                                   phylogenetic analysis */
} MSA;

#define NCHARS 256
#define MAX_LINE_LEN 10000
#define OUTPUT_LINE_LEN 70
/** Default alphabet, assumed throughout PHAST */
#define DEFAULT_ALPHABET "ACGT"
/** Gap character, assumed throughout PHAST */
#define GAP_CHAR '-'
/** Default missing data character */
#define DEFAULT_MDATA_CHARS "*N"

/** Format types */
typedef enum {PHYLIP,           /**< PHYLIP format  */
              MPM,              /**< Format used by MultiPipMaker and
                                   some of Webb Miller's older
                                   tools  */
              FASTA,            /**< Standard FASTA format */
              SS,               /**< "Sufficient statistics" format,
                                   in which each unique alignment
                                   column (or tuple of columns) is
                                   represented only once, and a count
                                   is maintained of how many times it
                                   occurs */
              LAV,              /**< lav format, used by BLASTZ */
              MAF               /**< Multiple Alignment Format (MAF)
                                   used by MULTIZ and TBA  */
} msa_format_type;

#define STRIP_ALL_GAPS -2 
#define STRIP_ANY_GAPS -1
#define NO_STRIP 0


/** Coordinate map, defined by a sequence/alignment pair.  Allows fast
    conversion between the coordinate frame of a multiple alignment
    and the coordinate frame of one of the sequences in the
    alignment.  */
typedef struct {
  List *seq_list;               /**< list of indexes in sequence
                                   immediately following gaps,
                                   expressed in the coordinate frame
                                   of the sequence (starting with
                                   position 1) */
  List *msa_list;               /**< list of corresponding indices in
                                   the frame of the MSA */
  int seq_len;                  /**< length of sequence */
  int msa_len;                  /**< length of alignment */
} msa_coord_map;


MSA *msa_new(char **seqs, char **names, int nseqs, int length, 
             char *alphabet);
MSA *msa_new_from_file(FILE *F, msa_format_type format, char *alphabet);
MSA *msa_create_copy(MSA *msa, int suff_stats_only);
MSA *msa_read_fasta(FILE *F, char *alphabet);
void msa_print(FILE *F, MSA *msa, msa_format_type format, int pretty_print);
void msa_print_to_file(const char *filename, MSA *msa, msa_format_type format, 
		       int pretty_print);
void msa_free_categories(MSA *msa);
void msa_free_seqs(MSA *msa);
void msa_free(MSA *msa);
void reduce_to_4d(MSA *msa, CategoryMap *cm);
void msa_update_length(MSA *msa);
void msa_strip_gaps(MSA *msa, int gap_strip_mode);
void msa_project(MSA *msa, int refseq);
MSA* msa_sub_alignment(MSA *msa, List *seqlist, int include, int start_col, 
                       int end_col);
msa_coord_map* msa_build_coord_map(MSA *msa, int refseq);
void msa_coord_map_print(FILE *F, msa_coord_map *map);
int msa_map_seq_to_msa(msa_coord_map *map, int seq_pos);
int msa_map_msa_to_seq(msa_coord_map *map, int pos);
void msa_map_free(msa_coord_map *map);
void msa_label_categories(MSA *msa, GFF_Set *gff, CategoryMap *cm);
int msa_get_seq_idx(MSA *msa, const char *name);
void msa_map_gff_coords(MSA *msa, GFF_Set *set, int from_seq, int to_seq, 
                        int offset, CategoryMap *cm);
int msa_map_seq_to_seq(msa_coord_map *from_map, msa_coord_map *to_map, 
                       int coord);

int msa_add_seq(MSA *msa, char *name);
void msa_add_seq_ss(MSA *msa, int new_nseq);

/* Returns complement of a single DNA base.  Leaves all bases but A,C,G,T unchanged. */
char msa_compl_char(char c);

/* Reverse complements a DNA sequence. */
void msa_reverse_compl_seq(char *seq, int length);

/* Reverse complements a segment of a sequence */
void msa_reverse_compl_seq_segment(char *seq, int start, int end);

/* Reverse complements an entire MSA. */
void msa_reverse_compl(MSA *msa);

/* Reverse complements a segment of an MSA. */
void msa_reverse_compl_segment(MSA *msa, int start, int end);

/* Reverse complements all segments of an MSA corresponding to "groups"
   in a GFF that appear to be completely on the reverse strand.  The
   GFF is partitioned into "groups" using gff_partition_by_group, and
   groups are tested using gff_reverse_strand_only (see docs for these
   functions). */
void msa_reverse_compl_gff(MSA *msa, GFF_Set *gff, int *aux_data);

void msa_reverse_compl_feats(MSA *msa, GFF_Set *feats, int *aux_data);

void msa_partition_by_category(MSA *msa, List *submsas, List *cats_to_do, 
                               int tuple_size);


void msa_print_stats(MSA *msa, FILE *F, char *label, int header, int start,
                     int end);
Vector *msa_get_base_freqs(MSA *msa, int start, int end);
void msa_get_base_freqs_tuples(MSA *msa, Vector *freqs, int k, int cat);
int msa_seqlen(MSA *msa, int seqidx);
int msa_num_gapped_cols(MSA *msa, int gap_strip_mode, int start, int end);
unsigned int msa_ninformative_sites(MSA *msa, int cat);
GFF_Set *msa_get_informative_feats(MSA *msa, int min_informative,
				   List *speclist, int refseq, 
				   int gaps_are_informative);
void msa_index_cols(MSA *msa, int order);
String *msa_read_seq_fasta(FILE *F);
int msa_codon_clean(MSA *msa, const char *refseq, char strand);
int msa_coding_clean(MSA *msa, int refseq, int min_ncodons, 
                     String *errstr, int keep_stop_codons);
MSA *msa_concat_from_files(List *msa_fname_list, msa_format_type input_format, 
                           List *seqnames, char *alphabet);
void msa_concatenate(MSA *aggregate_msa, MSA *source_msa);
void msa_indel_clean(MSA *msa, int indel_border, int min_nbases, 
                     int min_nseqs, int tuple_size, char mdata_char);
void msa_permute(MSA *msa);
void msa_reorder_rows(MSA *msa, List *target_order);
char msa_get_char(MSA *msa, int seq, int pos);
msa_format_type msa_str_to_format(const char *str);
msa_format_type msa_format_for_suffix(char *fname);
char *msa_suffix_for_format(msa_format_type t);
void msa_remove_N_from_alph(MSA *msa);
void msa_find_noaln(MSA *msa, int refseqidx, int min_block_size, int *noaln);
int msa_missing_col(MSA *msa, int ref, int pos);
List *msa_seq_indices(MSA *msa, List *seqnames);
void msa_mask_macro_indels(MSA *msa, int k, int refseq);
void msa_set_informative(MSA *msa, List *not_informative);
void msa_reset_alphabet(MSA *msa, char *newalph);
void msa_missing_to_gaps(MSA *msa, int refseq);
int msa_alph_has_lowercase(MSA *msa);
void msa_toupper(MSA *msa);
void msa_delete_cols(MSA *msa, int *delete_cols);
void msa_realloc(MSA *msa, int new_length, int new_alloclen, int do_cats, int store_order);
double msa_fraction_pairwise_diff(MSA *msa, int idx1, int idx2, 
				  int ignore_missing, int ignore_gaps);

#endif
