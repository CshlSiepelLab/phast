/* $Id: msa.h,v 1.7 2004-06-25 07:58:37 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

/** \file msa.h
   Multiple sequence alignments.
   Reading and writing are supported in a few common formats
   (currently FASTA, PHYLIP, and the format used in the "Zoo"
   project), and a sort of "grab bag" of auxiliary functionality is
   provided, including extraction of sub-alignments (by sequence or by
   column), stripping of columns with gaps (columns having either all
   gaps or any gaps), reporting of simple statistics on gap content,
   converting between coordinate frames of the alignment and of
   individual sequences, and computing simple measures of per-column
   sequence divergence (entropy and all-pairs identity; can be done
   directly or with a sliding window). 
   \ingroup msa
*/

/* To do:

      - Although code is provided to convert between coordinate
      systems, most functions use the coordinates of the entire
      alignment.  It might be useful to generalize each routine to use
      either the full-MSA coordinates or the coordinates of any
      individual sequence in the alignment.

      - There's some awkwardness throughout with indexing, stemming
      from the use of indices that start with 0 for storage, and the
      need to produce output using indices that start with 1.  The
      conversion between the two conventions should be made more
      consistent throughout.

      - The code that claims to read PHYLIP format does not actually
      address all of the idiosyncracies of that format.  For example,
      it will not allow for the possibility that a sequence name is
      not supported by whitespace from the sequence, and it will not
      read "interleaved" sequences.

      - I have a few scripts for converting the output of
      report_gap_stats into a more directly usable form.  I should
      clean them up and make them part of this package.

      - Need option to read and write "categories" in PAML format

      - Would a "view" of an alignment be useful?  Would allow storage
      only once, multiple ways of accessing.  Could be used eventually
      for fitting models to categories (or sets of categories), for
      extraction of coding regions

      - Should support MSF format (compatibility with EMBOSS).

      - Properly scale entropy.  Also, better handling of gaps in
      entropy, PW identity

      - Should really be built on top of new List and String objects,
      for auto-memory management, easy "appends" (of sequences or
      columns).  This will require almost a complete rewrite, however.

      - Should automatically recognize alignment format when reading
      (not hard with new regex code).
*/

#ifndef MSA_H
#define MSA_H

#include <stdio.h>
#include "gff.h"
#include "category_map.h"
#include "markov_matrix.h"

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
  char **names;
  char **seqs;
  int *categories;
  int *col_types[MSA_MAX_ORDER];
  int ncol_types[MSA_MAX_ORDER];
  struct msa_ss_struct *ss;
  int ncats;
  int alloc_len;
  int idx_offset;
} MSA;

#define NCHARS 256
#define MAX_NAME_LEN 256
#define MAX_LINE_LEN 10000
#define OUTPUT_LINE_LEN 70
/** Default alphabet, assumed throughout PHAST */
#define DEFAULT_ALPHABET "ACGTN"
/** Gap character, assumed throughout PHAST */
#define GAP_CHAR '-'
/** Missing data character, assumed throughout PHAST */
#define MDATA_CHAR 'N'

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
void msa_free(MSA *msa);
void strip_gaps(MSA *msa, int gap_strip_mode);
void project(MSA *msa, int refseq);
MSA* msa_sub_alignment(MSA *msa, List *seqlist, int include, int start_col, 
                       int end_col);
void report_gap_stats(MSA *msa, char *fname_root);
msa_coord_map* msa_build_coord_map(MSA *msa, int refseq);
void msa_coord_map_print(FILE *F, msa_coord_map *map);
int msa_map_seq_to_msa(msa_coord_map *map, int seq_pos);
int msa_map_msa_to_seq(msa_coord_map *map, int pos);
void msa_map_free(msa_coord_map *map);
void msa_compute_entropy_per_column(MSA *msa, double *entropy);
void msa_compute_pwid_per_column(MSA *msa, double *ave_pw_id);
void msa_apply_sliding_window(int window_size, double *col_scores, int len,
                              double *smoothed_scores);

void msa_label_categories(MSA *msa, GFF_Set *gff, CategoryMap *cm);
int msa_get_seq_idx(MSA *msa, String *name);
void msa_map_gff_coords(MSA *msa, GFF_Set *set, int from_seq, int to_seq, 
                        int offset, CategoryMap *cm);
int msa_map_seq_to_seq(msa_coord_map *from_map, msa_coord_map *to_map, 
                       int coord);

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

int msa_read_category_labels(MSA *msa, FILE *F);

void msa_partition_by_category(MSA *msa, List *submsas, List *cats_to_do, 
                               int tuple_size);


void msa_print_stats(MSA *msa, FILE *F, char *label, int header, int start,
                     int end);
gsl_vector *msa_get_base_freqs(MSA *msa, int start, int end);
void msa_get_base_freqs_tuples(MSA *msa, gsl_vector *freqs, int k, int cat);
int msa_num_gapped_cols(MSA *msa, int gap_strip_mode, int start, int end);
unsigned int msa_ninformative_sites(MSA *msa, int cat);
void msa_index_cols(MSA *msa, int order);
String *msa_read_seq_fasta(FILE *F);
int msa_coding_clean(MSA *msa, int refseq, int min_ncodons, 
                     String *errstr);
MSA *msa_concat_from_files(List *msa_fname_list, msa_format_type input_format, 
                           List *seqnames, char *alphabet);
void msa_concatenate(MSA *aggregate_msa, MSA *source_msa);
void msa_indel_clean(MSA *msa, int indel_border, int min_nbases, 
                     int min_nseqs, int tuple_size, char mdata_char);
void msa_permute(MSA *msa);
void msa_scores_as_samples(MSA *msa, FILE *F, double *scores, 
                           char *chrom, char *name,
                           double mult_fact, double threshold, int refseq, 
                           int coord_offset);
void msa_reorder_rows(MSA *msa, List *target_order);
char msa_get_char(MSA *msa, int seq, int pos);
msa_format_type msa_str_to_format(char *str);
msa_format_type msa_format_for_suffix(char *fname);
char *msa_suffix_for_format(msa_format_type t);
void msa_remove_N_from_alph(MSA *msa);
int msa_all_gaps_but_ref(MSA *msa, int pos, int ref);
void msa_find_noaln(MSA *msa, int refseqidx, int min_block_size, int *noaln);
List *msa_seq_indices(MSA *msa, List *seqnames);

#endif
