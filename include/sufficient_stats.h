/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file sufficient_stats.h 
    Representation of multiple alignments in terms of their sufficient statistics 
    @ingroup msa
*/

#ifndef MSA_SS_H
#define MSA_SS_H

#include "hashtable.h"
#include "lists.h"
#include "msa.h"
#include "external_libs.h"

/** Sufficient Statistics object for an alignment. 
  @note For now, allow only one tuple_size per object */
struct msa_ss_struct {
  int tuple_size;               /**< Number of adjacent columns to
                                   consider as a 'column tuple' */
  int ntuples;                  /**< Number of distinct tuples */
  char **col_tuples;            /**< The actual column tuples;
                                   col_tuples[i] is string of length
                                   ss->nseqs * tuple_size */
  int *tuple_idx;               /**< Defines order of column tuples in
                                   alignment; tuple_idx[i] is the
                                   index in col_tuples of the tuple
                                   that appears at position i of the
                                   MSA */
  double *counts;		/** Number of times each tuple appears
                                    within the sequence data */
  double **cat_counts;		/** Counts per category  */
  MSA *msa;                     /** Parent alignment */
  int alloc_len, alloc_ntuples; /** for ss_realloc */
};

/** Alignment sufficient statistics.
    @note Completes incomplete declaration from msa.h */
typedef struct msa_ss_struct MSA_SS; 

/** Set of multiple alignments, whose statistics are "pooled".
   @note Serves well as a training set for a combined phylogenetic and hidden
   Markov model.  
   @note Sufficient statistics for the individual alignments will remain in the 'ss'
     attributes of the source MSAs
*/
typedef struct {
  MSA *pooled_msa;              /**< Contains pooled statistics from multiple MSAs */
  List *source_msas;            /**< List of individual (unpooled) MSAs */
  int *lens;                    /**< Array of int indicating length of each MSA */
  int **tuple_idx_map;          /**< Mapping from tuple indices of
                                   source msas to those of pooled msa */
} PooledMSA;

/** \name Calculate (populate) Sufficient Statistics functions */

/**  Calculate Sufficient Statistics for an MSA.
   @param[in,out] msa Uses sequence data in this MSA to generate sufficient statistics and adds results to msa->ss
   @param[in] tuple_size Number of columns that make up a tuple  (e.g., if tuple_size==1, then each column is considered individually, and if tuple_size==2, then each is considered wrt its predecessor)
   @param[in] store_order Whether to store the order in which columns appear in the source MSAs
   @param[in] cats_to_do (Optional) List of category numbers to generate sufficient statistics for; defaults to all
   @param[in,out] source_msa (Optional) Statistics for this alignment will be *added* to the statistics of 'msa' (used for combining SS)
   @param[in,out] existing_hash (Optional) Existing running hash from previous call to this function. This is used when source_msa is non-NULL
   @param[in] idx_offset (Optional) Set to -1 when not in use. If idx_offset is
   non-negative, then it must be true that store_order == 1 and
   source_msa != NULL.  If idx_offset is non-negative, then msa->ss is
   assumed to be pre-allocated (offsets complicate reallocs)
   @param[in] non_overlapping If TRUE and tuple_size > 1, only use non_overlapping tuples in SS

   @code
   //e.g. Read in FASTA file and compute Sufficient Statistics
     char *filename = "sequences.fasta"
     char *alphabet = "ACGT";
     int tuple_size = 3; 

     //Read multiple alignment fasta file
     MSA *msa = msa_new_from_file(phast_fopen(filename, "r"), alphabet); 

     //Generate Sufficient Statistics
     //Tuple Size of 3, Preserve order
     ss_from_msas(msa, tuple_size, TRUE, NULL, NULL, NULL, -1, 0);

     //Sufficient Statistics are stored here 
     msa->ss

   @endcode

 Given a multiple alignment object, create a representation based on
   its sufficient statistics -- i.e., the distinct columns that it
   includes, the number of times each one appears, and (if store_order
   == 1), the order in which they appear.  An 'MSA_SS' object will be
   created and linked to the provided alignment object.  If the
   argument 'source_msa' is non-NULL, then the stats for the specified
   alignment will be *added* to the MSA_SS object of 'msa'.  This
   option can be used in repeated calls to create aggregate
   representations of sets of alignments (see ss_pooled_from_msas and
   ss_aggregate_from_files below).  In this case, a running hash table
   should also be passed in (existing_hash).  The argument
   'tuple_size' determines what size of column tuple to consider
   (e.g., if tuple_size==1, then each column is considered
   individually, and if tuple_size==2, then each is considered wrt its
   predecessor).  If store_order==1, then the order of the columns
   will be stored.

   If source_msa != NULL, then any sequences in the 'msa' object will
   be ignored, but 'msa' must be initialized with the appropriate
   sequence names, alphabet, and number of sequences (length will be
   adjusted as needed).  A specified source alignment must have the
   same number of sequences as 'msa' and the sequences must appear in
   corresponding order.  

  If msa->ncats > 0 and msa->categories != NULL then
   category-by-category counts will be maintained.  If in addition
   source_msas != NULL, then each source alignment will be expected to
   have non-NULL 'categories' vectors containing values of at most
   msa->ncats.  The optional list 'cats_to_do' (expected to be a list
   of integers if non-NULL) allows consideration to be limited to
   certain categories.  All of this does not apply when source
   alignments are represented by (unordered) sufficient statistics
   only, in which case category-by-category counts are maintained iff
   they exist for source alignments.  

   If non_overlapping == TRUE and tuple_size > 1, then will only
   convert non_overlapping tuples into the SS structure.  The first
   tuple will start with position 1 and be of length tuple_size, the
   second will start at tuple_size+1 and go to tuple_size*2, etc.
   (This is useful for converting coding sequence into codons).

   'idx_offset' parameter for use when storing
   order and source alignments refer to local segments of a long
   alignment.  Will be added to coordinates in source_msa when setting
   tuple_idx for msa.  Set to -1 when not in use.  If idx_offset is
   non-negative, then it must be true that store_order == 1 and
   source_msa != NULL.  If idx_offset is non-negative, then msa->ss is
   assumed to be pre-allocated (offsets complicate reallocs). 


*/
void ss_from_msas(MSA *msa, int tuple_size, int store_order, 
                  List *cats_to_do, MSA *source_msa, 
                  Hashtable *existing_hash, int idx_offset,
		  int non_overlapping);

/** Pool multiple MSAs into a single object of type PooledMSA.  
   @pre All msas have same names, nseqs, and alphabet (it uses those from the first MSA in the list) 
   @param source_msas List of MSA objects to pool together 
   @param tuple_size Number of columns that make up a tuple  (e.g., if tuple_size==1, then each column is considered individually, and if tuple_size==2, then each is considered wrt its predecessor)
   @param ncats Number of categories
   @param cats_to_do (Optional) List of category numbers to generate sufficient statistics for; defaults to all
   @param non_overlapping If TRUE, then for tuples of size > 1, only collect sufficient statistics for non-overlapping tuples (useful for codons).
   @result MSAs pooled into a single PooledMSA object
   @note The source_msas parameter will be referenced in the resulting PooledMSA object so please don't delete it
*/
PooledMSA *ss_pooled_from_msas(List *source_msas, int tuple_size, 
                               int ncats, List *cats_to_do, 
			       int non_overlapping);

/** Create an aggregate MSA from a list of MSA filenames.
   @param fnames List of filenames of MSAs
   @param format File format of the MSAs i.e. FASTA, Phylip, etc. Must all be the same
   @param seqnames List of sequences (by name) you wish to use from the MSA files. Also defines the order of sequences in the aggregate MSA
   @param tuple_size  Number of columns that make up a tuple  (e.g., if tuple_size==1, then each column is considered individually, and if tuple_size==2, then each is considered wrt its predecessor)
   @param cats_to_do (Optional) List of category numbers to generate sufficient statistics for; defaults to all
   @param cycle_size (Optional) If cycle_size > 0, site categories will be labeled 1,2,...,<cycle_size>,...,1,2,...,<cycle_size>. 
   @note Missing sequences will be replaced with missing data.  
   @note All source MSAs must share the same alphabet, and each must contain a subset of the names in 'seqnames'.      
   @note No direct representation of the source or tuple order is retained 
   
 */
MSA *ss_aggregate_from_files(List *fnames, 
                             List *seqnames, int tuple_size, 
                             List *cats_to_do, int cycle_size);

/** \} */

/** Reconstruct sequences from sufficient statistics.  
   @pre Requires length and nseqs to be correct
   @pre Requires seqs to be NULL
   @pre Requires categories to be NULL
   @param msa MSA containing Sufficient Statistics but not original sequences
   @note Only uses right-most column in tuple
   @note Does not set category labels  
   @warning This will be wrong if sufficient statistics only describe certain
   categories and a higher-order model is used 
*/
void ss_to_msa(MSA *msa);

/** Return a sequence data for a single sequence from MSA.
   @param msa MSA containing sequence data
   @param spec Index of sequence to get
   @result String of sequence data for specified sequence
*/
char *ss_get_one_seq(MSA *msa, int spec);

/** Read multiple MSAs in .axt file format into a single MSA
    @param[out] Initialized MSA
    @param axt_fnames Filenames of MSA files in .axt file format
*/
void msa_read_AXT(MSA *msa, List *axt_fnames);

/** \name Sufficient Statistics file access functions 
\{*/

/** Write MSA to file as sufficient statistics.
    @param msa MSA to save as sufficient statistics
    @param F File descriptor to save to
    @param show_order Keep track of tuple order
*/
void ss_write(MSA *msa, FILE *F, int show_order);

/** Read MSA from file as sufficient statistics.
    @param F File descriptor to read sufficient statistics from
    @param alphabet Alphabet of MSA being read in
    @result MSA reconstructed from sufficient statistics
*/
MSA* ss_read(FILE *F, char *alphabet);

/** \} */

/**  Update category count according to 'categories' attribute of MSA
   object.  
   @pre Requires ordered sufficient statistics.
   @param msa MSA containing category counts to update
   @note Will allocate and initialize cat_counts if necessary.
*/
void ss_update_categories(MSA *msa);
/** \name Sufficient Statistics allocation/cleanup functions
\{ */
/** Creates a new sufficient statistics object (SS are not calculated) and links it to the
   specified alignment.  
   @param msa MSA to attach new sufficient statistics object to
   @param tuple_size Number of columns that make up a tuple  (e.g., if tuple_size==1, then each column is considered individually, and if tuple_size==2, then each is considered wrt its predecessor)
   @param max_ntuples Allocates sufficient space for max_ntuples distinct column tuples
*/
void ss_new(MSA *msa, int tuple_size, int max_ntuples, int do_cats, 
            int store_order);
/** Expand memory allocated to Sufficient Statistics to handle additional
  sequence data or tuples.
   @param msa MSA containing Sufficient Statistics
   @param tuple_size NOT USED
   @param max_ntuples Amount of tuples the Sufficient Statistics should be able to support
   @param do_cats Whether categories are being used on MSA to re-allocate
   @param store_order Whether tuple order is being stored on MSA to re-allocate
 */
void ss_realloc(MSA *msa, int tuple_size, int max_ntuples, int do_cats, 
                int store_order);

/** Free category counts.
    @param ss Sufficient Statistics containing category counts
*/
void ss_free_categories(MSA_SS *ss);

/** Free Sufficient Statistics object.
    @param ss Sufficient Statistics object to free 
*/
void ss_free(MSA_SS *ss);

/** Free a PooledMSA object and all of its contents 
  @param pmsa PooledMSA to free
*/
void ss_free_pooled_msa(PooledMSA *pmsa);

/** \} */

/** Shrinks Sufficient Statistics to size ss->ntuples (to get rid of excess allocated memory)
    @param ss Sufficient Statistics
*/
void ss_compact(MSA_SS *ss);

/** Create copy of MSA with different tuple size.
   @param orig_msa Original MSA to be copied from
   @param new_tuple_size Tuple size for the new MSA being returned
   @param store_order Whether the MSA stores order
   @param col_offset Allows category labels to be shifted *left*
   by the specified amount
   @result Copied MSA with new tuple size
   @note Sequence data and names are shared between MSA objects but new MSA will have a new sufficient statistics object.
   @note Left shifting with 'col_offset' is is convenient when estimating
   conditional probabilities with a two-pass approach. */
MSA* ss_alt_msa(MSA *orig_msa, int new_tuple_size, int store_order, 
                int col_offset);

/** Extract a sub-alignment from an alignment, in terms of ordered
   sufficient statistics.
   @param msa Alignment to extract sub-alignment from
   @param new_names List of names for sequence that will be in the sub-alignment
   @param include_list List of indexes specifying which sequences should be included in sub-alignment
   @param start_col Only copy data in columns at or after start_col
   @param end_col Only copy data in columns before end_col
   @see msa_sub_alignment
   @note The new alignment will represent the interval [start_col, end_col).
*/
MSA *ss_sub_alignment(MSA *msa, char **new_names, List *include_list, 
                      int start_col, int end_col);

/** \name Sufficient Statistics modification functions
\{ */

/** Adjust sufficient statistics to reflect the reverse complement of
   an alignment.
   @param msa MSA containing Sufficient Statistics
   @see msa_reverse_compl 
*/
void ss_reverse_compl(MSA *msa);

/** Change sufficient statistics to reflect reordered rows of an alignment.
   @param msa MSA containing Sufficient Statistics
   @param new_to_old Array of integers mapping the new row order of the alignment to the old row order
   @param new_nseqs Number of sequences in MSA
   @see msa_reorder_rows. 
*/
void ss_reorder_rows(MSA *msa, int *new_to_old, int new_nseqs);

/** Eliminate all tuples with counts of zero.  
    @pre Main counts must reflect category counts (category counts won't be checked).  
    @pre Tuples with counts of zero are assumed not to appear in tuple_idx.
    @param msa Multiple Alignment to remove tuples with counts of zero from
*/
void ss_remove_zero_counts(MSA *msa);

/** Ensures all tuples are unique.  
    @param msa MSA containing Sufficient Statistics
    @note If non-unique tuples are found, their counts are combined and indices remapped
*/
void ss_unique(MSA *msa);

/** Convert all missing data characters to the default missing data
    character (msa->missing[0]).  
    @param msa MSA containing Sufficient Statistics
    @param do_gaps Whether to also convert gap
    characters.  
    @note Can be useful in reducing number of tuples 
*/
void ss_collapse_missing(MSA *msa, int do_gaps);

/** Strip gaps from Sufficient Statistics.
   @param msa Multiple Alignment containing Sufficient Statistics
   @param gap_strip_mode How to handle gaps; If
   gap_strip_mode is STRIP_ALL_GAPS or STRIP_ANY_GAPS, removes all
   columns with ALL or ANY gaps, respectively.  Otherwise, assumes a
   *projection* is desired onto the sequence whose index is
   gap_strip_mode (indexing starts with 1).  
   @note Changes are made to original alignment.  
   @note Gaps are expected to be represented by GAP_CHAR.  
   @note If msa->categories is non-NULL, it will be adjusted accordingly.
   @see msa_strip_gaps
*/
void ss_strip_gaps(MSA *msa, int gap_strip_mode);
/** \} */

/** Strip columns consisting of missing data in all
   sequences but the reference sequence.
   @param refseq Index of the reference sequence (starts with 1)
   @note If msa->categories is non-NULL, it will be adjusted accordingly 
*/
void ss_strip_missing(MSA *msa, int refseq);

/** Return the number of non-gap characters for a particular sequence.
    @param msa Multiple Alignment
    @param seqidx Index of sequence to use when counting non-gap characters
    @result Number of non-gap characters for specified sequence
*/
int ss_seqlen(MSA *msa, int seqidx);

/** Check if a column of is 4 fold degenerate.
  @param msa Multiple Alignment that contains Sufficient Statistics
  @param tuple Tuple index to check for 4 fold degenerate 
*/
int ss_is_4d(MSA *msa, int tuple);



/** Reduce the tuple size of an MSA.
   @param msa containing Sufficient Statistics
   @param new_tuple_size New tuple size
*/
void ss_reduce_tuple_size(MSA *msa, int new_tuple_size);

/* Not found in implementation */
void ss_add_seq(MSA *msa, int new_nseq);

/** Retrieve tuple index from hash table using tuple as key
  @param[in] coltuple_str Key used for hash table 
  @param[in] tuple_hash Hash table of tuple indexes (integers)
  @param[in] msa Alignment object
  @result -1 if not found, otherwise Value stored in tuple_hash at key coltuple_str
 */
int ss_lookup_coltuple(char *coltuple_str, Hashtable *tuple_hash, MSA *msa);

/** Add tuple index to hash table using tuple as key
  @param[in] coltuple_str Key used for hash table
  @param[in] val Value (tuple index) to add to hash table
  @param[in] tuple_hash Hash table of tuple indices
  @param[in] msa Multiple Alignment object
 */
void ss_add_coltuple(char *coltuple_str, void *val, Hashtable *tuple_hash, MSA *msa);

/** Impose an artificial ordering on tuples if they aren't already ordered.
  @param msa Multiple Alignment to order
*/
void ss_make_ordered(MSA *msa);

/** Produce a string representation of an alignment column tuple.
   @param[out] str Representation of alignment column (must be allocated externally to size
   msa->nseqs * (tuple_size) + 1)
   @param[in] msa Multiple alignment containing column to represent
   @param[in] col Column index to represent
   @param[in] tuple_size Size of the tuples to be used in str
   @warning Any sequence which contains only values of msa->missing[0] are removed if
   There are not any sequences with greater indices and non-missing data
   @note The missing data side effect enables
   the string to be the same regardless of whether sequences containing
   only missing data have been added to the msa yet or not
 */
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

/** Given a string representation of a column tuple, return the
   character corresponding to the specified sequence and column. 
  
   @param msa NOT USED
   @param str Column Tuple
   @param seqidx Sequence Index 
   @param tuple_size Size of tuples
   @param col_offset Column Offset relative to the last column in the tuple.   Specifically, if
   col_offset == 0 then the last column is considered, if col_offset
   == -1 then the preceding one is considered, and so on. 
   @note sequence indexing begins with 0.  

*/
static PHAST_INLINE 
char col_string_to_char(MSA *msa, char *str, int seqidx, int tuple_size, int col_offset) {
  return str[tuple_size*seqidx + tuple_size - 1 + col_offset];
                                /* FIXME: WRONG */
}

/** Set a single character within a column tuple.
    @param msa NOT USED
    @param Column tuple
    @param seqidx Index of sequence to modify
    @param tuple_size Size of the tuples in str
    @param col_offset Column Offset relative to the last column in the tuple. Specifically, if
   col_offset == 0 then the last column is considered, if col_offset
   == -1 then the preceding one is considered, and so on. 
   @note sequence indexing begins with 0.
    @param c Value to set */
static PHAST_INLINE
void set_col_char_in_string(MSA *msa, char *str, int seqidx, int tuple_size, int col_offset, char c) {
  str[tuple_size*seqidx + tuple_size - 1 + col_offset] = c;
}

/** \name Get character of tuple from alignment 
\{ */

/** Return a specific base from alignment given sequence, tuple, and column offset of tuple
   @param msa Multiple Alignment
   @param tupleidx Index specifying which tuple to retrieve base from
   @param seqidx Index specifying which sequence to retrieve base from
   @param col_offset Identifies which column within the tuple to retrieve base from. Specifically, if
   col_offset == 0 then the last column is considered, if col_offset
   == -1 then the preceding one is considered, and so on. 
   @result Base character at sequence seqidx, tuple tupleidx, and column offset col_offset
 */
static PHAST_INLINE
char ss_get_char_tuple(MSA *msa, int tupleidx, int seqidx, 
                       int col_offset) {
  return col_string_to_char(msa, msa->ss->col_tuples[tupleidx], seqidx, 
                            msa->ss->tuple_size, col_offset);
}

/** Return character for specified sequence at specified alignment
   position; 
  @pre Requires representation of column order 
  @param msa Multiple Alignment
  @param position Position of tuple within MSA
  @param seqidx Sequence index
  @param col_offset Identifies which column within the tuple to retrieve base from. Specifically, if
   col_offset == 0 then the last column is considered, if col_offset
   == -1 then the preceding one is considered, and so on. 
  @result Base char from tuple in alignment
*/
static PHAST_INLINE
char ss_get_char_pos(MSA *msa, int position, int seqidx,
                     int col_offset) {
  if (msa->ss->tuple_idx == NULL)
    die("ERROR ss_get_char_pos: msa->ss->tuple_idx is NULL\n");
  return col_string_to_char(msa, 
                            msa->ss->col_tuples[msa->ss->tuple_idx[position]], 
                            seqidx, msa->ss->tuple_size, col_offset);
}
/** \} */
/** Produce a printable representation of the specified tuple. 
  @param[out] str String to hold representation, must be externally allocated to size (nseqs * (tuple_size + 1) + 1 )
  @param[in] msa Multiple Alignment
  @param[in] tupleidx Tuple index to print
  @note Strings representing each column will be separated by spaces 
  i.e. 'ATA AAA AAT AAT' represents 4 sequences at specified tuple
*/
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

/** Get tuple as string specified by tuple & sequence indices 
   Fill out 'tuplestr' with tuple of characters present in specified
   sequence and column tuple; 
   @pre Tuplestr must be allocated to at least
   tuple_size (null terminator will not be added, but if present will
   be left unchanged).
   @param[in] msa Multiple Alignment
   @param[in] tupleidx Index of tuple to copy
   @param[in] seqidx Index of sequence to copy
   @param[out] tuplestr  String of tuple at tupleidx, seqidx
 */
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
