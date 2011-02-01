/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: maf.h,v 1.7 2009-01-09 22:01:00 mt269 Exp $ */

/** \file maf.h 
    Reading of alignments from MAF ("Multiple Alignment Format")
    files to generate Sufficient Statistics.
    (See http://www.bx.psu.edu/miller_lab.)  These functions are primarily
    concerned with extracting sufficient statistics from a MAF file
    (see sufficient_stats.c).  They avoid representing the alignment
    explicitly, and as a result allow large MAF files (e.g., spanning
    whole mammalian chromosomes) to be read and stored fairly
    efficiently.  A "reference sequence" alignment is currently
    assumed, with the reference sequence appearing first in each
    alignment block, and always in the positive (rather than reverse
    complemented) orientation (this is the convention with MULTIZ).
    \ingroup msa
*/

#ifndef MAF_H
#define MAF_H

#include "stdio.h"
#include "msa.h"
#include "hashtable.h"
#include "gff.h"

/** Hold data for a single block within a MAF file */
typedef struct {
  String *text;  /**< Maf Block contents, text[i] contains line i of the block */
  int numline, /**< Number of lines in Maf Block */
  seqlen;      /**< Sequence length of sequence in Maf Block */

  int *specmap;  /**< Maps Maf Block contents (by line) to species.
		      specmap[i]=j if line[i] of block contains information
                    about species j.  If line i isn't species-specific 
		    (ie, the first line), specmap[i] = -1 */
  char **seq;  /**< Sequence data for all species within Maf Block.
		seq[i] points to sequence for species i (pointer to
		  position in text array) */
  char **quality;  /**< Quality scores for all species within Maf Block.
			quality[i] points to string of quality scores for
                      species[i], or NULL if no data */
  char **qline,  /**< Pointer to lines starting with 'q' for each species  */
	**sline, /**< Pointer to lines starting with 's' for each species  */
	**iline, /**< Pointer to lines starting with 'i' for each species  */
	**eline; /**< Pointer to lines starting with 'e' for each species  */


} MAF_BLOCK;
                 
/** @name MAF File reading 
   \{ */
/** Read subset of an Alignment from a MAF file; subset selected by sequence and feature names. 
   @pre The MAF file must be sorted with respect to the reference sequence.  
   @param[in] F MAF file
   @param[in] REFSEQ (Optional) reference sequence.  If
                                        non-NULL, the indicated file will
                                        be used to define bases in regions
                                        of no alignment (not represented in
                                        the MAF).  File format is expected
                                        to be FASTA.  Ignored if
                                        store_order == FALSE.  If NULL and
                                        store_order == TRUE, then bases in
                                        reference seq not present in MAF
                                        are represented as Ns 
   @param[in] tuple_size tuple size for sufficient statistics 
   @param[in] alphabet (Optional) alphabet for alignment; if NULL, DEFAULT_ALPABET is assumed 
   @param[in] gff (Optional) GFF_Set.  If non-NULL,
                                        category-specific counts will be
                                        collected (cm must be non-NULL
                                        also).  The gff is assumed to use
                                        the indexing system of the
                                        reference sequence (sequence 1).
                                        Currently, a non-NULL gff implies
                                        gap_strip_mode == 1 (projection
                                        onto reference sequence).
   @param[in] cm Used for category-specific counts, ignored otherwise
   @param[in] cycle_size Label site categories
    					12...<cycle_size>...12...<cycle_size>
                                        instead of using gff and cm.
                                        Useful when stats are to be
                                        collected for non-overlapping
                                        tuples.  Use -1 to ignore.
   @param[in] store_order Whether to store order in which
                                        tuples occur.  Storing order
                                        requires more memory and larger
                                        files, and is not necessary in many
                                        cases.
   @param[in] reverse_groups  Tag defining groups in gff;
                      indicates groups on negative strand
                      should be reverse complemented.
                      Ignored if NULL.  Useful when
                      collecting counts for
                      strand-specific categories.  Can't
                      be used if store_order == TRUE
   @param[in] gap_strip_mode Gap stripping mode.  Currently,
                      if store_order == 1, may only have
                      value NO_STRIP or 1 (indicating
                      projection onto sequence 1, the
                      reference sequence).  This is
                      simply to avoid some complexity in
                      coordinate mapping. 
   @param[in] keep_overlapping If TRUE, keep overlapping blocks,
                      otherwise keep only first instance.
                      Must be FALSE if store_order ==
                      TRUE or gff != NULL or 
                      cycle_size != -1 
   @param[in] cats_to_do If non-NULL, only loads elements of the
                      alignment with features in this list.  
                      Requires gff != NULL and cm != NULL
   @param[in] seqnames If non-NULL, this is the list of sequence names
                             to keep (if seq_keep==TRUE) or discard 
                             (if seq_keep==FALSE)
   @param[in] seq_keep Only used if seqnames != NULL, determines whether
                             to keep only or discard seqnames
   @note NOTE: for now, if a GFF is defined, then all blocks are projected
     onto the reference seq (gff != NULL -> gap_strip_mode == 1).  This
     simplifies things somewhat, and it's rare that you want
     category-specific counts without projecting (because gaps in the
     reference sequence make it difficult to assign sites to categories
     rationally). 
   @note The alignment won't be constructed explicitly; instead, a sufficient-statistics representation will be extracted directly from the MAF.  
   @warning Any blocks falling out of order, or which are redundant with previous blocks, will be discarded.
 */
MSA *maf_read_cats_subset(FILE *F, FILE *REFSEQF, int tuple_size, char *alphabet,
		   GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order, 
		   char *reverse_groups, int gap_strip_mode, int keep_overlapping,
			  List* cats_to_do, List *seqnames, int seq_keep);

/** Read a subset of an Alignment from a MAF file subset selected by feature names.
   @pre The MAF file must be sorted with respect to the reference sequence.  
   @param[in] F MAF file
   @param[in] REFSEQ (Optional) reference sequence.  If
                                        non-NULL, the indicated file will
                                        be used to define bases in regions
                                        of no alignment (not represented in
                                        the MAF).  File format is expected
                                        to be FASTA.  Ignored if
                                        store_order == FALSE.  If NULL and
                                        store_order == TRUE, then bases in
                                        reference seq not present in MAF
                                        are represented as Ns 
   @param[in] tuple_size tuple size for sufficient statistics 
   @param[in] alphabet (Optional) alphabet for alignment; if NULL, DEFAULT_ALPABET is assumed 
   @param[in] gff (Optional) GFF_Set.  If non-NULL,
                                        category-specific counts will be
                                        collected (cm must be non-NULL
                                        also).  The gff is assumed to use
                                        the indexing system of the
                                        reference sequence (sequence 1).
                                        Currently, a non-NULL gff implies
                                        gap_strip_mode == 1 (projection
                                        onto reference sequence).
   @param[in] cm Used for category-specific counts, ignored otherwise
   @param[in] cycle_size Label site categories
    					12...<cycle_size>...12...<cycle_size>
                                        instead of using gff and cm.
                                        Useful when stats are to be
                                        collected for non-overlapping
                                        tuples.  Use -1 to ignore.
   @param[in] store_order Whether to store order in which
                                        tuples occur.  Storing order
                                        requires more memory and larger
                                        files, and is not necessary in many
                                        cases.
   @param[in] reverse_groups  Tag defining groups in gff;
                      indicates groups on negative strand
                      should be reverse complemented.
                      Ignored if NULL.  Useful when
                      collecting counts for
                      strand-specific categories.  Can't
                      be used if store_order == TRUE
   @param[in] gap_strip_mode Gap stripping mode.  Currently,
                      if store_order == 1, may only have
                      value NO_STRIP or 1 (indicating
                      projection onto sequence 1, the
                      reference sequence).  This is
                      simply to avoid some complexity in
                      coordinate mapping. 
   @param[in] keep_overlapping If TRUE, keep overlapping blocks,
                      otherwise keep only first instance.
                      Must be FALSE if store_order ==
                      TRUE or gff != NULL or 
                      cycle_size != -1 
   @param[in] cats_to_do If non-NULL, only loads elements of the
                      alignment with features in this list.  
                      Requires gff != NULL and cm != NULL
   @note NOTE: for now, if a GFF is defined, then all blocks are projected
     onto the reference seq (gff != NULL -> gap_strip_mode == 1).  This
     simplifies things somewhat, and it's rare that you want
     category-specific counts without projecting (because gaps in the
     reference sequence make it difficult to assign sites to categories
     rationally). 
   @note The alignment won't be constructed explicitly; instead, a sufficient-statistics representation will be extracted directly from the MAF.  
   @warning Any blocks falling out of order, or which are redundant with previous blocks, will be discarded.
 */
MSA *maf_read_cats(FILE *F, FILE *REFSEQF, int tuple_size, char *alphabet,
		   GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order, 
		   char *reverse_groups, int gap_strip_mode, int keep_overlapping,
		   List* cats_to_do);

/** Read An Alignment from a MAF file. 
   @pre The MAF file must be sorted with respect to the reference sequence.  
   @param[in] F MAF file
   @param[in] REFSEQ (Optional) reference sequence.  If
                                        non-NULL, the indicated file will
                                        be used to define bases in regions
                                        of no alignment (not represented in
                                        the MAF).  File format is expected
                                        to be FASTA.  Ignored if
                                        store_order == FALSE.  If NULL and
                                        store_order == TRUE, then bases in
                                        reference seq not present in MAF
                                        are represented as Ns 
   @param[in] tuple_size tuple size for sufficient statistics 
   @param[in] alphabet (Optional) alphabet for alignment; if NULL, DEFAULT_ALPABET is assumed 
   @param[in] gff (Optional) GFF_Set.  If non-NULL,
                                        category-specific counts will be
                                        collected (cm must be non-NULL
                                        also).  The gff is assumed to use
                                        the indexing system of the
                                        reference sequence (sequence 1).
                                        Currently, a non-NULL gff implies
                                        gap_strip_mode == 1 (projection
                                        onto reference sequence).
   @param[in] cm Used for category-specific counts, ignored otherwise
   @param[in] cycle_size Label site categories
    					12...<cycle_size>...12...<cycle_size>
                                        instead of using gff and cm.
                                        Useful when stats are to be
                                        collected for non-overlapping
                                        tuples.  Use -1 to ignore.
   @param[in] store_order Whether to store order in which
                                        tuples occur.  Storing order
                                        requires more memory and larger
                                        files, and is not necessary in many
                                        cases.
   @param[in] reverse_groups  Tag defining groups in gff;
                      indicates groups on negative strand
                      should be reverse complemented.
                      Ignored if NULL.  Useful when
                      collecting counts for
                      strand-specific categories.  Can't
                      be used if store_order == TRUE
   @param[in] gap_strip_mode Gap stripping mode.  Currently,
                      if store_order == 1, may only have
                      value NO_STRIP or 1 (indicating
                      projection onto sequence 1, the
                      reference sequence).  This is
                      simply to avoid some complexity in
                      coordinate mapping. 
   @param[in] keep_overlapping If TRUE, keep overlapping blocks,
                      otherwise keep only first instance.
                      Must be FALSE if store_order ==
                      TRUE or gff != NULL or 
                      cycle_size != -1 
   @note NOTE: for now, if a GFF is defined, then all blocks are projected
     onto the reference seq (gff != NULL -> gap_strip_mode == 1).  This
     simplifies things somewhat, and it's rare that you want
     category-specific counts without projecting (because gaps in the
     reference sequence make it difficult to assign sites to categories
     rationally). 
   @note The alignment won't be constructed explicitly; instead, a sufficient-statistics representation will be extracted directly from the MAF.  
   @warning Any blocks falling out of order, or which are redundant with previous blocks, will be discarded.
 */
MSA *maf_read(FILE *F, FILE *REFSEQF, int tuple_size, char *alphabet,
	      GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order, 
	      char *reverse_groups, int gap_strip_mode, int keep_overlapping);

/** Read a subset Alignment from a MAF file; subset selected by feature name; not necessarily read in order wrt reference sequence.
   @pre The MAF file must be sorted with respect to the reference sequence.  
   @param[in] F MAF file
   @param[in] REFSEQ (Optional) reference sequence.  If
                                        non-NULL, the indicated file will
                                        be used to define bases in regions
                                        of no alignment (not represented in
                                        the MAF).  File format is expected
                                        to be FASTA.  Ignored if
                                        store_order == FALSE.  If NULL and
                                        store_order == TRUE, then bases in
                                        reference seq not present in MAF
                                        are represented as Ns 
   @param[in] tuple_size tuple size for sufficient statistics 
   @param[in] alphabet (Optional) alphabet for alignment; if NULL, DEFAULT_ALPABET is assumed 
   @param[in] gff (Optional) GFF_Set.  If non-NULL,
                                        category-specific counts will be
                                        collected (cm must be non-NULL
                                        also).  The gff is assumed to use
                                        the indexing system of the
                                        reference sequence (sequence 1).
                                        Currently, a non-NULL gff implies
                                        gap_strip_mode == 1 (projection
                                        onto reference sequence).
   @param[in] cm Used for category-specific counts, ignored otherwise
   @param[in] cycle_size Label site categories
    					12...<cycle_size>...12...<cycle_size>
                                        instead of using gff and cm.
                                        Useful when stats are to be
                                        collected for non-overlapping
                                        tuples.  Use -1 to ignore.
   @param[in] store_order Whether to store order in which
                                        tuples occur.  Storing order
                                        requires more memory and larger
                                        files, and is not necessary in many
                                        cases.
   @param[in] reverse_groups  Tag defining groups in gff;
                      indicates groups on negative strand
                      should be reverse complemented.
                      Ignored if NULL.  Useful when
                      collecting counts for
                      strand-specific categories.  Can't
                      be used if store_order == TRUE
   @param[in] gap_strip_mode Gap stripping mode.  Currently,
                      if store_order == 1, may only have
                      value NO_STRIP or 1 (indicating
                      projection onto sequence 1, the
                      reference sequence).  This is
                      simply to avoid some complexity in
                      coordinate mapping. 
   @param[in] keep_overlapping If TRUE, keep overlapping blocks,
                      otherwise keep only first instance.
                      Must be FALSE if store_order ==
                      TRUE or gff != NULL or 
                      cycle_size != -1 
   @param[in] cats_to_do If non-NULL, only loads elements of the
                      alignment with features in this list.  
                      Requires gff != NULL and cm != NULL
   @note NOTE: for now, if a GFF is defined, then all blocks are projected
     onto the reference seq (gff != NULL -> gap_strip_mode == 1).  This
     simplifies things somewhat, and it's rare that you want
     category-specific counts without projecting (because gaps in the
     reference sequence make it difficult to assign sites to categories
     rationally). 
   @note The alignment won't be constructed explicitly; instead, a sufficient-statistics representation will be extracted directly from the MAF.  
   @warning Any blocks falling out of order, or which are redundant with previous blocks, will be discarded.
 */
MSA *maf_read_unsorted(FILE *f, FILE *REFSEQF, int tuple_size, char *alphabet,
		  GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order,
		  char *reverse_groups, int gap_strip_mode, int keep_overlapping,
		  List *cats_to_do);

/** Read An Alignment from a MAF file; not necessarily read in order wrt reference sequence.
   @pre The MAF file must be sorted with respect to the reference sequence.  
   @param[in] F MAF file
   @param[in] REFSEQ (Optional) reference sequence.  If
                                        non-NULL, the indicated file will
                                        be used to define bases in regions
                                        of no alignment (not represented in
                                        the MAF).  File format is expected
                                        to be FASTA.  Ignored if
                                        store_order == FALSE.  If NULL and
                                        store_order == TRUE, then bases in
                                        reference seq not present in MAF
                                        are represented as Ns 
   @param[in] tuple_size tuple size for sufficient statistics 
   @param[in] alphabet (Optional) alphabet for alignment; if NULL, DEFAULT_ALPABET is assumed 
   @param[in] gff (Optional) GFF_Set.  If non-NULL,
                                        category-specific counts will be
                                        collected (cm must be non-NULL
                                        also).  The gff is assumed to use
                                        the indexing system of the
                                        reference sequence (sequence 1).
                                        Currently, a non-NULL gff implies
                                        gap_strip_mode == 1 (projection
                                        onto reference sequence).
   @param[in] cm Used for category-specific counts, ignored otherwise
   @param[in] cycle_size Label site categories
    					12...<cycle_size>...12...<cycle_size>
                                        instead of using gff and cm.
                                        Useful when stats are to be
                                        collected for non-overlapping
                                        tuples.  Use -1 to ignore.
   @param[in] store_order Whether to store order in which
                                        tuples occur.  Storing order
                                        requires more memory and larger
                                        files, and is not necessary in many
                                        cases.
   @param[in] reverse_groups  Tag defining groups in gff;
                      indicates groups on negative strand
                      should be reverse complemented.
                      Ignored if NULL.  Useful when
                      collecting counts for
                      strand-specific categories.  Can't
                      be used if store_order == TRUE
   @param[in] gap_strip_mode Gap stripping mode.  Currently,
                      if store_order == 1, may only have
                      value NO_STRIP or 1 (indicating
                      projection onto sequence 1, the
                      reference sequence).  This is
                      simply to avoid some complexity in
                      coordinate mapping. 
   @param[in] keep_overlapping If TRUE, keep overlapping blocks,
                      otherwise keep only first instance.
                      Must be FALSE if store_order ==
                      TRUE or gff != NULL or 
                      cycle_size != -1 
   @note NOTE: for now, if a GFF is defined, then all blocks are projected
     onto the reference seq (gff != NULL -> gap_strip_mode == 1).  This
     simplifies things somewhat, and it's rare that you want
     category-specific counts without projecting (because gaps in the
     reference sequence make it difficult to assign sites to categories
     rationally). 
   @note The alignment won't be constructed explicitly; instead, a sufficient-statistics representation will be extracted directly from the MAF.  
   @warning Any blocks falling out of order, or which are redundant with previous blocks, will be discarded.
 */
MSA *maf_read_old(FILE *F, FILE *REFSEQF, int tuple_size, char *alphabet,
              GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order, 
              char *reverse_groups, int gap_strip_mode, int keep_overlapping);

/** Read a block from an MAF file.
   @param[in] F File containing MAF data to read
   @param[out] mini_msa MAF Block is stored here
   @param[out] name_hash Hash table mapping sequence names to sequence indices (prefix of name wrt '.' character)
   @param[out] start_idx stating coord of reference sequence
   @param[out] length Length of reference sequence
   @param[in] do_toupper Make sequences read in all upper case 
   @note Allocates memory for sequences if they are NULL (as with first block)
   @note Reads to next "a" line or EOF
   @result 0 if successful, EOF if no more blocks available
*/
int maf_read_block(FILE *F, MSA *mini_msa, Hashtable *name_hash,
                   int *start_idx, int *length, int do_toupper);

/** Add sequences from a MAF file to an existing MAF block
   @param[in] F File containing MAF data to read
   @param[in,out] mini_msa MAF Block is stored here
   @param[out] name_hash Hash table mapping sequence names to sequence indices (prefix of name wrt '.' character)
   @param[out] start_idx stating coord of reference sequence
   @param[out] length Length of reference sequence
   @param[in] do_toupper Make sequences read in all upper case 
   @param[in] skip_new_species If mini_msa already contains a given species, don't add it from MAF file
   @note Allocates memory for sequences if they are NULL (as with first block)
   @note Reads to next "a" line or EOF
   @result 0 if successful, EOF if no more blocks available
*/
int maf_read_block_addseq(FILE *f, MSA *mini_msa, Hashtable *name_hash,
			  int *start_idx, int *length, int do_toupper,
			  int skip_new_species);

/** Get partial list of sequence names (roots of names only) and length of refseq.
   @param[in] F File containing MAF Block contents
   @param[out] names List of sequence names
   @param[out] name_hash Hash Table mapping names to sequence index
   @param[out] nseqs Number of sequences in MAF Block
   @param[out] refseqlen Length of Reference Sequence
   @param[in] add_seqs Whether to add sequences to names and name_hash
   @warning May re-order names to put refseq first even if add_seqs == 0
 */
void maf_quick_peek(FILE *f, char ***names, Hashtable *name_hash,
		    int *nseqs, int *refseqlen, int add_seqs);

/** Get complete list of sequence names (roots of names only) and length of refseq.
   @param[in] F File containing MAF Block contents
   @param[out] names List of sequence names
   @param[out] name_hash Hash Table mapping names to sequence index (Hash Table must be preallocated)
   @param[out] nseqs Number of sequences in MAF Block
   @param[out] map (Optional) Coordinate map for the reference sequence (map must be pre-allocated)
   @param[out] redundant_blocks List containing indices of blocks that overlap previous blocks in 'redundant_blocks' list
   @param[in] keep_overlapping Whether to keep overlapping blocks. If 0 overlapping blocks discarded, otherwise retained
   @param[out] refseqlen Length of Reference Sequence
   @warning May re-order names to put refseq first even if add_seqs == 0
*/
void maf_peek(FILE *F, char ***names, Hashtable *name_hash, 
              int *nseqs, msa_coord_map *map, List *redundant_blocks,
              int keep_overlapping, int *refseqlen);
/** \} */

/** Extracts features from gff relevant to a specified interval.
   @pre sub_gff is allocated, with an empty feature list
   @pre gff is sorted
   @param[out] sub_gff Subset of features in gff object
   @param[in] gff Feature Set take subset from
   @param[in] start_idx Starting of interval for which to extract features
   @param[in] end_idx Ending of interval for which to extract features  
   @param[in] gff_idx Feature number to start extracting from in GFF
   @param[in] cm Category Map used for category ranges of features
   @param[in] reverse_compl Whether to reverse complement '-' strands
   @param[in] tuple_size Size of tuples, used in reverse complementing
   @note Shifts all coords such that start_idx is position 1 
*/
void maf_block_sub_gff(GFF_Set *sub_gff, GFF_Set *gff, int start_idx, 
                       int end_idx, int *gff_idx, CategoryMap *cm,
                       int reverse_compl, int tuple_size);

#endif
