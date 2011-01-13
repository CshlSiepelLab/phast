/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file local_alignment.h
   Functions for dealing with pairwise local alignments and coordinate transformations.
   Primarily used for pairwise local alignments as
   produced by BLASTZ.  This code is somewhat experimental.  

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

   @ingroup msa
*/

#ifndef LOC_ALN
#define LOC_ALN

#include "lists.h"
#include "stringsplus.h"
#include "msa.h"
#include "hashtable.h"

/** Keeps track of locations for Reference sequence and sequence to align */
typedef struct {
  int query_beg, /**< Beginning of reference sequence  */
   query_end, /**< Ending of reference sequence */
   target_beg, /**< Beginning of Sequence to align */
   target_end; /**< Ending of Sequence to align */
  /** @note Eventually might keep track of percent id, number of mismatches,
     number of matches, as allowed (variously) in PSL, LAV formats */
 } GaplessAlignment;

/** Information to align a sequence */
typedef struct {
  int query_beg, /**< Beginning of reference sequence  */
   query_end, /**< Ending of reference sequence */
   target_beg, /**< Beginning of sequence to align */
   target_end; /**< Ending of sequence to align */
  double score; /**< Alignment score (Optional) */
  String *target_seq; /**< Target sequence to align (Sequence Data) */
  List *gapless_alns; /**< List of type GaplessAlignment  */
} AlignmentBlock;

/** Information for a local Pairwise Alignment */
typedef struct{
  String *query_name,  /**< Name of reference sequence */
  *target_name;       /**< Name of sequence to align */
  String *query_seq,  /**< Sequence data of reference sequence */
  *target_seq;       /**< Sequence data of sequence to align */
  int query_len,     /**< Total length of reference sequence */
  target_len;    	     /**< Total length of sequence to be aligned */
  List *alignment_blocks;  /**< List of type AlignmentBlock */
} LocalPwAlignment;

/** Defines direction in which to adjust coordinates when
    transforming via a pairwise alignment.

    @see la_get_target_coord
*/
typedef enum {
	ADJUSTLEFT, /**< Adjust coords left when transforming via PwAlign */
	ADJUSTRIGHT /**< Adjust coords right when transforming via PwAlign */
        } adjust_dir;

/** Create a new local pairwise alignment object.
    @result new pairwise alignment object
 */
LocalPwAlignment *la_new();

/** Create a new local pairwise alignment object from a *.lav file
   @param F LAV file containing alignment information
   @param read_seqs Whether to read in sequences referred to in LAV file
   @result new pairwise alignment object with data from file
*/
LocalPwAlignment *la_read_lav(FILE *F, int read_seqs);

/** Create a new gap-less alignment object
   @param query_beg Beginning of Reference sequence
   @param query_end Ending of Reference sequence 
   @param target_beg Beginning of Sequence to align
   @param target_end Ending of Sequence to align 
   @result new gap-less alignment with specified sequence locations
 */
GaplessAlignment *la_new_gapless_aln(int query_beg, int query_end, 
                                     int target_beg, int target_end);

/** Create a new Alignment block object 
   @param query_beg Beginning of Reference sequence
   @param query_end Ending of Reference sequence 
   @param target_beg Beginning of Sequence to align
   @param target_end Ending of Sequence to align 
   @param score Alignment Score
   @param target_seq Sequence data for sequence to align
   @param force_global Causes all of target seq to be represented also.
   @result new alignment block with specified parameters
*/
AlignmentBlock *la_new_alignment_block(int query_beg, int query_end, 
                                       int target_beg, int target_end, 
                                       double score, String *target_seq);

/** Cast local pairwise alignment to multiple sequence alignment object.
    @param lpwa Local pairwise alignment to cast
    @param force_global Causes all of target seq to be represented also.
    @result MSA object containing data from local pairwise alignment
*/
MSA* la_to_msa(LocalPwAlignment *lpwa, int force_global);

/**
   Estimate the coordinate in the target sequence corresponding 
   to the specified coordinate in the query sequence.  
   
   @pre Alignment blocks, and gap-less alignments within them, are
   in sorted order (wrt query sequence).  
   @param lpwa Local pairwise alignment object 
   @param query_coord Specifies the coordinate in the query sequence to map to target sequence coordinates
   @param adjust Adjusts the mapping of coords i.e. ADJUSTRIGHT
   @result Estimate of coordinate in target sequence corresponding
    to coordinate in query sequence OR -1 is returned if no reasonable
    estimate is possible. 
   @note Currently does linear search (should use binary search).  
   @warning This routine should only be used with relatively close,
   generally orthologous sequences, having good synteny.
*/
int la_get_target_coord(LocalPwAlignment *lpwa, int query_coord, 
                        adjust_dir adjust);

/**  Transform the coordinates of all features in a GFF according to a
   local alignment.  

   Each feature in the original GFF will be replaced
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
   between gap-less alignments of the same block are ignored, so that a
   transformed feature may contain some bases that do not directly
   align to the query sequence.  The rationale is that these
   insertions should generally be small, and should reflect
   small-scale events that do not radically disrupt the local
   properties of the sequence.
   @param lpwa Local pairwise alignment to have coords set
   @param gff Feature Set used to create coords in lpwa
 */
void la_gff_transform(LocalPwAlignment *lpwa, GFF_Set *gff);

#endif
