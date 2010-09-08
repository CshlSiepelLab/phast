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
    files, as produced by MULTIZ and TBA.
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


typedef struct {
  String *text;  //text[i] contains line i of the block
  int numline, seqlen;

  int *specmap;  /* specmap[i]=j if line[i] of block contains information
                    about species j.  If line i isn't species-specific 
		    (ie, the first line), specmap[i] = -1 */
  char **seq;  /* seq[i] points to sequence for species i (pointer to
		  position in text array) */
  char **quality;  /* quality[i] points to string of quality scores for
                      species[i], or NULL if no data */
  char **qline, **sline, **iline, **eline;  /* pointers to lines
                                               starting with q, s, i, e for
					       each species (or NULL if
					       no data for species i) */
} MAF_BLOCK;
                 


MSA *maf_read_cats_subset(FILE *F, FILE *REFSEQF, int tuple_size, char *alphabet,
		   GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order, 
		   char *reverse_groups, int gap_strip_mode, int keep_overlapping,
			  List* cats_to_do, List *seqnames, int seq_keep);

MSA *maf_read_cats(FILE *F, FILE *REFSEQF, int tuple_size, char *alphabet,
		   GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order, 
		   char *reverse_groups, int gap_strip_mode, int keep_overlapping,
		   List* cats_to_do);
MSA *maf_read(FILE *F, FILE *REFSEQF, int tuple_size, char *alphabet,
	      GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order, 
	      char *reverse_groups, int gap_strip_mode, int keep_overlapping);

MSA *maf_read_unsorted(FILE *f, FILE *REFSEQF, int tuple_size, char *alphabet,
		  GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order,
		  char *reverse_groups, int gap_strip_mode, int keep_overlapping,
		  List *cats_to_do);

MSA *maf_read_old(FILE *F, FILE *REFSEQF, int tuple_size, char *alphabet,
              GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order, 
              char *reverse_groups, int gap_strip_mode, int keep_overlapping);

int maf_read_block(FILE *F, MSA *mini_msa, Hashtable *name_hash,
                   int *start_idx, int *length, int do_toupper);

int maf_read_block_addseq(FILE *f, MSA *mini_msa, Hashtable *name_hash,
			  int *start_idx, int *length, int do_toupper,
			  int skip_new_species);

void maf_quick_peek(FILE *f, char ***names, Hashtable *name_hash,
		    int *nseqs, int *refseqlen, int add_seqs);

void maf_peek(FILE *F, char ***names, Hashtable *name_hash, 
              int *nseqs, msa_coord_map *map, List *redundant_blocks,
              int keep_overlapping, int *refseqlen);

void maf_block_sub_gff(GFF_Set *sub_gff, GFF_Set *gff, int start_idx, 
                       int end_idx, int *gff_idx, CategoryMap *cm,
                       int reverse_compl, int tuple_size);

#endif
