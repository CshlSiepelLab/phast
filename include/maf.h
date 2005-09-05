/* $Id: maf.h,v 1.4 2005-09-05 23:03:54 acs Exp $
   Written by Adam Siepel, 2003
   Copyright 2003, Adam Siepel, University of California */

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

MSA *maf_read(FILE *F, FILE *REFSEQF, int tuple_size, char *alphabet,
              GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order, 
              char *reverse_groups, int gap_strip_mode, int keep_overlapping);

int maf_read_block(FILE *F, MSA *mini_msa, Hashtable *name_hash,
                 int *start_idx, int *length);

void maf_peek(FILE *F, char ***names, Hashtable *name_hash, 
              int *nseqs, msa_coord_map *map, List *redundant_blocks,
              int keep_overlapping, int *refseqlen);

void maf_block_sub_gff(GFF_Set *sub_gff, GFF_Set *gff, int start_idx, 
                       int end_idx, int *gff_idx, CategoryMap *cm,
                       int reverse_compl, int tuple_size);

#endif
