#ifndef MAF_H
#define MAF_H

#include "stdio.h"
#include "msa.h"
#include "hashtable.h"
#include "gff.h"

MSA *maf_read(FILE *F, FILE *REFSEQF, char *alphabet, int tuple_size, 
              GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order, 
              int rev_compl, int gap_strip_mode);
int maf_read_block(FILE *F, MSA *mini_msa, Hashtable *name_hash,
                 int *start_idx, int *length);
MSA *maf_read_next_msa(FILE *F);
void maf_peek(FILE *F, char ***names, Hashtable *name_hash, 
              int *nseqs, msa_coord_map *map, List *redundant_blocks,
              int *refseqlen);
void maf_block_sub_gff(GFF_Set *sub_gff, GFF_Set *gff, int start_idx, 
                       int end_idx, int *gff_idx, CategoryMap *cm,
                       int reverse_compl, int tuple_size);

#endif
