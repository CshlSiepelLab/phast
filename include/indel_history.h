/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: indel_history.h,v 1.4 2008-11-12 02:07:59 acs Exp $ */

#ifndef IND_HIST
#define IND_HIST

#include <stdio.h>
#include <trees.h>
#include <msa.h>

#define NINDEL_CHARS 3
typedef enum {INS, DEL, BASE} indel_char; /* note: order is used in places */

typedef struct {
  indel_char type;
  int start;
  int len;
} Indel;

typedef struct {
  TreeNode *tree;
  int ncols;
  List **indels;
} CompactIndelHistory;

typedef struct {
  TreeNode *tree;
  int ncols;
  char **indel_strings;            /* make bin vector later */
} IndelHistory;

IndelHistory *ih_new(TreeNode *tree, int ncols);
void ih_free(IndelHistory *ih);
CompactIndelHistory *ih_new_compact(TreeNode *tree, int ncols);
void ih_free_compact(CompactIndelHistory *cih);
IndelHistory *ih_expand(CompactIndelHistory *cih);
CompactIndelHistory *ih_compact(IndelHistory *ih);
void ih_print(IndelHistory *ih, FILE *outf, char *msa_name, char *prog_name);
void ih_print_compact(CompactIndelHistory *cih, FILE *outf, char *msa_name, 
                      char *prog_name);
MSA *ih_as_alignment(IndelHistory *ih, MSA *msa);
CompactIndelHistory *ih_read_compact(FILE *inf);
IndelHistory *ih_new_from_file(FILE* inf);
IndelHistory *ih_extract_from_alignment(MSA *msa, TreeNode *tree);
IndelHistory *ih_reconstruct(MSA *msa, TreeNode *tree);
void ih_convert_ia_names(MSA *msa, TreeNode *tree);

#endif
