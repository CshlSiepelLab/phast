/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* Code relating to the patterns of gaps and nongaps that occur at
   each site in an alignment, which are the basis of the phylo-hmm's
   indel model */

#ifndef GAPPAT_H
#define GAPPAT_H

#include <trees.h>
#include <msa.h>

typedef struct {
  int ncats;
  int ngap_cats;
  int ngap_patterns;
  int nbranches;
  TreeNode *topology;
  char **pattern;
  int *gapcat_to_cat;
  int *gapcat_to_pattern;
  int **cat_x_pattern_to_gapcat;
  int *pattern_to_node;         /* mapping from gap patterns to nodes
                                   in tree (branches above) */
  int *node_to_branch;          /* mapping from nodes in tree
                                   (branches above) to branch indices
                                   based on in-order traversal */
} GapPatternMap;

/* types of gap patterns */
typedef enum {NULL_PATTERN, DELETION_PATTERN, 
              INSERTION_PATTERN, COMPLEX_PATTERN} pattern_type;

#define GP_BASE 'b'

GapPatternMap *gp_create_gapcats(CategoryMap *cm, List *indel_cats, 
                                 TreeNode *topology, int rooted);
void gp_free_map(GapPatternMap *gpm);
void gp_set_phylo_patterns(GapPatternMap *gpm, int *patterns, MSA *msa);
pattern_type gp_pattern_type(GapPatternMap *gpm, int pattern);
void gp_tuple_matches_pattern(GapPatternMap *gpm, MSA *msa, int pattern, 
                              int *matches);

#endif
