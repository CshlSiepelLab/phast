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

GapPatternMap *gp_create_gapcats(CategoryMap *cm, List *indel_cats, 
                                 TreeNode *topology);
void gp_free_map(GapPatternMap *gpm);
void gp_set_phylo_patterns(GapPatternMap *gpm, int *patterns, MSA *msa);
pattern_type gp_pattern_type(GapPatternMap *gpm, int pattern);

#endif
