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
  int *gapcat_to_cat;
  int *gapcat_to_pattern;
  int **cat_x_pattern_to_gapcat;
} GapPatternMap;

GapPatternMap *gp_create_gapcats(CategoryMap *cm, List *indel_cats, int nseqs);
void gp_free_map(GapPatternMap *gpm);
void gp_set_phylo_patterns(int *patterns, MSA *msa, TreeNode *tree);

#endif
