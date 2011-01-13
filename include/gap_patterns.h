/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/**
  @file gap_patterns.h
   Code relating to the patterns of gaps and non-gaps that occur at
   each site in an alignment, which are the basis of the phylo-hmm's
   indel model 
   \ingroup phylo_hmm
*/

#ifndef GAPPAT_H
#define GAPPAT_H

#include <trees.h>
#include <msa.h>

/** Information about a pattern of either gaps or non-gaps */
typedef struct {
  int ncats;			/**< Number of categories from data file */
  int ngap_cats;		/**< Number of gap categories */
  int ngap_patterns;		/**< Number of phylogenetic gap patterns (2 per branch + 1 for complex gap patterns) */
  int nbranches;		/**< Number of branches from data file  */
  TreeNode *topology;		/**< Topology of nodes from data file  */
  char **pattern;		
  int *gapcat_to_cat;		/**< Mapping of gap categories to (spooled) categories */
  int *gapcat_to_pattern;	/**< Mapping of gap categories to gap patterns (all initially -1) */
  int **cat_x_pattern_to_gapcat; /**< Mapping from x gap patterns to gap cats */
  int *pattern_to_node;         /**< mapping from gap patterns to nodes
                                   in tree (branches above) */
  int *node_to_branch;          /**< mapping from nodes in tree
                                   (branches above) to branch indices
                                   based on in-order traversal */
} GapPatternMap;

/** Types of gap patterns */
typedef enum {NULL_PATTERN, /**< Pattern is not related to changes*/
 	      DELETION_PATTERN,  /**< Pattern is related to deletions */
              INSERTION_PATTERN, /**< Pattern is related to insertions */
 	      COMPLEX_PATTERN  /**< Pattern is related to changes more complex than indels */
} pattern_type;

#define GP_BASE 'b'

/** \name Gap Pattern allocation function
 \{ */

/** Infer gap cats 
  @param[in,out] cm Category Map containing initial categories
  @param[in] indel_cats List of categories associated with insertions/deletions
  @param[in] topology Topology of the tree associated with the categories
  @param[in] rooted If the Topology is rooted
  @result Object that defines the key mappings between the original
   (spooled) categories and the new "gapped categories", defined as
   category x gap pattern pairs.
  @note Expands the set of categories in the specified cat map such that
   each designated "indel category" is replaced by ngap_patterns
   categories, where ngap_patterns is a function of the tree topology.
   In addition, for category ranges of size greater than one,
   categories corresponding to non-empty gap patterns are "conditioned
   on" categories corresponding to gapless sites, and the catmap's
   unspooler is redefined accordingly (allows for "cyclic" gaps, as in
   coding regions).
  @warning Category Map (cm) passed in is modified
*/
GapPatternMap *gp_create_gapcats(CategoryMap *cm, List *indel_cats, 
                                 TreeNode *topology, int rooted);

/** \name Gap Pattern cleanup function
 \{ */
/** Free gap cats object */
void gp_free_map(GapPatternMap *gpm);

/** \name Gap Pattern misc. functions
 \{ */

/** Set a mapping between phylogenetic patterns and multiple sequence aligned data.
   Given a multiple alignment and a gap pattern map, fill an array
   with indices describing a "phylogenetic gap pattern" at each
   position.  The value for each column will be an integer i such that
   i == 0 if there are no gaps, 1 <= i <= 2N-3 (where N is the number
   of leaves in the tree) if gaps occur that can be explained by a
   deletion on branch i, 2N-2 <= i <= 4N-6 if gaps occur that can be
   explained by an insertion on branch i, and i == 4N-5 if gaps occur
   that require more than one indel event to explain (so-called
   "complex" gap patterns).  The branches of the tree are numbered
   from 1 to 2N-3 in an in-order traversal (see gp_create_gapcats).
   The branch between the root and its right child is ignored, because
   insertions on one branch below the root cannot be distinguished
   from deletions on the other.
  @pre Alignment is not 100% gap characters
  @param[in,out] gpm Gap Pattern Map to add to
  @param[in] patterns Patterns to map to msa at least length msa->length
  @param[in] msa Sequence data to map patterns to
*/
void gp_set_phylo_patterns(GapPatternMap *gpm, int *patterns, MSA *msa);

/** Translating an integer to pattern type */
pattern_type gp_pattern_type(GapPatternMap *gpm, int pattern);

/** Retrieve tuples in sequence data that match a pattern.
  @param[in] gpm Gap Pattern Map to map between sequence and pattern
  @param[in] msa Multiple Sequence Alignment data
  @param[in] pattern Pattern to look for
  @param[out] matches Locations where matches are found
*/
void gp_tuple_matches_pattern(GapPatternMap *gpm, MSA *msa, int pattern, 
                              int *matches);

/** \} */
#endif
