/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: category_map.h,v 1.8 2008-11-12 02:07:59 acs Exp $ */

/** \file category_map.h
    Mapping between feature types and site categories.     

    Category maps allow arbitrary 'features' (e.g., CDSs, UTRs,
   introns, ancestral repeats, CpG islands), as specified in an
   annotations file (e.g., GFF or BED), to be mapped to numbered site
   "categories" which can be treated separately in various analyses.
   They also allow a mapping backwards from site categories to
   features.  Given a category map and a set of annotations, the sites
   of an alignment can unambiguously be labeled by category, and
   conversely, given a category map and an alignment whose sites are
   labeled by category, a set of annotations can be produced
   unambiguously.  

   Category maps work as follows.  Each feature is defined by a 'type'
   (a string), and each category by an integer.  Category numbers must
   run consecutively from 1 to the total number of categories (but
   they need not be defined in order).  Categories fall into 'category
   ranges', each of which spans a sequence of consecutive integers.  A
   'simple' feature type has a category range of size one, and a
   'cyclic' feature type has a range of size greater than one.
   Multiple types may map to the same category number, but if they are
   cyclic, the ranges must coincide exactly.  The mapping from feature
   types to category numbers is unambiguous.  To allow mapping in the
   other direction, however, one type associated with a given category
   is designated as 'primary'.  By convention, the first type to be
   defined is considered the primary type.  Category 0 is reserved as the
   'background' or 'default' category (all other categories stand in
   relief to this category).

   The system is best understood by considering an example of a
   category file:

   <pre>
   NCATS = 6

   CDS         1-3
   intron      4
   5'UTR       5
   3'UTR       6
   start_codon 1-3
   stop_codon  1-3
   5'splice    4       
   3'splice    4
   </pre>

   The NCATS line must appear first; it defines the number of
   categories. Blank lines are ignored throughout.  The remaining
   lines define category ranges.  Each one must have two or three
   columns separated by whitespace.  The first column gives the
   feature type and the second the range of category numbers.  Strings
   of the form 'a-b' indicate cyclic feature types, and single
   integers indicate simple feature types.  Here CDS, start_codon, and
   stop_codon are cyclic feature types, mapped to the same range, 1-3
   (CDS is the primary feature type).  The others are all simple
   types.

   This example illustrates how a category map can be used to
   "project" a relatively large set of features onto a small set of
   categories.  It's sometimes useful to generate an extensive set of
   annotations, then to have multiple category maps that "collapse"
   categories in various ways.  For example, suppose we wanted to
   collapse the features above even further, to consider just coding
   and noncoding sites.  We could use:

   <pre>
   NCATS = 1

   CDS         1
   start_codon 1
   stop_codon  1
   </pre>

   Here we take advantage of the fact that all undefined feature types
   will map to the 'background' category (category 0).

   \ingroup feature
*/

#ifndef CAT_MAP_H
#define CAT_MAP_H

#include "hmm.h"
#include "gff.h"
#include "lists.h"
#include "stringsplus.h"

/* Node in unspooler */
typedef struct {
  List *children;
  int oldstate;
  int newstate;
} UnspoolNode;

/** Allows "unspooling" of CategoryMap, i.e., a mapping from the
    original categories to a larger set that allows for categories to
    be "conditioned on" other categories.  For example, an intron
    category might be conditioned on a CDS category, meaning that a
    different version of the intron category should be instantiated
    for each possible category in the CDS range that the intron
    follows. */
typedef struct {
  int nstates_spooled;          /**< number of original, "spooled" categories */
  int nstates_unspooled;        /**< number of "unspooled" categories (always larger) */
  int *unspooled_to_spooled;    /**< mapping from unspooled category
                                   numbers to spooled category
                                   numbers  */
  UnspoolNode **spooled_to_unspooled;
                                /**< data structure allowing for
                                   mapping from spooled category
                                   numbers to unspooled category
                                   numbers (see
                                   cm_spooled_to_unspooled) */
} Unspooler;


/** Represents range of categories for a given type */
typedef struct {
  List *feature_types;          /**< list of feature types mapped to
                                   this category range; elements of
                                   the List are to be of type
                                   String* */
  int start_cat_no;             /**< first category number in the range */
  int end_cat_no;               /**< last category noumber in the range */
} CategoryRange;

/** Category map object */
typedef struct {
  int ncats;                    /**< total number of categories */
  CategoryRange **ranges;       /**< array of CategoryRange objects.
                                   The ith element is the range for
                                   the ith category.  In the case of a
                                   range category, there will be
                                   multiple pointers to the same
                                   CategoryRange object  */
  int *labelling_precedence;    /**< array defining order of
                                   precedence to use when the same
                                   position is labeled by multiple
                                   feature types.  The ith element
                                   gives the precedence level for the
                                   ith category.  0 means highest
                                   precedence; larger values indicate
                                   lower precedence.  If unassigned,
                                   default is -1. */
  int *fill_precedence;         /**< array defining order of precedence
                                   when filling unlabelled gaps
                                   between features (see function
                                   msa_label_categories in msa.c).
                                   Conventions are the same as for
                                   labelling_precedence.  */
  List **conditioned_on;
  Unspooler *unspooler;
} CategoryMap;

/** Name of background category */
#define BACKGD_CAT_NAME "background"

CategoryMap *cm_new(int new_ncats); 

CategoryMap *cm_create_copy(CategoryMap *src);

CategoryMap* cm_create_trivial(int ncats, char *feature_prefix);

CategoryMap* cm_new_from_features(GFF_Set *feats);

CategoryMap* cm_new_string_or_file(const char *optarg);

void cm_realloc(CategoryMap *cm, int ncats);

void cm_free(CategoryMap *cm);

CategoryMap *cm_read(FILE *F);

int cm_get_category(CategoryMap *cm, String *type);

List *cm_get_category_list(CategoryMap *cm, List *names, int ignore_missing);

List *cm_get_category_str_list(CategoryMap *cm, List *names, int ignore_missing);

String *cm_get_feature(CategoryMap *cm, int cat);

String *cm_get_feature_unique(CategoryMap *cm, int cat);

List *cm_get_features(CategoryMap *cm, List *catnos);

void cm_print(CategoryMap *cm, FILE *F);

void cm_add_feature_type(CategoryMap *cm, String *type, 
                         int cycle_size); /* Size of cycle; use
                                             cycle_size=1 for a simple
                                             feature type  */

GFF_Set *cm_labeling_as_gff(CategoryMap *cm, int *path, int length, int *path_to_cat,
                            int *reverse_compl, char *seqname, char *source, 
                            List *frame_cats, char *grouproot, char *idpref);

CategoryRange* cm_new_category_range(String *type, int start_cat_no,
                                     int end_cat_no);

CategoryRange* cm_category_range_create_copy(CategoryRange *src);

void cm_free_category_range(CategoryRange *cr);

void cm_spooled_to_unspooled(CategoryMap *cm, int *path, int pathlen);

void cm_unspooled_to_spooled(CategoryMap *cm, int *path, int pathlen);

int cm_unspooled_to_spooled_cat(CategoryMap *cm, int spooled_cat);

Unspooler *cm_create_unspooler(int nstates_spooled, List **conditioned_on);

UnspoolNode *cm_new_unspool_node(int oldstate);

void cm_free_unspool_node(UnspoolNode *n);

void cm_free_unspooler(Unspooler *unsp);

int cm_get_unspooled_state(CategoryMap *cm, int spooled_state,
                           List *predecessors);

CategoryMap* cm_read_from_fname(char *fname);

List *cm_get_unspooled_list(CategoryMap *cm, List *spooled);

#endif
