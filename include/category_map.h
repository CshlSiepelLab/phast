/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file category_map.h
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

   @ingroup feature
*/

#ifndef CAT_MAP_H
#define CAT_MAP_H

#include "gff.h"
#include "lists.h"
#include "stringsplus.h"

/** Struct defines node in unspooler */
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


/** Represents range of categories for a given feature. */
typedef struct {
  List *feature_types;          /**< list of feature types mapped to
                                   this category range; elements of
                                   the List are to be of type
                                   String* */
  int start_cat_no;             /**< first category number in the range */
  int end_cat_no;               /**< last category number in the range */
} CategoryRange;

/** Category map object. */
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
                                   when filling unlabeled gaps
                                   between features (see function
                                   msa_label_categories in msa.c).
                                   Conventions are the same as for
                                   labelling_precedence.  */
  List **conditioned_on;
  Unspooler *unspooler;
} CategoryMap;

/** Name of background category. */
#define BACKGD_CAT_NAME "background"

/** \name Category Map allocation functions 
 \{ */

/** Create a new Category Map with an initial amount of categories.
    @param new_ncats Number of initial categories, can be expanded with cm_realloc.
    @result Newly allocated Category Map that can hold new_ncats categories  */
CategoryMap *cm_new(int new_ncats); 

/** Create a trivial Category Map, and populate it with # ncats categories.
   Category names are equal to category numbers (plus an optional prefix) 
   i.e. cm_create_trivial(3, "cat") yields a new Category Map
   with categories (cat1, cat2, cat3) with range sizes of one.  
  @param ncats Number of categories to create
  @param feature_prefix Prefix added to the beginning of each created category name
  @result Newly created Category Map with populated categories
*/
CategoryMap* cm_create_trivial(int ncats, char *feature_prefix);

/** Create a Category Map with a category for each feature type in a GFF_Set.  
    Category numbers are assigned in order of appearance of types.
    @param feats Features from GFF_Set to base new categories on
    @result newly created Category Map with populated categories
*/
CategoryMap* cm_new_from_features(GFF_Set *feats);

/** Copy an existing Category Map.
  @param Category Map to copy
  @result Copied category map */
CategoryMap *cm_create_copy(CategoryMap *src);

/** \} \name Category Map read/save from file functions 
 \{ */

/** Read a Category Map from a file.
   @param F File descriptor to read category map from
   @result Category Map created from file */
CategoryMap *cm_read(FILE *F);

/** Read a Category Map from a file.
  @param fname Full path to file containing Category Map
  @result Category Map created from file
*/
CategoryMap* cm_read_from_fname(char *fname);

/** Create a new Category Map from a string (filename or inlined category map)
    @param optarg filename or brief "inlined" category map e.g., "NCATS = 3 ; CDS 1-3"
    @note Useful for command-line arguments
    @result Newly created Category Map populated with specified data
*/
CategoryMap* cm_new_string_or_file(const char *optarg);

/** Save a Category Map to a file
  @param cm Category Map to save to file
  @param F File descriptor of file to save to */
void cm_print(CategoryMap *cm, FILE *F);


/** \} \name Category Map cleanup functions 
 \{ */

/** Free memory associated with category map.
   @param cm Category Map to free*/
void cm_free(CategoryMap *cm);

/** Reallocate a category map to allow for the specified size. 
    @param cm Existing Category Map that needs to be expanded
    @param ncats Number of categories the Category Map should now support
*/
void cm_realloc(CategoryMap *cm, int ncats);


/** \} \name Category Map retrieve data functions 
 \{ */

/** Get the first (base) category for a given feature type.
    @param cm Category Map containing feature type
    @param type Name of feature type to get first category for
    @result 'base' category for a given feature type OR 0 (background) if no match.
 */
int cm_get_category(CategoryMap *cm, String *type);

/** Map a list of category numbers corresponding to a given list of
   category names and or numbers.
   @param cm (Optional) CategoryMap object. May be NULL if categories are specified by number rather than by name.  (In that case, this function will reduce to conversion of strings to integers.)
   @param names List of categories.  May be specified by name or number (useful when accepting input from users)
   @param ignore_missing Whether to ignore unrecognized types.  If FALSE, then function will abort when it encounters an unrecognized type.
   @result List of Category Numbers corresponding to list of names (or numbers)
 */
List *cm_get_category_list(CategoryMap *cm, List *names, int ignore_missing);


/** Return a list of category names corresponding to a given list of
   category names and or numbers.  Doesn't allocate new names,
   just pointers to Strings in the Category Map object or the
   provided List
   @param cm (Optional) CategoryMap object.  May be NULL if categories are specified by name rather than by number.  (In that case, this function will create a copy of names with pointers to its elements)
   @param names  List of categories specified by name or number
   @param ignore_missing  Whether to ignore unrecognized types.  If FALSE, then function will abort when it encounters an unrecognized type.
   @result List of category names corresponding to a given list of category names and or numbers
   @note Useful when accepting input from users
   */
List *cm_get_category_str_list(CategoryMap *cm, List *names, int ignore_missing);

/** Retrieve the (primary) feature type associated with the specified
   category number.
   @param cm Category Map
   @param cat Category Number to get feature type associated with
   @param Feature type associated with category number (cat)
   @result Name of feature associated with category number
   @note return value is passed by reference -- do not alter. */
String *cm_get_feature(CategoryMap *cm, int cat);

/** Obtain a unique name based on the (primary) feature type
   associated with the specified category number.
   @param cm Category Map
   @param cat Category Number to get unique name for
   @result Unique name based on the (primary) feature type for given category number (cat)
   @note Return value is newly allocated string
   */
String *cm_get_feature_unique(CategoryMap *cm, int cat);

/** Get list of category names corresponding to list of category
   numbers.
   @param cm Category Map
   @param catnos List of category numbers to get names for
   @result List of category names
    */
List *cm_get_features(CategoryMap *cm, List *catnos);

/** \} \name Category Map modification functions 
 \{ */

/** Add new feature to Category Map.
    @pre Feature of the specified type must not already exist.
    @param cm Category Map to add feature to
    @param type Feature type specified by name
    @param cycle_size Size of cycle (use 1 for a a simple feature type)
*/
void cm_add_feature_type(CategoryMap *cm, String *type, 
                         int cycle_size);

/** \} \name Category Map Spooling/Unspooling functions 
 \{ */

/**
  Create an unspooler given a spooled state
   @param nstates_spooled Number of states to be unspooled
   @param conditioned_on How each spooled element is conditioned
   @result New unspooler
   @note Conditioned_on must be an array of integer lists; specifically, the
   ith element must be the list of state numbers on which the ith
   state is conditioned. */
Unspooler *cm_create_unspooler(int nstates_spooled, List **conditioned_on);

/** Create a node of the unspooler.
  @param oldstate state of spooled node
  @result new unspooled node
*/
UnspoolNode *cm_new_unspool_node(int oldstate);


/** Maps category numbers from the spooled space to unspooled space.
   @pre cm->unspooler must exist
   @param[in] cm Category Map
   @param[in,out] path Sequence (array) of spooled category numbers to map to unspooled space (using current unspooler)
   @param[in] pathlen Number of category numbers in path
   @warning Original sequence (in path) is overwritten
*/
void cm_spooled_to_unspooled(CategoryMap *cm, int *path, int pathlen);

/** Maps category numbers from the unspooled space to spooled space.
   @pre cm->unspooler must exist
   @param[in] cm Category Map
   @param[in,out] path Sequence (array) of unspooled category numbers to map to spooled space (using current unspooler)
   @param[in] pathlen Number of category numbers in path
   @warning Original sequence (in path) is overwritten
*/
void cm_unspooled_to_spooled(CategoryMap *cm, int *path, int pathlen);

/** Translate unspooled category number to spooled category number
   @param cm Category map
   @param spooled_cat Unspooled category number
   @result Spooled category number */
int cm_unspooled_to_spooled_cat(CategoryMap *cm, int spooled_cat);


/** 
  Gets the unspooled state of a spooled state
  @param cm Category Map
  @param spooled_state Index of spooled state to get as unspooled
  @param predecessors List of states (nodes) preceding the node called spooled_state
  @result unspooled version of state specified
  @note Last item in predecessors is assumed to be the most recently visited
*/
int cm_get_unspooled_state(CategoryMap *cm, int spooled_state,
                           List *predecessors);

/** Translate list of spooled category names/numbers, to a list of
   corresponding unspooled category numbers 
   @param cm Category Map
   @param spooled List of spooled numbers for Category Names
   @result list of unspooled category numbers*/
List *cm_get_unspooled_list(CategoryMap *cm, List *spooled);

/** Free a node of the unspooler.
  @param n node to free
*/
void cm_free_unspool_node(UnspoolNode *n);

/** Free the entire unspooler.
  @param unsp unspooler to free
*/
void cm_free_unspooler(Unspooler *unsp);

/** \} \name Category Map misc. functions 
 \{ */

/** Create a GFF_Set from a sequence of category/state numbers, using
   a specified category map and mapping from raw state numbers to
   category numbers.
   @param cm Category Map to use in mapping
   @param path Raw sequence of state/category numbers
   @param length Length of sequence
   @param path_to_cat Mapping from raw numbers to category numbers
   @param reverse_compl Array of boolean values indicating whether each raw-sequence value corresponds to the reverse strand
   @param seqname Char string to use as 'seqname' in generated GFF_Set
   @param source Char string to use as 'source' in GFF_Set
   @param frame_cats Categories for which to obtain frame information (by name)
   @param grouptag Tag to use to define groups in GFF_Set (e.g., "transcript_id")
   @param idpref (Optional) Prefix for ids of predicted elements. Can be used to ensure ids are unique
   @result Newly created Feature Set
    */
GFF_Set *cm_labeling_as_gff(CategoryMap *cm, int *path, int length, int *path_to_cat, int *reverse_compl, char *seqname, char *source, List *frame_cats, char *grouptag,char *idpref);

/** \} \name Category Range Allocation functions 
 \{ */

/** Create a range of Category numbers with a given feature type.
  @param type Feature type of Category Range
  @param start_cat_no Starting category number
  @param end_cat_no Ending category number
  @result New Category Range as specified [start_cat_no to end_cat_no]
*/
CategoryRange* cm_new_category_range(String *type, int start_cat_no,
                                     int end_cat_no);

/**  Copy an existing Category range. 
  @param src Existing Category Range to copy from
  @result Copy of initial Category Range
*/
CategoryRange* cm_category_range_create_copy(CategoryRange *src);

/** \} \name Category Range Cleanup function 
\{ */

/** Free a Category Range.
   @param cr Category Range to free*/
void cm_free_category_range(CategoryRange *cr);

/** \} */

#endif
