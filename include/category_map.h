/* category_map - Data structures and supporting functions to map between feature types and site categories. */

/* $Id: category_map.h,v 1.2 2004-06-09 17:10:29 acs Exp $
   Written by Adam Siepel, Summer 2002
   Copyright 2002, Adam Siepel, University of California 

   Each feature is defined by a 'type' (a string), and each category
   by an integer.  Category numbers must run consecutively from 1 to
   the total number of categories (but they need not be defined in
   order).  Categories fall into 'category ranges', each of which
   spans a sequence of consecutive integers.  A 'simple' feature type
   has a category range of size one, and a 'cyclic' feature type has a
   range of size greater than one.  Multiple types may map to the same
   category number, but if they are cyclic, the ranges must coincide
   exactly.  The mapping from feature types to category numbers is
   unambiguous.  To allow mapping in the other direction, however, one
   type associated with a given category must be designated as 'primary'.
   By convention, the first type to be defined is considered the
   primary type.  Two types of precedence may be defined among
   category numbers: 'labelling precedence', which determines which
   category should be used when a given site is described by multiple
   feature types, and 'fill precedence', which describes which
   category should expand to fill an unlabeled gap, when such behavior
   is desired.  Category 0 is reserved as the 'background' or
   'default' category (all other categories stand in relief to this
   category).

   The system is best understood by considering an example of a
   category file:

   NCATS = 6

   cds         1-3
   intron      4       1,2,3
   5'UTR       5
   3'UTR       6
   start       1-3
   stop        1-3
   5'splice    4       
   3'splice    4

   LABELLING_PRECEDENCE = 1,2,3,4
   FILL_PRECEDENCE = 4,5,6
   
   The NCATS line must appear first; it defines the number of
   categories.  The LABELLING_PRECEDENCE and FILL_PRECEDENCE lines are
   optional and may appear anywhere.  The categories listed must
   appear in order of priority.  Categories not listed are assumed to
   have low priority.  Blank lines are ignored throughout.  The
   remaining lines define category ranges.  Each one must have two or
   three columns separated by whitespace.  The first column gives the
   feature type and the second the range of category numbers.  Strings
   of the form 'a-b' indicate cyclic feature types, and single
   integers indicate simple feature types.  The third column is a list
   of categories to 'condition on'.  These are used in the definition
   of a higher-order HMM (see hmm.h).

   Here cds, start, and stop are cyclic feature types, mapped to the
   same range, 1-3 (cds is the primary feature type).  The others are
   all simple types.  The 'intron' category is conditioned on the
   'cds' categories, meaning that an unspooled HMM will have different
   states for intron states following state 1, state 2, and state 3.
   The rest is fairly self-explanatory.
*/

#ifndef CAT_MAP_H
#define CAT_MAP_H

#include "hmm.h"
#include "gff.h"
#include "lists.h"
#include "stringsplus.h"

/* ADD COMMENTS */
typedef struct {
  List *children;
  int oldstate;
  int newstate;
} UnspoolNode;

typedef struct {
  int nstates_spooled;
  int nstates_unspooled;
  int *unspooled_to_spooled;
  UnspoolNode **spooled_to_unspooled;
} Unspooler;


typedef struct {
  List *feature_types;          /* list of feature types mapped to
                                   this category range; elements of
                                   the List are to be of type
                                   String* */
  int start_cat_no;             /* first category no. in the range */
  int end_cat_no;               /* last category no. in the range */
} CategoryRange;

typedef struct {
  int ncats;                    /* total number of categories */
  CategoryRange **ranges;       /* array of CategoryRange objects.
                                   The ith element is the range for
                                   the ith category.  In the case of a
                                   range category, there will be
                                   multiple pointers to the same
                                   CategoryRange object  */
  int *labelling_precedence;    /* array defining order of
                                   precedence to use when the same
                                   position is labeled by multiple
                                   feature types.  The ith element
                                   gives the precedence level for the
                                   ith category.  0 means highest
                                   precedence; larger values indicate
                                   lower precedence.  If unassigned,
                                   default is -1. */
  int *fill_precedence;         /* array defining order of precedence
                                   when filling unlabelled gaps
                                   between features (see function
                                   msa_label_categories in msa.c).
                                   Conventions are the same as for
                                   labelling_precedence.  */
  List **conditioned_on;
  Unspooler *unspooler;
  List **feat_ext_lst;          /* list of feature types for each
                                   category range to be used when
                                   generating types from category
                                   numbers (see cm_labeling_as_gff);
                                   allows default feature types to be
                                   overridden  */
} CategoryMap;


#define BACKGD_CAT_NAME "background"

/* Create new category map.
   Returns newly-allocated CategoryMap object, with initial capacity
   for ncats categories. */
CategoryMap *cm_new(int new_ncats); 

CategoryMap *cm_create_copy(CategoryMap *src);

CategoryMap* cm_create_trivial(int ncats, char *feature_prefix);

/* Reallocate a category map to allow for the specified size. */
void cm_realloc(CategoryMap *cm, int ncats);

/* Free memory associated with category map. */
void cm_free(CategoryMap *cm);

/* Read a category map from a file.
   Returns newly allocated CategoryMap on success, NULL on failure. */
CategoryMap *cm_read(FILE *F);

/* Retrieve the category number for the specified feature type.

   Returns category number or 0 (background) if no match.  In the case
   of a cyclic category, the first number in the range is returned
   (start_cat_no) */
int cm_get_category(CategoryMap *cm, String *type);

/* Map a list of strings representing either category names or
   category numbers to a list of category numbers.  Designed for use
   with programs that accept lists of categories as input.  Does not
   allow "null" categories (that is, non-background categories that
   map to background) -- when one is encountered, aborts unless
   ignore_missing == 1.

   Returns newly allocated list of integers.
*/
List *cm_get_category_list(CategoryMap *cm, List *names, int ignore_missing);

/* Retrieve the (primary) feature type associated with the specified
   category number. 
   Returns a pointer to a String representing the type.

   Warning: Do not alter the returned String object.  */
String *cm_get_feature(CategoryMap *cm, int cat);

/* Obtain a unique name based on the (primary) feature type associated
   with the specified category number.
   Returns a pointer to a newly allocated String.*/
String *cm_get_feature_unique(CategoryMap *cm, int cat);

/* Print the category map to the specified stream, in the standard format. */
void cm_print(CategoryMap *cm, FILE *F);

/* Add new feature to category map. 
   Assumes a feature does not already exist of the specified type.
   Currently does not allow labelling or fill precedence to be defined
   for the new feature.

   Warning: this function performs several reallocs and is likely to
   be somewhat inefficient (designed only for occasional use).
*/
void cm_add_feature_type(CategoryMap *cm, String *type, 
                         int cycle_size); /* Size of cycle; use
                                             cycle_size=1 for a simple
                                             feature type  */

GFF_Set *cm_labeling_as_gff(CategoryMap *cm, int *path, int length, int *path_to_cat,
                            int *reverse_compl, char *seqname, char *source, 
                            char defaultstrand, List *frame_cats, char *grouproot);

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
