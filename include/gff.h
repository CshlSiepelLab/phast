/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file gff.h
    Reading and writing of sequence features in General Feature Format (GFF). 
    Obeys file specification at
    http://web.archive.org/web/20070823232227/http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml
    @ingroup feature
*/


#ifndef GFF_H
#define GFF_H

#include <stringsplus.h>

#include <stdio.h>
#include <stdlib.h>

/** GFF_Feature structure.  Simply mirrors file-format spec (see
   http://web.archive.org/web/20070823232227/http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml) */ 
typedef struct {
  String *seqname;              /**< name of sequence (a single GFF file
                                   may describe multiple sequences) */
  String *source;               /**< source of feature -- usually
                                   program or database */
  String *feature;              /**< feature type; so far appears only
                                   semi-standardized.  One suggestion
                                   is to use the EMBL/DDBJ/GenBank
                                   feature table as a standard */
  int start,					/**< Start position.  Convention
									is to start numbering with 1 */
  end;              		 	/**< End positions.
                                   Convention is make the range
                                   inclusive. */
  double score;                 /**< arbitrary floating-point score.  If
                                   null, set score_is_null (see
                                   below) */
  char strand;                  /**< one of '+', '-', and '.' */
  int frame;                    /**< reading frame.  Should be 0-2 or
                                   GFF_NULL_FRAME (represented as '.'
                                   in files).  NOTE: internal storage
                                   is not the same as the GFF standard
                                   (conversion is performed at input
                                   and output time).  Internal storage
                                   is offset from start of codon mod
                                   3, i.e., the positions of a codon
                                   are assigned frames 0, 1, 2.  The
                                   GFF standard is no. positions to
                                   start of next codon, i.e.,
                                   positions are assigned frames 0, 2,
                                   1 */
  String *attribute;            /**< string describing auxiliary data in
                                   tag-value format.  See
                                   http://web.archive.org/web/20070823232227/http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml. */
  int score_is_null;            /**< whether score is null (if null,
                                   represented as '.' in files).  This
                                   extra field necessary because all
                                   real numbers are potentially
                                   legitimate scores */
} GFF_Feature;


/** GFF_Set: a set of GFF_Feature objects, generally as appear
     together in a file.  Simply consists of a list of GFF_Feature
     objects and some optional meta-data extracted from file comments
     (see
     http://web.archive.org/web/20070823232227/http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml).
     NOTE: "type" meta-data is currently ignored */
typedef struct {
  List *features;               /**< list of GFF_Feature objects */
  String *gff_version;          /**< version of GFF in use; generally 2 */
  String *source;               /**< program used to generate file */
  String *source_version;       /**< version of program used to generate file */
  String *date;                 /**< date of generation */
  List *groups;                 /**< used when grouping features by attribute */
  String *group_tag;            /**< tag defining groups */
} GFF_Set;

/* Group of features:  used by gff_group and related functions */
typedef struct {
  String *name;
  List *features;
  int start, end;
} GFF_FeatureGroup;

/** total number of columns */
#define GFF_NCOLS 9             
/** minimum allowable number of columns */
#define GFF_MIN_NCOLS 5         
/** starting number of features */
#define GFF_SET_START_SIZE 1000	

/** tags that identify meta-data in comments */
#define GFF_VERSION_TAG "gff-version"
#define GFF_SOURCE_VERSION_TAG "source-version"
#define GFF_DATE_TAG "date"

#define GFF_DEFAULT_VERSION 2

/** use when frame is null */
#define GFF_NULL_FRAME -1       

/** default feature types */
#define GFF_CDS_TYPE "CDS"
#define GFF_EXON_TYPE "exon"
#define GFF_INTRON_TYPE "intron"
#define GFF_START_TYPE "start_codon"
#define GFF_STOP_TYPE "stop_codon"
#define GFF_UTR5_TYPE "5'UTR"
#define GFF_UTR3_TYPE "3'UTR"
#define GFF_SPLICE5_TYPE "5'splice"
#define GFF_SPLICE3_TYPE "3'splice"



/** \name Feature allocation functions
 \{ */

/** Create new GFF_Feature object with specified attributes, parameters
    are copied by reference.  
    @param seqname Sequence Name
    @param source Source of feature, usually a program or database
    @param feature Feature type
    @param start Starting column for the feature
    @param end Ending column for the feature
    @param score Floating-point score. if null, set score_is_null
    @param strand One of '+', '-', and '.' 
    @param frame Reading frame. (should be 0-2 or GFF_NULL_FRAME)
    @param attribute String describing auxiliary data in tag-value format
    @param score_is_null If score in above parameter is null or valid
    @result Returns newly allocated GFF_Feature object.
*/
GFF_Feature *gff_new_feature(String *seqname, String *source, String *feature,
                             int start, int end, double score, char strand, 
                             int frame, String *attribute, int score_is_null);


/** Create an exact copy of a GFF_Feature object.
    @param orig Feature to be copied
    @result Exact copy of original
 */
GFF_Feature *gff_new_feature_copy(GFF_Feature *orig);


/** Create new GFF_Feature object with specified attributes, parameters
    are copied by value.  
    @param seqname Sequence Name
    @param source Source of feature, usually a program or database
    @param feature Feature type
    @param start Starting column for the feature
    @param end Ending column for the feature
    @param score Floating-point score. if null, set score_is_null
    @param strand One of '+', '-', and '.' 
    @param frame Reading frame. (should be 0-2 or GFF_NULL_FRAME)
    @param attribute String describing auxiliary data in tag-value format
    @param score_is_null If score in above parameter is null or valid
    @result Returns newly allocated GFF_Feature object.

 */
GFF_Feature *gff_new_feature_copy_chars(const char *seqname, const char *source,
					const char *feature,
					int start, int end, double score, 
					char strand, int frame, 
					const char *attribute,
					int score_is_null);

/** Create a new GFF_Feature from a genomic position string of the
   type used in the UCSC genome browser. 
   @param position Location of the feature i.e. "chr10:102553847-102554897"
   @param source Source of feature, usually a program or database
   @param feature Feature type
   @param score arbitrary floating-point score. if null, set score_is_null
   @param frame Reading frame. (should be 0-2 or GFF_NULL_FRAME)
   @param attribute String describing auxiliary data in tag-value format
   @param score_is_null If score in above parameter is null or valid
   @note A trailing '+' or '-' will be interpreted as the strand; otherwise the null strand is used.  
   @result NULL is returned if the string can't be parsed, otherwise new Feature  */
GFF_Feature *gff_new_feature_genomic_pos(String *position, String *source, 
                                         String *feature, double score, 
                                         int frame, String *attribute,
                                         int score_is_null);

/** \name Feature Set allocation functions
 \{ */

/** Create new GFF_Set object with empty string attributes. 
 */
GFF_Set *gff_new_set();

/** Create new GFF_Set object with the number of rows specified.
    @param len Number of features new Feature Set can hold
    @result new Feature Set of a specific size
 */
GFF_Set *gff_new_set_len(int len);

/** Create new GFF_Set object, using typical defaults and other
   parameters as specified.  
   Sets gff version to '2' and date to current date, and sets 
   source and source version as specified. 
    
   @param source Where the feature came from
   @param source_version Which version of the source
*/
GFF_Set *gff_new_set_init(char *source, char *source_version);

/** Create new GFF_Set object, initializing with same version, source,
    etc. as a template GFF_Set object
    @param gff Original feature set
    @result Copied version of gff input
 */
GFF_Set *gff_new_from_template(GFF_Set *gff);

/** \name GFF cleanup functions
 \{ */

/** Clear all features from a GFF_Set and free associated memory.
    @param gff Feature Set to empty, freeing anything emptied
 */
void gff_clear_set(GFF_Set *gff);

/** Free resources associated with GFF_Set object (including all
   features and the set object itself). 
   @param set Feature Set to free
*/
void gff_free_set(GFF_Set *set);

/** Free resources associated with GFF_Feature object. 
    @param feat Feature to free
 */
void gff_free_feature(GFF_Feature *feat);

/** \name GFF read/save access file functions
 \{ */

/** Read a set of features from a file and return a newly allocated
   GFF_Set object.  
   Function reads until end-of-file is encountered or
   error occurs (aborts on error).  Comments and blank lines are
   ignored and special "meta-data" comments are parsed (see
   http://web.archive.org/web/20070823232227/http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml).  Only
   the first five columns of feature lines are considered required
   ('name', 'source', 'feature', 'start', and 'end'); subsequent
   fields are optional ('score', 'strand', 'frame', and 'attribute').
   Default value for score, strand, and frame is null ('.') and for
   attribute is the empty string ('').  Columns must be separated by
   tabs.  
   @param F File descriptor of file containing Feature Set
   @result Feature Set created from File descriptor
   @see gff_read_from_fname
*/
GFF_Set* gff_read_set(FILE *F);

/** Read a gff from file given a filename.
    For more specifics or to read from a file descriptor see gff_read_set.
    @param fname Full path and filename of file to create Feature Set from
    @result Feature Set created from file
    @see gff_read_set
*/
GFF_Set* gff_read_from_fname(char *fname);

/** Output a GFF_Set to the specified stream in GFF.
    @param F File to save to
    @param set Feature Set to save to F
 */
void gff_print_set(FILE *F, GFF_Set *set);

/** Print an individual GFF_Feature object as a GFF line. 
    @param F File to save to
    @param feat Feature to save to F
*/
void gff_print_feat(FILE *F, GFF_Feature *feat);

/** \} \name GFF Grouping functions
 \{ */


/** Group features by value of specified tag.  All features with
    undefined values will be placed in a single group. 
    @param set Feature Set to group by tag
    @param tag String to group by
*/
void gff_group(GFF_Set *set, char *tag);

/** Group features by feature type.  
    All features with undefined values will be placed 
    in a single group. 
    @param set Feature Set to group by feature type
*/
void gff_group_by_feature(GFF_Set *set);

/** Group features by sequence name
    @param set Feature set
 */
void gff_group_by_seqname(GFF_Set *set);


/** Group contiguous features, e.g., an exon and adjacent splice
    sites.  If features have already been grouped (e.g., by transcript
    id), then subgroups are created by adding new tags. New tag values
    will be composed of the tag value for the "outer" group (e.g., the
    transcript_id) and a unique suffix indicating the "inner" group.
    In all cases, feature are sorted as a side effect, in a way that
    reflects the initial grouping, not the grouping into exons. 
   @param set Feature Set to process
   @param tag Tag to use for new groups (e.g. "exon_id") 
*/
void gff_exon_group(GFF_Set *set, char *tag);

/** Identify overlapping groups and remove all but the first
   one encountered.  
  @pre Features must already be grouped. 
  @param gff Feature Set to process
  @param discards_f If non-NULL, discarded features will be written here
*/
void gff_remove_overlaps(GFF_Set *gff, FILE *discards);

/** Remove all groups whose names are not in the specified list.
    @param feats Feature Set to filter
    @param groups list of groups to keep 
 */
void gff_filter_by_group(GFF_Set *feats, List *groups);

/** Find the name of the group a feature belongs to.
   @param feats Feature Set which contains feature f
   @param f Feature whos group we are looking for
   @result Name of group feature f belongs to or NULL if f belongs to no group
   @warning Not efficient, must search through all groups until feature is found */
String *gff_group_name(GFF_Set *feats, GFF_Feature *f);

/** Copies a GFF set but not any groups.
    @param orig Original Feature Set to copy
    @result New Feature Set identical to original, except for groups
 */
GFF_Set *gff_copy_set_no_groups(GFF_Set *orig);

/** Remove grouping of features.
    @param set Feature Set to ungroup
 */
void gff_ungroup(GFF_Set *set);

/** \} \name GFF Flatten functions
 \{ */

/** Merges overlapping or adjacent features of same type.  
    @pre Features are sorted 
    @pre Features frames must be null (GFF_NULL_FRAME)
    @note Attributes are ignored in merge
    @note Groups are ignored in merge and removed unless no features are changed
     */
void gff_flatten(GFF_Set *feats);

/** Merges overlapping or adjacent features of same type, if they
    are they are in the same group.  Sorts the features within groups.
    When two features are merged, scores are summed, but attributes are 
    ignored.  Will not merge if 'frame' is non-null.  */
void gff_flatten_within_groups(GFF_Set *feats);

/** Flatten a GFF without regard to strand, score, feature type (though do merge these
    appropriately).  
    Does pay attention to seqname (so elements on different chromosomes are not merged).
    @param feats Feature Set to merge
    @result Number bases after merge
 */
long gff_flatten_mergeAll(GFF_Set *feats);

/** \} \name GFF Comparator functions
 \{ */

/** Compare two features to determine which should be first in a list, used for sorting features.
    @param ptr1 Pointer to first feature
    @param ptr2 Pointer to second feature
    @result How much lower ptr1 is than ptr2
 */
int gff_feature_comparator(const void* ptr1, const void* ptr2);

/** Compares two groups to see which should be first in a list, used for sorting groups. 
    @param ptr1 Pointer to first group
    @param ptr2 Pointer to second group
    @result How much lower ptr1 is than ptr2
*/
int gff_group_comparator(const void* ptr1, const void* ptr2);

/** \} \name GFF Subset functions
 \{ */

/** Create a new GFF_Set representing the features fully within a particular
    coordinate range.  
    Every site of the feature must be within the range for it to be included.
    @param set Feature Set to create subset from
    @param startcol All subset features must start at or after this column number
    @param endcol All subset features must end at or before this column number
    @param reset_indices Used to set indices of features in result relative to startcol
    @result new GFF_Set with features within startcol to endcol
 */
GFF_Set *gff_subset_range(GFF_Set *set, int startcol, int endcol, 
                          int reset_indices);

/** Create a new GFF_Set representing the features partially or fully in a
    particular  coordinate range. 
    As long as any part of the feature is within the range, it is included 
    @param set Feature Set to create subset from
    @param startcol All subset features must have one or more sites at or after this column number
    @param endcol All subset features must have one or more sites at or before this colum number
    @result new GFF_Set with features partially or fully within startcol to endcol
 */
GFF_Set *gff_subset_range_overlap(GFF_Set *set, int startcol, int endcol);

/** Create a new GFF_Set representing the features partially or fully in a
    particular coordinate range, starting at a specific index in GFF. 
  
    @pre GFF must be sorted by start position of feature.
    @param set Feature Set to create subset from
    @param startcol All subset features must have one or more sites at or after this column number
    @param endcol All subset features must have one or more sites at or before this column number
    @param startSearchIdx Column to start looking for features within GFF
    @result new GFF_Set with features partially or fully within startcol to endcol
    @note For features to be included, they must not overlap startSearchIdx
    @note startSearchIdx set to first match if any found
*/
GFF_Set *gff_subset_range_overlap_sorted(GFF_Set *set, int startcol, int endcol, 
					 int *startSearchIdx);


/** \} \name GFF Add extra features by type
 \{ */


/** Add a gene_id tag, along with whatever other tags are in use.
    @pre Feats must contain non-NULL groups
    @param feats Features Set containing features to add gene id attributes to
    @note Required in output by some programs 
*/
void gff_add_gene_id(GFF_Set *feats);

/** Adds 5'UTR and 3'UTR features for cases in which exon features
    extend beyond cds features of the same group.
    @pre Feats must contain non-NULL groups
    @param feats Features Set containing groups and features to add utr features to    
*/
void gff_create_utrs(GFF_Set *feats);

/** Adds intron features between exons of the same group.
    @pre feats must contain non-NULL groups
    @param feats Features Set containing groups and features to add intron features to.
 */
void gff_create_introns(GFF_Set *feats);

/** Adds signal features for start and stop codons and 5' and 3' splice sites.  
    @pre Feats must contain non-NULL groups
    @param feats Features Set containing groups and features to add start/stop codon features to
    @note Splice sites in UTR will only be added if UTR is annotated (see function above) */
void gff_create_signals(GFF_Set *feats);

/** \} \name GFF Misc.
 \{ */

/** Discard any feature whose feature type is not in the specified
    list.
    @param gff Feature Set to process
    @param types Feature types to include (List of strings)
    @param exclude Exclude rather than include specified types
    @param discards_f Discarded features will be written here (if non-NULL)
 */
void gff_filter_by_type(GFF_Set *gff, List *types, int exclude, 
                        FILE *discards_f);

/** Test whether a set of GFF_Feature objects refers to the reverse
   strand.  
   @param features Features List to test for reference to reverse strand
   @result 1 if no features have strand equal to '+' and at
   least one has strand equal to '-'; otherwise returns 0. */
int gff_reverse_strand_only(List *features);

/** Adjust coordinates and strand of GFF_Feature objects to reflect
   reverse complement of given interval of sequence.  Also
   reverses order of appearance of features.  
   @pre Features, start_range, and end_range must all use the same
        coordinate frame. 
   @param features List of GFF_Feature objects
   @param start_range First coordinate of internal (inclusive, 1-based indexing, as in features) 
   @param end_range Last coordinate of interval (inclusive, 1-based indexing, as in features)
*/
void gff_reverse_compl(List *features, int start_range, int end_range);

/** Sort features primarily by start position, and secondarily by end
    position (ascending order).  
    @param set Features to sort
    @note If features are grouped (see gff_group), then they will 
          be sorted within groups, and groups
          will be sorted by start position of first feature */
void gff_sort(GFF_Set *set);

/** Adjust coords of CDS features such that start codons are included
    and stop codons are excluded, as required in GTF2.  
    @pre GFF must be grouped such that at most one start codon 
         and at most one stop codon occur per group. 
    @param gff Feature Set containing CDS features
*/
void gff_fix_start_stop(GFF_Set *gff);

/** Adjust coords of features of "primary" types (e.g., CDS) to
    include any features of "helper" types (e.g., start_codon).
    @pre Features must be grouped and sorted
    @pre feats must contain non-NULL groups
    @param feats Features Set containing primary features  
    @param primary_types List of types of features in feats that should be adjusted
    @param helper_types List of types of features that are 'helpers' and should be absorbed
*/
void gff_absorb_helpers(GFF_Set *feats, List *primary_types, 
                        List *helper_types);



/** Partition features of a GFF by feature type.  
    @param[in] feats Feature Set
    @param[out] types List of feature types corresponding to subsets
    @param[out] subsets List of lists of features corresponding to types
*/
void gff_partition_by_type(GFF_Set *feats, List *types, List *subsets);


/** Adds offset to start and end of each feature.
   @param gff Features to modify
   @param offset Amount to add to start and end of each feature
   @param maxCoord All features that start or end above maxCoord are removed
   @note If new end coordinate < 1, or if maxCoord > 0 and new
   start coordinate > maxCoord, then feature is removed.
   @note If part of feature is out of bounds, start and end
   may be truncated to be in the range [1,maxCoord] (or
   [1, infinity) if maxCoord < 0 
*/
void gff_add_offset(GFF_Set *gff, int offset, int maxCoord);

/**
  Create a new GFF by overlaps two existing GFFs and specifying the amount of overlap desired.
  The amount of overlap is defined by the parameters passed in.
  @param First feature set to overlap with the second feature set
  @param Second feature set to overlap with the first feature set
  @param numbase_overlap Minimum number of overlapping bases to consider fragments overlapping (or -1 to use percentOverlap only)
  @param percentOverlap Minimum percent of bases in a single feature of filter_gff that must overlap with gff to consider fragments overlapping (or use -1 to use numbase_overlap only).
  @param nonOverlapping If TRUE, return non-overlapping fragments instead of overlapping.
  @param overlappingFragments If FALSE, entire records of gff are selected based on whether they meet overlap criteria with filterGff.  Otherwise, only overlapping portion of gff fragments are returned.  In this case, the same fragments may be returned multiple times if they overlap with multiple elemnts of gffFilter.
  @param[output] overlapping_frag (Only used when overlappingFragments is TRUE).  If non-NULL, should be an (already allocated) gff object.  This object will be cleared, and filled with the elements from filter_gff which selected the fragments in the return value.  There will be a one-to-one corrrespondence between the elements in the return value and overlapping_frag.  
  @result New GFF that meets the criteria described by the above parameters.
  @note nonOverlapping and overlappingFragments cannot both be TRUE
  @note At least one of numbaseOverlap and percentOverlap needs to be > 0
*/
GFF_Set *gff_overlap_gff(GFF_Set *gff, GFF_Set *filter_gff,
			 int numbaseOverlap, 
			 double percentOverlap,
			 int nonOverlapping,
			 int overlappingFragments,
			 GFF_Set *overlapping_frag);

/** 
  Create a new GFF containing the features outside a region.
  @param gff GFF containing region0, will be used to find features outside region0
  @param region0 GFF containing features information of what to exclude from gff
  @result New GFF excluding features from region0
*/
GFF_Set *gff_inverse(GFF_Set *gff, GFF_Set *region0);

/**
  Create a new GFF where features are split.  
  @param maxlen Length to split by.  Array will be recycled to match the number of features in the gff.
  @param nmaxlen Length of maxlen array.
  @param drop Whether to include splits that are not exactly the size of maxlen
  @param splitFromRight Array of integer [0-1] for each split, TRUE indicates to split from the right, FALSE means split from the left.  Array will be recycled to match the number of features in GFF.
  @param splitFromRightLen Length of splitFromRight array.
  @result New GFF with split features
*/
GFF_Set *gff_split(GFF_Set *gff, int *maxlen, int nmaxlen,
		   int drop, int *splitFromRight, int splitFromRight_len);

/** 
  Create a GFF by thresholding an array of scores.
  @param seqname Sequence name to use in GFF
  @param firstIdx The coordinate of the first score (1-based)
  @param scores The array of scores to be thresholded; should represent a contiguous region starting with position firstIdx.
  @param numscore The length of the scores array
  @param threshold The threshold to use.  All coordinates with scores >= threshold will be included in resulting GFF.
  @param src The source string to use in the GFF
  @param featureName The name to give each feature
  @result A GFF_Set with features representing all regions with scores >= threshold.
*/
GFF_Set *gff_from_wig_threshold(char *seqname, int firstIdx, double *scores,
				int numscore, double threshold, char *src,
				char *featureName);

/** \} */
#endif
