/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: gff.h,v 1.16 2008-11-12 02:07:59 acs Exp $ */

/** \file gff.h
    Reading and writing of sequence features in General Feature Format (GFF). 
    Obeys file specification at
    http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml.
    \ingroup feature
*/


#ifndef GFF_H
#define GFF_H

#include <stringsplus.h>

#include <stdio.h>
#include <stdlib.h>

/** GFF_Feature structure.  Simply mirrors file-format spec (see
   http://www.sanger.ac.uk/resources/software/gff/spec.html) */ 
typedef struct {
  String *seqname;              /**< name of sequence (a single GFF file
                                   may describe multiple sequences) */
  String *source;               /**< source of feature -- usually
                                   program or database */
  String *feature;              /**< feature type; so far appears only
                                   semi-standardized.  One suggestion
                                   is to use the EMBL/DDBJ/GenBank
                                   feature table as a standard */
  int start, end;               /**< start and end positions.
                                   Convention is to start numbering
                                   with 1, and to make the range
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
                                   http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml. */
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
     http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml).
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

/* total number of columns */
#define GFF_NCOLS 9             
/* minimum allowable number of columns */
#define GFF_MIN_NCOLS 5         
/* starting number of features */
#define GFF_SET_START_SIZE 1000	

/* tags that identify meta-data in comments */
#define GFF_VERSION_TAG "gff-version"
#define GFF_SOURCE_VERSION_TAG "source-version"
#define GFF_DATE_TAG "date"

#define GFF_DEFAULT_VERSION 2

/* use when frame is null */
#define GFF_NULL_FRAME -1       

/* default feature types */
#define GFF_CDS_TYPE "CDS"
#define GFF_EXON_TYPE "exon"
#define GFF_INTRON_TYPE "intron"
#define GFF_START_TYPE "start_codon"
#define GFF_STOP_TYPE "stop_codon"
#define GFF_UTR5_TYPE "5'UTR"
#define GFF_UTR3_TYPE "3'UTR"
#define GFF_SPLICE5_TYPE "5'splice"
#define GFF_SPLICE3_TYPE "3'splice"


GFF_Set* gff_read_set(FILE *F);

GFF_Feature *gff_new_feature(String *seqname, String *source, String *feature,
                             int start, int end, double score, char strand, 
                             int frame, String *attribute, int score_is_null);

GFF_Feature *gff_new_feature_copy_chars(const char *seqname, const char *source,
					const char *feature,
					int start, int end, double score, 
					char strand, int frame, 
					const char *attribute,
					int score_is_null);

GFF_Feature *gff_new_feature_genomic_pos(String *position, String *source, 
                                         String *feature, double score, 
                                         int frame, String *attribute,
                                         int score_is_null);

GFF_Set *gff_new_set_len(int len);

GFF_Set *gff_new_set();

GFF_Set *gff_new_from_template(GFF_Set *gff);

GFF_Set *gff_new_set_init(char *source, char *source_version);

void gff_free_set(GFF_Set *set);

void gff_free_feature(GFF_Feature *feat);

void gff_print_set(FILE *F, GFF_Set *set);

void gff_print_feat(FILE *F, GFF_Feature *feat);

GFF_Feature *gff_new_feature_copy(GFF_Feature *orig);

GFF_Set *gff_copy_set_no_groups(GFF_Set *orig);

GFF_Set *gff_subset_range(GFF_Set *set, int startcol, int endcol, 
                          int reset_indices);

GFF_Set *gff_subset_range_overlap(GFF_Set *set, int startcol, int endcol);

GFF_Set *gff_subset_range_overlap_sorted(GFF_Set *set, int startcol, int endcol, 
					 int *startSearchIdx);

void gff_filter_by_type(GFF_Set *gff, List *types, int exclude, 
                        FILE *discards_f);

int gff_reverse_strand_only(List *features);

void gff_reverse_compl(List *features, int start_range, int end_range);

int gff_feature_comparator(const void* ptr1, const void* ptr2);

int gff_group_comparator(const void* ptr1, const void* ptr2);

void gff_sort(GFF_Set *set);

void gff_group(GFF_Set *set, char *tag);

void gff_group_by_feature(GFF_Set *set);

void gff_exon_group(GFF_Set *set, char *tag);

void gff_ungroup(GFF_Set *set);

GFF_Set* gff_read_from_fname(char *fname);

void gff_remove_overlaps(GFF_Set *gff, FILE *discards);

void gff_fix_start_stop(GFF_Set *gff);

void gff_absorb_helpers(GFF_Set *feats, List *primary_types, 
                        List *helper_types);

void gff_add_gene_id(GFF_Set *feats);

void gff_filter_by_group(GFF_Set *feats, List *groups);

void gff_create_utrs(GFF_Set *feats);

void gff_create_introns(GFF_Set *feats);

void gff_create_signals(GFF_Set *feats);

String *gff_group_name(GFF_Set *feats, GFF_Feature *f);

void gff_flatten(GFF_Set *feats);

void gff_flatten_within_groups(GFF_Set *feats);

//returns the total number of bases covered by feats,
//and flattens without regard to feature type, strand, or frame,
// but with regard to seqname only
int gff_flatten_mergeAll(GFF_Set *feats);

void gff_partition_by_type(GFF_Set *feats, List *types, List *subsets);

void gff_clear_set(GFF_Set *gff);

void gff_add_offset(GFF_Set *gff, int offset, int maxCoord);

GFF_Set *gff_overlap_gff(GFF_Set *gff, GFF_Set *filter_gff,
			 int numbase_overlap, 
			 double percentOverlap,
			 int nonOverlapping,
			 int overlappingFragments);

GFF_Set *gff_inverse(GFF_Set *gff, GFF_Set *region0);
GFF_Set *gff_split(GFF_Set *gff, int *maxlen, int nmaxlen,
		   int drop, int *splitFrom, int splitFrom_len);

#endif
