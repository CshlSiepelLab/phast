/* $Id: gff.h,v 1.3 2004-06-14 22:52:16 acs Exp $
   Written by Adam Siepel, Summer 2002
   Copyright 2002, Adam Siepel, University of California */

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
   http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml) */ 
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
                                   in files) */
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


GFF_Set* gff_read_set(FILE *F);

GFF_Feature *gff_new_feature(String *seqname, String *source, String *feature,
                             int start, int end, double score, char strand, 
                             int frame, String *attribute, int score_is_null);

GFF_Feature *gff_new_feature_genomic_pos(String *position, String *source, 
                                         String *feature, double score, 
                                         int frame, String *attribute,
                                         int score_is_null);

GFF_Set *gff_new_set();

GFF_Set *gff_new_from_template(GFF_Set *gff);

GFF_Set *gff_new_set_init(char *source, char *source_version);

void gff_free_set(GFF_Set *set);

void gff_free_feature(GFF_Feature *feat);

void gff_print_set(FILE *F, GFF_Set *set);

void gff_print_feat(FILE *F, GFF_Feature *feat);

GFF_Feature *gff_new_feature_copy(GFF_Feature *orig);

GFF_Set *gff_subset_range(GFF_Set *set, int startcol, int endcol, 
                          int reset_indices);

void gff_filter_by_type(GFF_Set *gff, List *include, FILE *discards_f);

int gff_reverse_strand_only(List *features);

void gff_reverse_compl(List *features, int start_range, int end_range);

void gff_sort(GFF_Set *set);

void gff_group(GFF_Set *set, char *tag);

void gff_exon_group(GFF_Set *set, char *tag);

void gff_ungroup(GFF_Set *set);

GFF_Set* gff_read_from_fname(char *fname);

void gff_remove_overlaps(GFF_Set *gff, FILE *discards);

void gff_fix_stops(GFF_Set *gff, char* cds_type, char *stop_type);

#endif
