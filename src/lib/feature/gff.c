/* $Id: gff.c,v 1.1.1.1 2004-06-03 22:43:11 acs Exp $
   Written by Adam Siepel, Summer 2002
   Copyright 2002, Adam Siepel, University of California */

/** \file gff.c
    Reading and writing of sequence features in General Feature Format
    (GFF).  Obeys file specification at
    http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml. 
    \ingroup feature
*/
   
#include <gff.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <hashtable.h>
#include <misc.h>
#include <ctype.h>
#include <bed.h>

/** Read a set of features from a file and return a newly allocated
   GFF_Set object.  Function reads until end-of-file is encountered or
   error occurs (aborts on error).  Comments and blank lines are
   ignored and special "meta-data" comments are parsed (see
   http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml).  Only
   the first five columns of feature lines are considered required
   ('name', 'source', 'feature', 'start', and 'end'); subsequent
   fields are optional ('score', 'strand', 'frame', and 'attribute').
   Default value for score, strand, and frame is null ('.') and for
   attribute is the empty string ('').  Columns must be separated by
   tabs.  */
GFF_Set* gff_read_set(FILE *F) {
  int start, end, frame, score_is_null, lineno, done_with_header = 0;
  double score;
  char strand;
  String *attr, *line;
  GFF_Feature *feat;
  GFF_Set *set;
  List *l, *substrs;
  static Regex *spec_comment_re = NULL;

  line = str_new(STR_LONG_LEN);
  set = gff_new_set();
  l = lst_new_ptr(GFF_NCOLS);
  substrs = lst_new_ptr(4);

  lineno = 0;
  while (str_readline(line, F) != EOF) {
    lineno++;

    str_double_trim(line);
    if (line->length == 0) continue;

    if (!done_with_header && str_starts_with_charstr(line, "##")) {
      int unrecognized = 0;
      if (spec_comment_re == NULL)
        spec_comment_re = str_re_new("^[[:space:]]*##[[:space:]]*([^[:space:]]+)[[:space:]]+([^[:space:]]+)([[:space:]]+([^[:space:]]+))?");
 
      if (str_re_match(line, spec_comment_re, substrs, 4) >= 0) {
        String *tag, *val1, *val2;
        tag = (String*)lst_get_ptr(substrs, 1);
        val1 = (String*)lst_get_ptr(substrs, 2);
        val2 =  lst_size(substrs) > 4 ? (String*)lst_get_ptr(substrs, 4) : NULL;
        
        if (str_equals_nocase_charstr(tag, GFF_VERSION_TAG))
          str_cpy(set->gff_version, val1);
        else if (str_equals_nocase_charstr(tag, GFF_SOURCE_VERSION_TAG) && 
                 val2 != NULL) {
          str_cpy(set->source, val1);
          str_cpy(set->source_version, val2); 
        }
        else if (str_equals_nocase_charstr(tag, GFF_DATE_TAG))
          str_cpy(set->date, val1);
        else unrecognized = 1;
      }
      else unrecognized = 1;

      if (unrecognized)
        fprintf(stderr, "WARNING: unrecognized meta-data: '%s'\n", 
                line->chars);
      lst_free_strings(substrs);

      continue;
    }

    else if (line->chars[0] == '#') continue; /* just skip ordinary comments */
    
    done_with_header = 1;

    str_split(line, "\t", l);

    /* if first record, check to see if the file's a BED; if there are
       between 3 and 8 or 12 columns, and if the 2nd and 3rd columns are
       integers, then we'll try reading it as a BED. */
    if (lst_size(set->features) == 0 && 
        ((lst_size(l) >= 3 && lst_size(l) <= 8) || lst_size(l) == 12) && 
        str_as_int(lst_get_ptr(l, 1), &start) == 0 && 
        str_as_int(lst_get_ptr(l, 2), &end) == 0) {
      rewind(F);
      gff_read_from_bed(set, F);
      break;
    }

    /* set defaults for optional fields */
    strand = '.';
    frame = GFF_NULL_FRAME;
    score_is_null = 1;            

    if (lst_size(l) < GFF_MIN_NCOLS) 
      die("ERROR at line %d (gff_read_set): minimum of %d columns are required.\n", 
          lineno, GFF_MIN_NCOLS);
  
    if (str_as_int(lst_get_ptr(l, 3), &start) != 0) 
      die("ERROR at line %d (gff_read_set): non-numeric 'start' value ('%s').\n", 
          lineno, ((String*)lst_get_ptr(l, 3))->chars);

    if (str_as_int((String*)lst_get_ptr(l, 4), &end) != 0) 
      die("ERROR at line %d (gff_read_set): non-numeric 'end' value ('%s').\n", 
          lineno, ((String*)lst_get_ptr(l, 4))->chars);

    if (lst_size(l) > 5) {
      String *score_str = (String*)lst_get_ptr(l, 5);
      if (! str_equals_charstr(score_str, ".")) {
        if (str_as_dbl(score_str, &score) != 0) 
          die( "ERROR at line %d (gff_read_set): non-numeric and non-null 'score' value ('%s').\n", 
                  lineno, score_str->chars);
        else 
          score_is_null = 0;
      }
    }

    strand = '.';
    if (lst_size(l) > 6) {
      String *tmp = (String*)lst_get_ptr(l, 6);
      if (tmp->length == 0 || tmp->length > 1 || 
          (tmp->chars[0] != '+' && tmp->chars[0] != '-' && 
           tmp->chars[0] != '.')) 
        die("ERROR at line %d: illegal 'strand' ('%s').\n", 
            lineno, tmp->chars);
      strand = tmp->chars[0];
    }

    frame = GFF_NULL_FRAME;
    if (lst_size(l) > 7) {
      String *tmp = (String*)lst_get_ptr(l, 7);
      if (! str_equals_charstr(tmp, ".") != 0) {
        if (str_as_int(tmp, &frame) != 0 || frame < 0 || frame > 2) 
          die("ERROR at line %d: illegal 'frame' ('%s').\n", 
              lineno, tmp->chars);
      }
    }

    if (lst_size(l) > 8) 
      attr = str_dup(lst_get_ptr(l, 8));
    else attr = str_new(0);
    
    feat = gff_new_feature(str_dup(lst_get_ptr(l, 0)), 
                           str_dup(lst_get_ptr(l, 1)),
                           str_dup(lst_get_ptr(l, 2)), start, end, score, 
                           strand, frame, attr, score_is_null);

    lst_push_ptr(set->features, feat);
    lst_free_strings(l);
  }

  str_free(line);
  lst_free(l);
  lst_free(substrs);
  return set;
}

/** Create new GFF_Feature object with specified attributes.  Warning:
   strings are copied by reference (for efficiency).
   Returns newly allocated GFF_Feature object. */
GFF_Feature *gff_new_feature(String *seqname, String *source, String *feature,
                             int start, int end, double score, char strand, 
                             int frame, String *attribute, 
                             int score_is_null) {
  GFF_Feature *feat = (GFF_Feature*)smalloc(sizeof(GFF_Feature));

  assert(seqname != NULL && source != NULL && feature != NULL && 
         attribute != NULL && 
         (strand == '+' || strand == '-' || strand == '.') &&
         (frame == GFF_NULL_FRAME || (0 <= frame && frame <=2)));

  feat->seqname = seqname;
  feat->source = source;
  feat->feature = feature;
  feat->start = start;
  feat->end = end;
  feat->score = score;
  feat->strand = strand;
  feat->frame = frame;
  feat->attribute = attribute;
  feat->score_is_null = score_is_null;
  return feat;
}

/** Create a new GFF_Feature representing a genomic position, based on
   a position string obeying the conventions of the UCSC genome
   browser, e.g., chr10:102553847-102554897.  If a trailing '+' or '-'
   appears in the string, then it is interpreted as the strand;
   otherwise the null strand is used.  The function returns NULL if
   the string cannot be parsed. */
GFF_Feature *gff_new_feature_genomic_pos(String *position, String *source, 
                                         String *feature, double score, 
                                         int frame, String *attribute,
                                         int score_is_null) {
  GFF_Feature *retval = NULL;
  List *substrs = lst_new_ptr(4);
  static Regex *posre = NULL;
  if (posre == NULL) 
    posre = str_re_new("(chr[_a-zA-Z0-9]+):([0-9]+)-([0-9]+)([-+])?");

  if (str_re_match(position, posre, substrs, 4) >= 3) {
    int start, end;
    char strand= '.';
    String *chr = str_dup(lst_get_ptr(substrs, 1)), 
      *tmpstr = lst_get_ptr(substrs, 4);
    str_as_int(lst_get_ptr(substrs, 2), &start);
    str_as_int(lst_get_ptr(substrs, 3), &end);
    if (tmpstr != NULL) strand = tmpstr->chars[0];
    return gff_new_feature(chr, source, feature, start, end, score, 
                           strand, frame, attribute, score_is_null);
  }
  lst_free_strings(substrs);
  lst_free(substrs);
  return retval;
}

/** Create new GFF_Set object.  Attributes must be set directly.
   Returns newly allocated GFF_Set object. */
GFF_Set *gff_new_set() {
  GFF_Set *set = (GFF_Set*)smalloc(sizeof(GFF_Set));
  set->features = lst_new_ptr(GFF_SET_START_SIZE);
  set->gff_version = str_new(STR_SHORT_LEN);
  set->source = str_new(STR_SHORT_LEN);
  set->source_version = str_new(STR_SHORT_LEN);
  set->date = str_new(STR_SHORT_LEN);
  set->groups = NULL;
  return set;
}

/** Create new GFF_Set object, initializing with same version, source,
    etc.\ as a template GFF_Set object */
GFF_Set *gff_new_from_template(GFF_Set *gff) {
  GFF_Set *retval = gff_new_set();
  str_cpy(retval->gff_version, gff->gff_version);
  str_cpy(retval->source, gff->source);
  str_cpy(retval->source_version, gff->source_version);
  str_cpy(retval->date, gff->date); 
  return retval;
}

/** Create new GFF_Set object, using ordinary defaults and other
   parameters as specified.  Sets gff version to '2' and date to current
   date, and sets source and source version as specified.
   Returns newly allocated GFF_Set object. */
GFF_Set *gff_new_set_init(char *source, char *source_version) {
  struct tm *tm;
  time_t t;
  GFF_Set *set = gff_new_set();
  str_append_int(set->gff_version, GFF_DEFAULT_VERSION);
  str_cpy_charstr(set->source, source);
  str_cpy_charstr(set->source_version, source_version);

  t = time(NULL);
  tm = localtime(&t);
  str_append_int(set->date, 1900 + tm->tm_year);
  str_append_charstr(set->date, "-");
  str_append_int(set->date, tm->tm_mon + 1);
  str_append_charstr(set->date, "-");
  str_append_int(set->date, tm->tm_mday);

  return set;
}

/** Free resources associated with GFF_Set object (including all
   features and the set object itself). */ 
void gff_free_set(GFF_Set *set) {
  int i;
  if (set->features != NULL) {
    for (i = 0; i < lst_size(set->features); i++)
      gff_free_feature((GFF_Feature*)lst_get_ptr(set->features, i));
    lst_free(set->features);
  }
  str_free(set->gff_version);
  str_free(set->source);
  str_free(set->source_version);
  str_free(set->date);
  if (set->groups != NULL) gff_ungroup(set);
  free(set);
}

/** Free resources associated with GFF_Feature object.  Strings will be
   freed even if they were allocated externally. */
void gff_free_feature(GFF_Feature *feat) {
  str_free(feat->seqname);
  str_free(feat->source);
  str_free(feat->feature);
  str_free(feat->attribute);
  free(feat);
}

/** Output a GFF_Set to the specified stream in GFF. */
void gff_print_set(FILE *F, GFF_Set *set) {
  int i;

  if (set->gff_version->length > 0)
    fprintf(F, "##%s %s\n", GFF_VERSION_TAG, set->gff_version->chars);

  if (set->source_version->length > 0)
    fprintf(F, "##%s %s %s\n", GFF_SOURCE_VERSION_TAG, set->source->chars, 
            set->source_version->chars);

  if (set->date->length > 0)
    fprintf(F, "##%s %s\n", GFF_DATE_TAG, set->date->chars);

  if (set->features != NULL)
    for (i = 0; i < lst_size(set->features); i++) 
      gff_print_feat(F, (GFF_Feature*)lst_get_ptr(set->features, i));
}

/** Print an individual GFF_Feature object as a GFF line. */
void gff_print_feat(FILE *F, GFF_Feature *feat) {
  char score_str[50], frame_str[50];

  if (feat->score_is_null) strcpy(score_str, ".");
  else sprintf(score_str, "%f", feat->score);

  if (feat->frame == GFF_NULL_FRAME) strcpy(frame_str, ".");
  else sprintf(frame_str, "%d", feat->frame);

  fprintf(F, "%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%s\n", 
          feat->seqname->chars, feat->source->chars, feat->feature->chars, 
          feat->start, feat->end, score_str, feat->strand, frame_str, 
          feat->attribute->chars);    
}

/** Create an exact copy of a GFF_Feature object */
GFF_Feature *gff_new_feature_copy(GFF_Feature *orig) {
  String *seqname = str_dup(orig->seqname);
  String *source = str_dup(orig->source);
  String *feature = str_dup(orig->feature);
  String *attribute = str_dup(orig->attribute);
  return gff_new_feature(seqname, source, feature, orig->start, orig->end, 
                         orig->score, orig->strand, orig->frame, attribute, 
                         orig->score_is_null);
}

/** Create a new GFF_Set representing the features in a particular
    coordinate range.  Keeps features such that feat->start >= startcol
    and feat->end <= endcol. */
GFF_Set *gff_subset_range(GFF_Set *set, int startcol, int endcol, 
                          int reset_indices) {
  GFF_Set *subset = gff_new_set();
  int i;

  str_cpy(subset->gff_version, set->gff_version);
  str_cpy(subset->source, set->source);
  str_cpy(subset->source_version, set->source_version);
  str_cpy(subset->date, set->date); /* make current date instead? */

  /* Note: uses linear search */
  for (i = 0; i < lst_size(set->features); i++) {
    GFF_Feature *feat = (GFF_Feature*)lst_get_ptr(set->features, i);
    if (feat->start >= startcol && feat->end <= endcol) {
      GFF_Feature *newfeat = gff_new_feature_copy(feat);
      if (reset_indices) {
        newfeat->start = newfeat->start - startcol + 1;
        newfeat->end = newfeat->end - startcol + 1;
      }
      lst_push_ptr(subset->features, newfeat);
    }
  }
  return subset;
}

/** Discard any features whose feature type is *not* in the specified
    list.  Objects in list should be of type String. */
void gff_filter_by_type(GFF_Set *gff, List *include) {
  List *newfeats = lst_new_ptr(lst_size(gff->features));
  int i;
  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *f = lst_get_ptr(gff->features, i);
    if (str_in_list(f->feature, include)) /* use linear search -- hash
                                             probably not worth the
                                             overhead */
      lst_push_ptr(newfeats, f);
    else gff_free_feature(f);
  }
  lst_free(gff->features);
  gff->features = newfeats;
}

/** Tests whether an entire GFF_Set appears to refer to the reverse
   strand.  Returns 1 if no features have strand equal to '+' and at
   least one has strand equal to '-'; otherwise returns 0. */
int gff_reverse_strand_only(GFF_Set *set) {
  int i, impossible = 0, possible = 0;
  for (i = 0; !impossible && i < lst_size(set->features); i++) {
    GFF_Feature *feat = (GFF_Feature*)lst_get_ptr(set->features, i);
    if (feat->strand == '-') 
      possible = 1;
    else if (feat->strand == '+')
      impossible = 1;
  }
  if (impossible || !possible)
    return 0;
  return 1;
}

/** Converts all features of specified GFF_Set to positive strand
   (assuming they currently represent the negative strand).  Adjusts
   positions and strand, and also reverses order of appearance.
   start_range and end_range should indicate the first and last
   coordinates of the range covered by the set (that is, the
   coordinates are inclusive).  This function assumes all feature
   coordinates are in the same frame of reference.
*/
void gff_reverse_compl(GFF_Set *set, int start_range, int end_range) {
  int i;
  for (i = 0; i < lst_size(set->features); i++) {
    GFF_Feature *feat = (GFF_Feature*)lst_get_ptr(set->features, i);
    int tmp = feat->start;
    feat->start = end_range - feat->end + start_range;
    feat->end = end_range - tmp + start_range;
    if (feat->strand == '-') feat->strand = '+';
  }
  /* also reverse order of features (will generally be in ascending order) */
  for (i = 0; i < lst_size(set->features)/2; i++) {
    GFF_Feature *tmp = (GFF_Feature*)lst_get_ptr(set->features, i);
    lst_set_ptr(set->features, i, 
                lst_get_ptr(set->features, lst_size(set->features)-i-1));
    lst_set_ptr(set->features, lst_size(set->features)-i-1, tmp);
  }
}

/* used by gff_sort (see below) */
int gff_feature_start_comparator(const void* ptr1, const void* ptr2) {
  GFF_Feature *feat1 = *((GFF_Feature**)ptr1);
  GFF_Feature *feat2 = *((GFF_Feature**)ptr2);
  if (feat1->start != feat2->start) 
    return (feat1->start - feat2->start);
  return (feat1->end - feat2->end);
                                /* note that this rule has the effect
                                   of putting short features that
                                   overlap the ends of longer ones in
                                   a sensible order */
}

/* used by gff_sort (see below) */
int gff_group_comparator(const void* ptr1, const void* ptr2) {
  GFF_FeatureGroup *group1 = *((GFF_FeatureGroup**)ptr1);
  GFF_FeatureGroup *group2 = *((GFF_FeatureGroup**)ptr2);
  GFF_Feature *feat1, *feat2;
  if (lst_size(group1->features) == 0 || 
      lst_size(group2->features) == 0) 
    return 0; 
                                /* in this case, order doesn't matter */

  feat1 = lst_get_ptr(group1->features, 0);
  feat2 = lst_get_ptr(group2->features, 0);
  if (feat1->start != feat2->start) 
    return (feat1->start - feat2->start);
  return (feat1->end - feat2->end);
}

/** Sorts features in ascending order by start position.  If features
    are grouped (see gff_group), then they will be sorted within
    groups, and groups will be sorted by start position of first
    feature */
void gff_sort(GFF_Set *set) {
  int i, j;
  if (set->groups == NULL)
    lst_qsort(set->features, gff_feature_start_comparator);
  else {
    for (i = 0; i < lst_size(set->groups); i++)
      lst_qsort(((GFF_FeatureGroup*)lst_get_ptr(set->groups, i))->features, 
                gff_feature_start_comparator);
    lst_qsort(set->groups, gff_group_comparator);
    /* now reorder features according to groups */
    lst_clear(set->features);
    for (i = 0; i < lst_size(set->groups); i++) {
      GFF_FeatureGroup *group = lst_get_ptr(set->groups, i);
      for (j = 0; j < lst_size(group->features); j++)
        lst_push_ptr(set->features, lst_get_ptr(group->features, j));
    }
  }
}

/** Group features by value of specified tag.  All features with
    undefined values will be placed in a single group. */
void gff_group(GFF_Set *set, String *tag) {
  char tmpstr[STR_SHORT_LEN];
  Regex *tag_re;
  List *l = lst_new_ptr(1);
  int est_no_groups = max(lst_size(set->features) / 10, 1);
  Hashtable *hash = hsh_new(est_no_groups);
  String *nullstr = str_new(1); /* empty string represents missing or
                                   null value for tag */
  int i;

  if (set->groups != NULL)
    gff_ungroup(set);

  set->groups = lst_new_ptr(est_no_groups);

  /* since we only use the 'attribute' field for grouping, we'll store
     it unparsed, and parse it only when we need to group */

  sprintf(tmpstr, ".*%s[[:space:]]+(\"[^\"]*\"|[^[:space:]])", tag->chars);
  tag_re = str_re_new(tmpstr);

  for (i = 0; i < lst_size(set->features); i++) {
    GFF_Feature *f = lst_get_ptr(set->features, i);
    String *val = nullstr;
    GFF_FeatureGroup *group;
    lst_clear(l);

    if (f->attribute->length > tag->length && /* avoid checking empty
                                                 or null attrs */
        str_re_match(f->attribute, tag_re, l, 1) >= 0) {
      val = lst_get_ptr(l, 1);
      str_remove_quotes(val);
    }

    if ((group = hsh_get(hash, val->chars)) == (void*)-1) {
      group = smalloc(sizeof(GFF_FeatureGroup));
      group->name = str_dup(val);
      group->features = lst_new_ptr(5);
      lst_push_ptr(set->groups, group);
      hsh_put(hash, val->chars, group);
    }
    lst_push_ptr(group->features, f);

    lst_free_strings(l);
  }

  lst_free(l);
  str_re_free(tag_re);
  str_free(nullstr);
  hsh_free(hash);
}

/** Remove grouping of features */
void gff_ungroup(GFF_Set *set) {
  int i;
  if (set->groups == NULL) return;
  for (i = 0; i < lst_size(set->groups); i++) {
    GFF_FeatureGroup *group = lst_get_ptr(set->groups, i);
    str_free(group->name);
    lst_free(group->features);
    free(group);
  }
  lst_free(set->groups);
  set->groups = NULL;
}

/** Groups features into sets that define contiguous blocks of
    coordinates, e.g., for an exon and adjacent splice sites within a
    gene.  If features have already been grouped (e.g., by transcript
    id), then creates subgroups by adding new tags, with tag name as
    specified, and tag value composed of the tag value for the "outer"
    group (e.g., the transcript_id) and a unique suffix indicating the
    "inner" group.  This function sorts the features as a side effect. */ 
void gff_exon_group(GFF_Set *set, String *tag) {
  List *groups;
  int i, j;
  char tmpstr[STR_SHORT_LEN];
  GFF_FeatureGroup *dummy = NULL;

  gff_sort(set);

  if (set->groups == NULL) {
    /* create dummy list consisting of one group, to simplify code
       below */
    dummy = smalloc(sizeof(GFF_FeatureGroup));
    dummy->name = NULL;
    dummy->features = set->features;
    groups = lst_new_ptr(1);
    lst_push_ptr(groups, dummy);
  }
  else 
    groups = set->groups;

  for (i = 0; i < lst_size(groups); i++) {
    GFF_FeatureGroup *group = lst_get_ptr(groups, i);

    int lastend = -1, idx = 0;
    for (j = 0; j < lst_size(group->features); j++) {
      GFF_Feature *f = lst_get_ptr(group->features, j);
      if (f->start > lastend + 1) idx++;

      if (f->attribute->length == 0 || str_equals_charstr(f->attribute, "."))
        str_clear(f->attribute);
      else 
        str_append_charstr(f->attribute, " ; ");

      if (group->name == NULL || group->name->length == 0)
        sprintf(tmpstr, "%s \"%d\"", tag->chars, idx);
      else 
        sprintf(tmpstr, "%s \"%s.%d\"", tag->chars, group->name->chars, idx);
      
      str_append_charstr(f->attribute, tmpstr);

      if (f->end > lastend)
        lastend = f->end;
    }
  }

  /* now regroup using tag for subgroups */
  gff_group(set, tag);

  if (dummy != NULL) {
    free(dummy);
    lst_free(groups);
  }
}

/** Creates a list of GFF_Sets, corresponding to the different "groups"
   in the given GFF_Set.  Two features are considered to belong to the
   same group if their "attribute" fields are identical.  Features
   with null attributes (denoted by empty strings or ".") are handled
   specially: they are considered part of a group iff they fall
   between two explicit members of the group (this is useful when
   groups are used to defined "ranges" of an alignment; see
   msa_reverse_compl_gff).  Order is important; calling code should
   sort if necessary.  Note: subset_list must already be initialized
   with lst_new_ptr.  This function creates new GFF_Set objects and
   new copies of all features.  */
void gff_partition_by_group(GFF_Set *gff, List *subset_list) {
  int i, j;
  String *thisgroup = NULL;
  GFF_Set *subset = NULL;
  List *nulls;

  lst_clear(subset_list);

  if (lst_size(gff->features) == 0) return;

  nulls = lst_new_ptr(max(lst_size(gff->features)/4, 10));

  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *thisfeat = 
      gff_new_feature_copy((GFF_Feature*)lst_get_ptr(gff->features, i));
    String *attr = thisfeat->attribute;

    if (attr->length == 0 || str_equals_charstr(attr, ".")) {
      lst_push_ptr(nulls, thisfeat); /* null group -- don't add
                                        immediately; wait to see if
                                        it's bracketed */
      continue;
    }

    /* feature must have non-null group */
    if (subset == NULL || thisgroup == NULL || 
        !str_equals(attr, thisgroup)) { /* new group */

      /* first put outstanding null features into their own subset */
      if (lst_size(nulls) > 0) {
        subset = gff_new_from_template(gff);      
        for (j = 0; j < lst_size(nulls); j++) 
          lst_push_ptr(subset->features, lst_get_ptr(nulls, j));
        lst_push_ptr(subset_list, subset);
        lst_clear(nulls);
      }

      /* now start a subset for the new group */
      subset = gff_new_from_template(gff);      
      lst_push_ptr(subset_list, subset);

      thisgroup = attr;
    }

    /* flush the buffer of null features if necessary, then add the
       new feature */
    if (lst_size(nulls) > 0) {
      for (j = 0; j < lst_size(nulls); j++) {
        GFF_Feature *nullfeat = lst_get_ptr(nulls, j);
        str_cpy(nullfeat->attribute, thisgroup);
        lst_push_ptr(subset->features, nullfeat);
      }
      lst_clear(nulls);
    }

    lst_push_ptr(subset->features, thisfeat);
  }

  /* finally, make a new group of any remaining null features */
  if (lst_size(nulls) > 0) {
    subset = gff_new_from_template(gff);      
    for (j = 0; j < lst_size(nulls); j++) 
      lst_push_ptr(subset->features, lst_get_ptr(nulls, j));
    lst_push_ptr(subset_list, subset);
  }

  lst_free(nulls);
}

/** Partition the features of a GFF such that each partition either has
   a single feature of the specified type and flanking "helper"
   features (e.g., a cds with accompanying splice sites or start/stop
   codons) or has no feature of the specified type.  The list
   'partitions' will be populated with lists of GFF_Features.  As
   above, order is important.  Here it is particularly important that
   a feature's helpers are adjacent to it in the list (as well as
   ordered appropriately wrt it). */
void gff_partition_by_feature(GFF_Set *gff, List *partitions, List *features,
                              List *helpers) {
  enum {FEATURE, HELPER, OTHER} thistype, lasttype;
  int *helper_map = smalloc(lst_size(gff->features) * sizeof(int));
  GFF_Feature *feat, *f;
  int i, j;
  List *thispart;

  lst_clear(partitions);

  /* first determine whether each (potential) helper is associated with a
     previous feature (-1), a subsequent feature (+1), or no feature (0) */
  for (i = 0; i < lst_size(gff->features); i++) helper_map[i] = 0;
  for (i = 0; i < lst_size(gff->features); i++) {
    int start, end;
    feat = lst_get_ptr(gff->features, i);
    if (str_in_list(feat->feature, features)) {
      /* extend in both directions as far as possible */
      start = feat->start;
      for (j = i - 1; j >= 0; j--) {
        f = lst_get_ptr(gff->features, j);
        if (!str_in_list(f->feature, helpers) || 
            (f->end != start - 1 && f->start != start)) 
                                /* allow for overlapping helpers */
          break;
        start = f->start;
        helper_map[j] = +1;
      }
      end = feat->end;
      for (j = i + 1; j < lst_size(gff->features); j++) {
        f = lst_get_ptr(gff->features, j);
        if (!str_in_list(f->feature, helpers) || 
            (f->start != end + 1 && f->end != end)) 
                                /* allow for overlapping helpers */
          break;
        end = f->end;
        helper_map[j] = -1;
      }
    }
  }

  /* now split */
  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *feat = lst_get_ptr(gff->features, i);

    if (str_in_list(feat->feature, features)) thistype = FEATURE;
    else if (str_in_list(feat->feature, helpers) && helper_map[i] != 0) 
      thistype = HELPER;        /* treat an unassociated HELPER like
                                   an OTHER */
    else thistype = OTHER;

    if (i == 0 || 
        (thistype == OTHER && lasttype != OTHER) ||
        (thistype == FEATURE && !(lasttype == HELPER && helper_map[i-1] == 1)) ||
        (thistype == HELPER && 
         (lasttype == OTHER ||
          (lasttype == FEATURE && helper_map[i] != -1) ||
          (lasttype == HELPER && helper_map[i] != helper_map[i-1])))) {
      /* start a new partition */
      thispart = lst_new_ptr(lst_size(features) + lst_size(helpers)); /* approx */
      lst_push_ptr(partitions, thispart);
    }

    lst_push_ptr(thispart, feat);    
    lasttype = thistype;
  }
  free(helper_map);
}

/** Simple wrapper for gff_read_set that opens specified filename and reads
   a GFF_Set from it; aborts with error message if unable to open
   the file.  Saves typing in mains for command-line programs */
GFF_Set* gff_read_from_fname(char *fname) {
  GFF_Set *gff = NULL;
  FILE *F;
  if ((F = fopen(fname, "r")) == NULL || 
      (gff = gff_read_set(F)) == NULL) {
    fprintf(stderr, "ERROR: cannot read annotations from %s.\n", fname);
    exit(1);
  }
  fclose(F);
  return gff;
}

/** Identify overlapping groups of features and remove all but the
   first instance.  Features of the same group are allowed to overlap
   (often done for convenience, e.g., with start/stop codons or splice
   sites, then resolved via labelling precedence).  Features with
   attribute "." are ignored (all are kept). */
void gff_remove_overlaps_by_group(GFF_Set *gff) {

  List *group_starts = lst_new_int(lst_size(gff->features));
  List *group_ends = lst_new_int(lst_size(gff->features));
  List *new_feats = lst_new_ptr(lst_size(gff->features));
  int last_gend = -1;
  int i, j, discard, gstart, gend;
  String *gname;
  GFF_Feature *feat;
  List *group = lst_new_ptr(50);
  
  if (lst_size(gff->features) == 0) return;

  feat = lst_get_ptr(gff->features, 0);
  for (i = 0; i < lst_size(gff->features); ) {
    discard = 0;

    gname = feat->attribute;

    if (str_equals_charstr(gname, ".")) {
      lst_push_ptr(new_feats, feat);
      feat = lst_get_ptr(gff->features, ++i);
      continue;
    }

    /* get boundaries of group */
    gstart = feat->start;
    gend = feat->end;
    lst_push_ptr(group, feat);
    for (i++; i < lst_size(gff->features); i++) {
      feat = lst_get_ptr(gff->features, i);
      if (!str_equals(feat->attribute, gname))
        break;
      if (feat->start < gstart) gstart = feat->start;
      if (feat->end > gend) gend = feat->end;
      lst_push_ptr(group, feat);
    }
    /* note: on exit, feat is always the first item in the *next*
       group, unless the end of the list has been reached */

    /* check for overlap */
    if (gstart > last_gend) {   /* common case, has to be safe */
      lst_push_int(group_starts, gstart);
      lst_push_int(group_ends, gend);
      last_gend = gend;
    }
    else {                  /* have to search list */
      int group_list_idx = lst_bsearch_int(group_starts, gstart);
      int prev_end = group_list_idx >= 0 ? 
        lst_get_int(group_ends, group_list_idx) : -1;
      int next_start = group_list_idx+1 < lst_size(group_starts) ?
        lst_get_int(group_starts, group_list_idx+1) : gend+1;         
      if (prev_end >= gstart || next_start <= gend) {
        fprintf(stderr, "WARNING: group '%s' (%d-%d) overlaps a previous group -- ignoring.\n", gname->chars, gstart, gend);
        discard = 1;
      }
      else {
        lst_insert_idx_int(group_starts, group_list_idx, gstart);
        lst_insert_idx_int(group_ends, group_list_idx, gend);
        if (gend > last_gend) last_gend = gend;
      }
    }
    
    for (j = 0; j < lst_size(group); j++) {
      GFF_Feature *f = lst_get_ptr(group, j);
      if (discard) 
        gff_free_feature(f);
      else
        lst_push_ptr(new_feats, f);
    }
    lst_clear(group);
  }

  lst_free(gff->features);
  gff->features = new_feats;

  lst_free(group);
  lst_free(group_starts);
  lst_free(group_ends);
}

/** Identify overlapping features and remove all but the first
   instance.  This function differs from the above in that features
   are not grouped.  If 'types' is non-NULL, then only features whose
   types are listed are checked for overlap (features of other types are
   always kept).  */
void gff_remove_overlaps(GFF_Set *gff, List *types) {
  int i, last_end = -1;
  List *starts = lst_new_int(lst_size(gff->features));
  List *ends = lst_new_int(lst_size(gff->features));
  List *new_feats = lst_new_ptr(lst_size(gff->features));

  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *feat = lst_get_ptr(gff->features, i);
    
    if (types != NULL && !str_in_list(feat->feature, types)) {
      lst_push_ptr(new_feats, feat);
      continue;
    }

    /* check for overlap */
    if (feat->start > last_end) {   /* common case, has to be safe */
      lst_push_int(starts, feat->start);
      lst_push_int(ends, feat->end);
      last_end = feat->end;
      lst_push_ptr(new_feats, feat);
    }
    else {                  /* have to search list */
      int list_idx = lst_bsearch_int(starts, feat->start);
      int prev_end = list_idx >= 0 ? 
        lst_get_int(ends, list_idx) : -1;
      int next_start = list_idx+1 < lst_size(starts) ?
        lst_get_int(starts, list_idx+1) : feat->end+1;         
      if (prev_end >= feat->start || next_start <= feat->end) {
        fprintf(stderr, "WARNING: feature '%s (%s)' (%d-%d) overlaps a previous feature -- ignoring.\n", feat->feature->chars, feat->attribute->chars, feat->start, feat->end);
        gff_free_feature(feat);
      }
      else {
        lst_insert_idx_int(starts, list_idx, feat->start);
        lst_insert_idx_int(ends, list_idx, feat->end);
        if (feat->end > last_end) last_end = feat->end;
        lst_push_ptr(new_feats, feat);
      }
    }
  }

  lst_free(gff->features);
  gff->features = new_feats;

  lst_free(starts);
  lst_free(ends);
}

/** Adjust coords of terminal cds exons such that stop codons are included */
void gff_fix_stops(GFF_Set *gff, String* cds_feat_type, String *stop_feat_type) {
  int i, j;
  List *partitions = lst_new_ptr(lst_size(gff->features));
  List *features = lst_new_ptr(1);
  List *helpers = lst_new_ptr(1);
  lst_push_ptr(features, cds_feat_type);
  lst_push_ptr(helpers, stop_feat_type);
  gff_partition_by_feature(gff, partitions, features, helpers);
  for (i = 0; i < lst_size(partitions); i++) {
    List *part = lst_get_ptr(partitions, i);
    GFF_Feature *cds = NULL, *stop = NULL;
    for (j = 0; j < lst_size(part); j++) {
      GFF_Feature *f = lst_get_ptr(part, j);
      if (str_equals(f->feature, cds_feat_type)) cds = f;
      else if (str_equals(f->feature, stop_feat_type)) stop = f;
    }
    if (cds != NULL && stop != NULL) {
      if (cds->strand == '+' && cds->end == stop->start-1)
        cds->end = stop->end;
      else if (cds->strand == '-' && cds->start == stop->end+1)
        cds->start = stop->start;
    }
  }
  lst_free(partitions); lst_free(features); lst_free(helpers);
}

