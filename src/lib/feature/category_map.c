/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: category_map.c,v 1.15 2008-11-12 02:07:59 acs Exp $ */

/** \file category_map.c
    Mapping between feature types and site categories.  See
    category_map.h for details.
   \ingroup feature
*/

#include "category_map.h"
#include "gff.h"
#include "stacks.h"
#include "misc.h"
#include <hashtable.h>
#include <unistd.h>

//this has a conflict with RPHAST
#undef prec

static int *prec;

/** Read a CategoryMap from a file */
CategoryMap *cm_read(FILE *F) {
  String *line, *name;
  List *l;
  int cat, cat2, lineno, i, cm_read_error;
  CategoryMap *cm = NULL;
  CategoryRange *existing_range;
  static Regex *cat_range_re = NULL;
  static Regex *ncats_re = NULL;
  static Regex *fill_re = NULL;
  static Regex *label_re = NULL;
  static Regex *extend_re = NULL;
  int has_dependencies = 0;

  line = str_new(STR_SHORT_LEN);
  l = lst_new_ptr(3);
  if (cat_range_re == NULL) {
    cat_range_re = str_re_new("^[[:space:]]*([^[:space:]]+)[[:space:]]+([[:digit:]]+)(-([[:digit:]]+))?([[:space:]]+([[:digit:]].*))?"); 
    ncats_re = str_re_new("^[[:space:]]*NCATS[[:space:]]*=[[:space:]]*([[:digit:]]+)");
    fill_re = str_re_new("^[[:space:]]*FILL_PRECEDENCE[[:space:]]*=[[:space:]]*(.*)$");
    label_re = str_re_new("^[[:space:]]*LABELLING_PRECEDENCE[[:space:]]*=[[:space:]]*(.*)$");
    extend_re = str_re_new("^[[:space:]]*FEATURE_EXTEND[[:space:]]*:[[:space:]]*(.+)[[:space:]]*\\((.+)\\)$");
  }

  lineno = 0;
  while ((str_readline(line, F)) != EOF) {
    lineno++;
    str_trim(line);
    if (str_equals_charstr(line, ""))
      continue;

    if (str_re_match(line, ncats_re, l, 1) >= 0) { 
                                /* NCATS line */
      int ncats;
      str_as_int(lst_get_ptr(l, 1), &ncats);
      cm = cm_new(ncats);

      /* 0th category is "background" */
      cm->ranges[0] = 
        cm_new_category_range(str_new_charstr(BACKGD_CAT_NAME), 0, 0);
    }

    else if (cm == NULL || cm->ncats == 0) 
      die("ERROR: NCATS line must appear first, and must specify a positive number of categories.\n");

    else if (str_re_match(line, label_re, l, 1) >= 0) {               
                                /* LABELLING_PRECEDENCE line */
      List *tmpl = lst_new_ptr(cm->ncats);
      int tmpi;
      str_split((String*)lst_get_ptr(l, 1), " ,", tmpl);
      for (i = 0; i < lst_size(tmpl); i++) {
        String *s = (String*)lst_get_ptr(tmpl, i);
        if (str_as_int(s, &tmpi) != 0 || tmpi < 0 || tmpi > cm->ncats) 
          die("ERROR: bad integer in LABELLING_PRECEDENCE.\n");
        cm->labelling_precedence[tmpi] = i;
        str_free(s);
      }
      lst_free(tmpl);
    }

    else if (str_re_match(line, fill_re, l, 1) >= 0) {               
                                /* FILL_PRECEDENCE line */
      List *tmpl = lst_new_ptr(cm->ncats);
      int tmpi;
      str_split(lst_get_ptr(l, 1), " ,", tmpl);
      for (i = 0; i < lst_size(tmpl); i++) {
        String *s = lst_get_ptr(tmpl, i);
        if (str_as_int(s, &tmpi) != 0 || tmpi < 0 || tmpi > cm->ncats) 
          die("ERROR: bad integer in FILL_PRECEDENCE.\n");
        cm->fill_precedence[tmpi] = i;
        str_free(s);
      }
      lst_free(tmpl);
    }

    else if (str_re_match(line, extend_re, l, 2) >= 0) {
                                /* FEATURE_EXTEND line */
      String *target = lst_get_ptr(l, 2);
      List *sources = lst_new_ptr(2);
      str_split(lst_get_ptr(l, 1), " ,", sources);

      if (cm == NULL || (cat = cm_get_category(cm, target)) == 0)
        die("ERROR: FEATURE_EXTEND target must be a previously-defined non-background feature type.\n");

      for (i = 0; i < lst_size(sources); i++) {
        if (cm_get_category(cm, lst_get_ptr(sources, i)) == 0) 
          die("ERROR: FEATURE_EXTEND source list must consist of previously-defined non-background feature types.\n");
      }
    }

    else {                      /* 'range' line */
      if (str_re_match(line, cat_range_re, l, 6) < 0) 
        die("ERROR at line %d: '%s'\n", lineno, line->chars);

      name = str_dup((String*)lst_get_ptr(l, 1));
      str_as_int((String*)lst_get_ptr(l, 2), &cat);
/*       if (cat == 0) { */
/*         fprintf(stderr, "ERROR at line %d: category number must be a positive integer.\n",  */
/*                 lineno); */
/*         return NULL; */
/*       } */

      cat2 = cat;
      if (lst_get_ptr(l, 4) != NULL) 
        str_as_int((String*)lst_get_ptr(l, 4), &cat2);

      if (cat < 0 || cat2 < cat || cat2 > cm->ncats)
        die("ERROR: Illegal category range.\n");

      /* check for existing definitions of the specified category
         range.  Either no such definition must exist, or one must
         exist that spans exactly the same category numbers */
      existing_range = NULL;
      cm_read_error = 0;
      for (i = cat; !cm_read_error && i <= cat2; i++) {
        if (cm->ranges[i] != NULL && existing_range == NULL) 
          existing_range = cm->ranges[i];
        else if (cm->ranges[i] != existing_range)
          cm_read_error = 1;
      }
      if (existing_range != NULL && (existing_range->start_cat_no != cat || 
                                     existing_range->end_cat_no != cat2)) 
        cm_read_error = 1;

      if (cm_read_error) 
        die("ERROR: Overlapping category ranges.\n");

      /* either add new category range, or add new type to existing one */
      if (existing_range != NULL) {
        lst_push_ptr(existing_range->feature_types, name);
      }
      else {
        CategoryRange *cr = cm_new_category_range(name, cat, cat2);
        for (i = cat; i <= cat2; i++) 
          cm->ranges[i] = cr;
      }

      /* now address "conditioned_on" dependencies, if they have been
         specified */
      if (lst_get_ptr(l, 6) != NULL) {
        if (existing_range != NULL) 
          fprintf(stderr, "WARNING: ignoring 'conditioned on' list for type '%s'\n", 
                  name->chars);
        else {
          List *tmpl = lst_new_ptr(cm->ncats);
          int tmpi;
	  if (cm->conditioned_on[cat] != NULL)
	    die("ERROR cm_read: cm->conditioned_on[%i] should be NULL\n",
		cat);

          str_split((String*)lst_get_ptr(l, 6), " ,", tmpl);
          cm->conditioned_on[cat] = lst_new_int(lst_size(tmpl));
          for (i = cat + 1; i <= cat2; i++)
            cm->conditioned_on[i] = cm->conditioned_on[cat];
          /* all categories in range point to
             same "conditioned on" list */

          for (i = 0; i < lst_size(tmpl); i++) {
            String *s = (String*)lst_get_ptr(tmpl, i);
            if (str_as_int(s, &tmpi) != 0 || tmpi < 0 || tmpi > cm->ncats) 
              die("ERROR: bad integer in 'conditioned on' list for type '%s'.\n", 
                      name->chars);
            lst_push_int(cm->conditioned_on[cat], tmpi);
            str_free(s);
          }
          lst_free(tmpl);
          
          has_dependencies = 1;
        }
      }
    }

    for (i = 0; i < lst_size(l); i++)
      if (lst_get_ptr(l, i) != NULL)
        str_free((String*)lst_get_ptr(l, i));
  }

  /* make sure every category has been specified */
  for (i = 0; i <= cm->ncats; i++) 
    if (cm->ranges[i] == 0) 
      die("ERROR: category %d has not been specified.\n", i);

  /* build unspooler, if necessary */
  if (has_dependencies)
    cm->unspooler = cm_create_unspooler(cm->ncats + 1, cm->conditioned_on);

  str_free(line);
  lst_free(l);
  return cm;
}

/* for use in cm_print (see below) */
int compare_prec(const void* ptr1, const void* ptr2) {
  int val1 = *((int*)ptr1);
  int val2 = *((int*)ptr2);
  return (prec[val1] - prec[val2]);
}

/** Print a CategoryMap to a file */
void cm_print(CategoryMap *cm, FILE *F) {
  int i, j, k;
  List *tmpl;
  fprintf(F, "NCATS = %d\n\n", cm->ncats);

  for (i = 1; i <= cm->ncats; i++) {
    CategoryRange *cr = cm->ranges[i];
    for (j = 0; j < lst_size(cr->feature_types); j++) {
      String *s = (String*)lst_get_ptr(cr->feature_types, j);
      fprintf(F, "%-15s %d", s->chars, cr->start_cat_no);
      if (cr->end_cat_no > cr->start_cat_no)
        fprintf(F, "-%d", cr->end_cat_no);
      if (cm->conditioned_on[i] != NULL) {
        fprintf(F, "\t");
        for (k = 0; k < lst_size(cm->conditioned_on[i]); k++)
          fprintf(F, "%d%s", lst_get_int(cm->conditioned_on[i], k),
                  k + 1 == lst_size(cm->conditioned_on[i]) ? "" : ",");
      }
      fprintf(F, "\n");
    }
    i = cr->end_cat_no;         /* avoid looking multiple times at the
                                   same range */
  }

  /* reconstruct precedence lists */
  tmpl = lst_new_int(cm->ncats + 1);
  for (i = 0; i <= cm->ncats; i++) 
    lst_push_int(tmpl, i);
  prec = cm->labelling_precedence;
  lst_qsort(tmpl, compare_prec);
  fprintf(F, "\nLABELLING_PRECEDENCE = ");
  for (i = 0; i <= cm->ncats; i++) {
    int cat = lst_get_int(tmpl, i);
    if (cm->labelling_precedence[cat] != -1)
      fprintf(F, "%d%s", cat, i < cm->ncats ? "," : "");
  }
  fprintf(F, "\n");

  lst_clear(tmpl);
  for (i = 0; i <= cm->ncats; i++) 
    lst_push_int(tmpl, i);
  prec = cm->fill_precedence;
  lst_qsort(tmpl, compare_prec);
  fprintf(F, "FILL_PRECEDENCE = ");
  for (i = 0; i <= cm->ncats; i++) {
    int cat = lst_get_int(tmpl, i);
    if (cm->fill_precedence[cat] != -1)
      fprintf(F, "%d%s", cat, i < cm->ncats ? "," : "");
  }
  fprintf(F, "\n");
  lst_free(tmpl);
}

/** Return the 'base' category for a given feature type (first
    category in range).  Returns 0 (background) if no match. */
int cm_get_category(CategoryMap *cm, String *type) {
  int i, j;
  if (cm->ncats == 0) return 0;
  for (i = 0; i <= cm->ncats; i++) {
    CategoryRange *cr = cm->ranges[i];
    if (cr == NULL) return 0;   /* possible if called in cm_read
                                   before object is complete */
    for (j = 0; j < lst_size(cr->feature_types); j++) 
      if (str_equals(type, (String*)lst_get_ptr(cr->feature_types, j)))
        return i;
    i = cr->end_cat_no;         /* avoid looking multiple times at the
                                   same range */
  }
  return 0;
}

/** Return a list of category numbers corresponding to a given list of
   category names and or numbers.  */
List *cm_get_category_list(CategoryMap *cm, 
                                /**< CategoryMap object.  May be NULL
                                   if categories are specified by
                                   number rather than by name.  (In
                                   that case, this function will
                                   reduce to conversion of strings to
                                   integers.) */
                           List *names, 
                                /**< List of categories.  May be
                                   specified by name or number (useful
                                   when accepting input from users) */
                           int ignore_missing
                                /**< Whether to ignore unrecognized
                                   types.  If FALSE, then function
                                   will abort when it encounters an
                                   unrecognized type. */
                           ) {
  int i, j, cat;
  List *retval = lst_new_int(lst_size(names));
  for (i = 0; i < lst_size(names); i++) {
    String *n = lst_get_ptr(names, i);
    if (str_as_int(n, &cat) == 0) {
      if (cat < 0 || (cm != NULL && cat > cm->ncats)) 
        die("ERROR: category number %d is out of bounds.\n", cat);
      lst_push_int(retval, cat);
    }
    else {
      if (cm == NULL) 
        die("ERROR: if categories are specified by name, a category map is required.\n");
      cat = cm_get_category(cm, n);
      if (cat == 0 && !ignore_missing && !str_equals(n, cm_get_feature(cm, 0))) 
        die("ERROR: illegal category name (\"%s\")\n", n->chars);
      for (j = cm->ranges[cat]->start_cat_no; 
           j <= cm->ranges[cat]->end_cat_no; j++)
        lst_push_int(retval, j);
   }
  }
  return retval;
}


/** Return a list of category names corresponding to a given list of
   category names and or numbers.  Doesn't allocate new names,
   just pointers to Strings in the CategoryMap object or the
   provided List */
List *cm_get_category_str_list(CategoryMap *cm, 
			       /**< CategoryMap object.  May be NULL
				  if categories are specified by
				  name rather than by number.  (In
				  that case, this function will create
				  a copy of names with pointers to
				  its elements) */
			       List *names, 
			       /**< List of categories.  May be
				  specified by name or number (useful
				  when accepting input from users) */
			       int ignore_missing
			       /**< Whether to ignore unrecognized
				  types.  If FALSE, then function
				  will abort when it encounters an
				  unrecognized type. */
                           ) {
  int i, cat;
  List *retval = lst_new_ptr(lst_size(names));
  for (i = 0; i < lst_size(names); i++) {
    String *n = lst_get_ptr(names, i);
    if (str_as_int(n, &cat) == 0) {
      if (cm == NULL)
	die("ERROR: if categories are specified by number, a category map is required\n");
      if (cat < 0 || (cm != NULL && cat > cm->ncats)) 
        die("ERROR: category number %d is out of bounds.\n", cat);
      lst_push_ptr(retval, cm_get_feature(cm, cat));
    }
    else {
      if (cm != NULL) {
	cat = cm_get_category(cm, n);
	if (cat == 0 && !ignore_missing && !str_equals(n, cm_get_feature(cm, 0))) {
	  die("ERROR: illegal category name (\"%s\")\n", n->chars);
	}
	//return pointers to cm if possible
	lst_push_ptr(retval, cm_get_feature(cm, cat));
      }
      //otherwise return pointers to strings in list
      else lst_push_ptr(retval, n);
    }
  }
  return retval;
}


/** Retrieve the (primary) feature type associated with the specified
   category number.  Note: return value is passed by reference -- do
   not alter. */
String *cm_get_feature(CategoryMap *cm, int cat) {
  if (!(cat >= 0 && cat <= cm->ncats))
    die("ERROR cm_get_feature: cat=%i, should be in [0, %i]\n",
	cat, cm->ncats);
  return lst_get_ptr(cm->ranges[cat]->feature_types, 0);
}

/** Obtain a unique name based on the (primary) feature type
   associated with the specified category number.  Returns a pointer
   to a newly allocated String.*/
String *cm_get_feature_unique(CategoryMap *cm, int cat) {
  String *retval = str_dup(lst_get_ptr(cm->ranges[cat]->feature_types, 0));
  if (! (cat >= 0 && cat <= cm->ncats))
    die("ERROR cm_get_feature_unique: cat=%i, should be in [0,%i]\n",
	cat, cm->ncats);
  if (cm->ranges[cat]->start_cat_no != cm->ranges[cat]->end_cat_no) {
    str_append_char(retval, '-');
    str_append_int(retval, cat - cm->ranges[cat]->start_cat_no + 1);
  }    
  return retval;
}

/* return list of category names corresponding to list of category
   numbers */
List *cm_get_features(CategoryMap *cm, List *catnos) {
  int mark[cm->ncats+1];
  List *retval = lst_new_ptr(lst_size(catnos));
  int i, cat;
  for (i = 0; i <= cm->ncats; i++) mark[i] = 0;
  for (i = 0; i < lst_size(catnos); i++) {
    cat = lst_get_int(catnos, i);
    if (!mark[cm->ranges[cat]->start_cat_no]) {
      lst_push_ptr(retval, cm_get_feature(cm, cat));
      mark[cm->ranges[cat]->start_cat_no] = 1;
    }
  }
  return retval;
}

CategoryMap *cm_new(int ncats) {
  int i;
  CategoryMap *cm = (CategoryMap*)smalloc(sizeof(CategoryMap));
  cm->ncats = ncats;
  cm->ranges = smalloc((cm->ncats + 1) * sizeof(CategoryRange*));
  cm->conditioned_on = smalloc((cm->ncats + 1) * sizeof(List*));
  cm->labelling_precedence = smalloc((cm->ncats + 1) * sizeof(int));
  cm->fill_precedence = smalloc((cm->ncats + 1) * sizeof(int));
  for (i = 0; i <= cm->ncats; i++) {
    cm->ranges[i] = NULL;
    cm->conditioned_on[i] = NULL;
    cm->labelling_precedence[i] = -1;
    cm->fill_precedence[i] = -1;
  }
  cm->unspooler = NULL;

  return cm;
}

CategoryMap *cm_create_copy(CategoryMap *src) {
  int i, has_dependencies = 0;
  CategoryMap *retval = cm_new(src->ncats);
  for (i = 0; i <= src->ncats; i++) {
    retval->ranges[i] = cm_category_range_create_copy(src->ranges[i]);
    if (retval->conditioned_on[i] != NULL) {
      retval->conditioned_on[i] = lst_new_int(lst_size(src->conditioned_on[i]));
      lst_cpy(retval->conditioned_on[i], src->conditioned_on[i]);
      has_dependencies = 1;
    }
    retval->labelling_precedence[i] = src->labelling_precedence[i];
    retval->fill_precedence[i] = src->fill_precedence[i];
  }
  if (has_dependencies) cm_create_unspooler(retval->ncats + 1, 
                                            retval->conditioned_on);
  return retval;
}

/** Create a trivial CategoryMap, with feature types equal to category
   numbers (plus an optional prefix) and ranges all of size one.  */
CategoryMap* cm_create_trivial(int ncats, char *feature_prefix) {
  int i;
  CategoryMap *retval = cm_new(ncats);
  for (i = 0; i <= ncats; i++) {
    String *type = str_new(STR_SHORT_LEN);
    if (feature_prefix != NULL) str_append_charstr(type, feature_prefix);
    str_append_int(type, i);
    retval->ranges[i] = cm_new_category_range(type, i, i);
  }
  return retval;
}

/** Create a category map with a category for each feature type in a
    GFF_Set.  Category numbers are assigned in order of appearance of
    types */
CategoryMap* cm_new_from_features(GFF_Set *feats) {
  int i;
  CategoryMap *retval;
  Hashtable *hash;
  List *types;

  /* first scan features for all types */
  hash = hsh_new(10);
  types = lst_new_ptr(10);
  for (i = 0; i < lst_size(feats->features); i++) {
    GFF_Feature *f = lst_get_ptr(feats->features, i);
    checkInterruptN(i, 10000);
    if (hsh_get(hash, f->feature->chars) == (void*)-1) {
      lst_push_ptr(types, f->feature);
      hsh_put_int(hash, f->feature->chars, 1);
    }
  }
  hsh_free(hash);

  /* now create a simple category map */
  retval = cm_new(lst_size(types));
  for (i = 0; i <= retval->ncats; i++) {
    String *type = i == 0 ? str_new_charstr(BACKGD_CAT_NAME) : 
      str_dup(lst_get_ptr(types, i-1));
    retval->ranges[i] = cm_new_category_range(type, i, i);
  }
  lst_free(types);
  return retval;
}

/** Create a new category map from a string that can
    either be a filename or a brief "inlined" category map, e.g.,
    "NCATS = 3 ; CDS 1-3".  Useful for command-line arguments. */
CategoryMap* cm_new_string_or_file(const char *optarg) {
  int i;
  char fname[STR_SHORT_LEN];
  String *str = str_new_charstr(optarg);
  FILE *F;
  CategoryMap *retval = NULL;

  str_double_trim(str);
  if (str_starts_with_charstr(str, "NCATS")) {
    /* replace semicolons with carriage returns */
    for (i = 0; i < str->length; i++)
      if (str->chars[i] == ';') str->chars[i] = '\n';

    /* we'll just dump a little tmp file and read it with cm_read */
    sprintf(fname, "cm.tmp.%d", getpid());
    F = fopen_fname(fname, "w+");
    fprintf(F, "%s", str->chars);
    fclose(F);
    retval = cm_read(fopen_fname(fname, "r"));
    unlink(fname);
  }
  else {
    F = fopen_fname(str->chars, "r");
    retval = cm_read(F);
    fclose(F);
  }

  str_free(str);
  return retval;
}

/** Reallocate a category map to allow for the specified size. */
void cm_realloc(CategoryMap *cm, int new_ncats) {
  int i;
  int old_ncats = cm->ncats;
  cm->ncats = new_ncats;
  cm->ranges = srealloc(cm->ranges, (cm->ncats + 1) * sizeof(CategoryRange*));
  cm->conditioned_on = srealloc(cm->conditioned_on,
                                (cm->ncats + 1) * sizeof(List*));
  cm->labelling_precedence = srealloc(cm->labelling_precedence,
                                      (cm->ncats + 1) * sizeof(int));
  cm->fill_precedence = srealloc(cm->fill_precedence,
                                 (cm->ncats + 1) * sizeof(int));
  for (i = old_ncats+1; i <= cm->ncats; i++) {
    cm->ranges[i] = NULL;
    cm->conditioned_on[i] = NULL;
    cm->labelling_precedence[i] = -1;
    cm->fill_precedence[i] = -1;
  }
}

/** Free memory associated with category map. */
void cm_free(CategoryMap *cm) {
  int i;
  for (i = 0; i <= cm->ncats; i++) {
    int len = 0;
    if (cm->ranges[i] != NULL) {
      len = cm->ranges[i]->end_cat_no - cm->ranges[i]->start_cat_no;
      cm_free_category_range(cm->ranges[i]);
    }
    if (cm->conditioned_on[i] != NULL)
      lst_free(cm->conditioned_on[i]);
    i += len;
  }
  sfree(cm->ranges);
  sfree(cm->conditioned_on);
  sfree(cm->labelling_precedence);
  sfree(cm->fill_precedence);

  if (cm->unspooler != NULL)
    cm_free_unspooler(cm->unspooler);

  sfree(cm);
  return;
}

CategoryRange* cm_new_category_range(String *type, int start_cat_no,
                                     int end_cat_no) {
  CategoryRange *cr = (CategoryRange*)smalloc(sizeof(CategoryRange));
  cr->feature_types = lst_new_ptr(1);
  lst_push_ptr(cr->feature_types, type);
  cr->start_cat_no = start_cat_no;
  cr->end_cat_no = end_cat_no;
  return cr;
}

CategoryRange* cm_category_range_create_copy(CategoryRange *src) {
  int i;
  CategoryRange *retval = 
    cm_new_category_range(str_dup(lst_get_ptr(src->feature_types, 0)),
                          src->start_cat_no, src->end_cat_no);
  for (i = 1; i < lst_size(src->feature_types); i++)
    lst_push_ptr(retval->feature_types, 
                 str_dup(lst_get_ptr(src->feature_types, i)));
  return retval;
}

void cm_free_category_range(CategoryRange *cr) {
  int i;
  for (i = 0; i < lst_size(cr->feature_types); i++) {
    String *s = (String*)lst_get_ptr(cr->feature_types, i);
    if (s != NULL) str_free(s);
  }
  lst_free(cr->feature_types);
  sfree(cr);
}

/** Add new feature to CategoryMap.  Assumes a feature of the
   specified type does not already exist. */
void cm_add_feature_type(CategoryMap *cm, String *type, int cycle_size) {
  int catstart = cm->ncats + 1;
  int catend = cm->ncats + cycle_size;
  int i;
  CategoryRange *cr = cm_new_category_range(type, catstart, catend);
  cm->ncats = catend;
  cm->ranges = (CategoryRange**)srealloc(cm->ranges, (cm->ncats + 1) *
                                        sizeof(CategoryRange*));
  cm->labelling_precedence = (int*)srealloc(cm->labelling_precedence, 
                                           (cm->ncats + 1) * sizeof(int));
  cm->fill_precedence = (int*)srealloc(cm->fill_precedence, 
                                      (cm->ncats + 1) * sizeof(int));
  cm->conditioned_on = (List**)srealloc(cm->conditioned_on, 
                                      (cm->ncats + 1) * sizeof(List*));
  for (i = catstart; i <= catend; i++) {
    cm->ranges[i] = cr;
    cm->labelling_precedence[i] = -1;
    cm->fill_precedence[i] = -1;
    cm->conditioned_on[i] = NULL;
  }
}

/** Create a GFF_Set from a sequence of category/state numbers, using
   a specified category map and mapping from raw state numbers to
   category numbers.  */
GFF_Set *cm_labeling_as_gff(CategoryMap *cm, 
                                /**< CategoryMap to use in mapping */
                            int *path, 
                                /**< Raw sequence of state/category numbers  */
                            int length, 
                                /**< Length of sequence */
                            int *path_to_cat, 
                                /**< Mapping from raw numbers to
                                   category numbers */
                            int *reverse_compl, 
                                /**< Array of boolean values indicating
                                   whether each raw-sequence value
                                   corresponds to the reverse strand  */
                            char *seqname, 
                                /**< char string to use as 'seqname' in
                                   generated GFF_Set  */
                            char *source, 
                                /**< char string to use as 'source' in
                                   generated GFF_Set  */
                            List *frame_cats, 
                                /**< Categories for which to obtain frame
                                   information (by name) */
                            char *grouptag,
                            /**< Tag to use to define groups in
                               GFF_Set (e.g., "transcript_id") */
                            char *idpref
                                /**< Prefix for ids of predicted
                                   elements (may be NULL).  Can be
                                   used to ensure ids are unique. */
                            ) {
  int beg, end, i, cat, frame, groupno;
  GFF_Set *gff = gff_new_set_init("PHAST", PHAST_VERSION);
  int do_frame[cm->ncats+1];
  char strand;
  char groupstr[STR_SHORT_LEN];
  int ignore_0 = str_equals_charstr(cm_get_feature(cm, 0), BACKGD_CAT_NAME);
                                /* ignore category 0 if background  */

  if (length <= 0) return gff;

  for (i = 0; i <= cm->ncats; i++) do_frame[i] = 0;
  if (frame_cats != NULL)
    for (i = 0; i < lst_size(frame_cats); i++) {
      int cat = cm_get_category(cm, lst_get_ptr(frame_cats, i));
      if (cat != 0)             /* ignore background or unrecognized name */
        do_frame[cat] = 1;
    }

  groupno = 1;
  if (idpref != NULL)
    sprintf(groupstr, "%s \"%s.%d\"", grouptag != NULL ? grouptag : "id", 
            idpref, groupno);
  else
    sprintf(groupstr, "%s \"%d\"", grouptag != NULL ? grouptag : "id", groupno);

  i = 0;
  while (i < length) {
    checkInterruptN(i, 10000);
    cat = cm->ranges[path_to_cat[path[i]]]->start_cat_no;
    strand = reverse_compl[path[i]] ? '-' : '+';
    frame = do_frame[cat] ? path_to_cat[path[i]] - cat : GFF_NULL_FRAME;

    /* scan ahead until enter new category range (or reach end of seq) */
    beg = i + 1;                /* begin of feature (GFF coords) */
    for (i++; i < length && 
           cm->ranges[path_to_cat[path[i]]]->start_cat_no == cat; i++);
    end = i;                    /* end of feature (GFF coords) */

    /* if minus strand, adjust frame to reflect end */
    if (strand == '-' && do_frame[cat]) 
      frame = path_to_cat[path[i-1]] - cat;

    /* if legitimate feature (non-background), then incorp into GFF_Set */
    if (cat != 0 || !ignore_0)  /* create new feature and add */
      lst_push_ptr(gff->features, 
                   gff_new_feature(str_new_charstr(seqname), 
                                   str_new_charstr(source), 
                                   str_dup(cm_get_feature(cm, cat)), 
                                   beg, end, 0, strand, frame, 
                                   str_new_charstr(groupstr), TRUE));

    if (cat == 0 && beg > 1) {
      groupno++;                /* increment group number each time a
                                   sequence of 0s is encountered  */
      if (idpref != NULL)
        sprintf(groupstr, "%s \"%s.%d\"", grouptag != NULL ? grouptag : "id", 
                idpref, groupno);
      else
        sprintf(groupstr, "%s \"%d\"", grouptag != NULL ? grouptag : "id", 
                groupno);
    }
  }

  return gff;
}

/** maps a sequence (array) of category numbers from the spooled space to
   the unspooled space, using the current unspooler.  Original
   sequence is overwritten */
void cm_spooled_to_unspooled(CategoryMap *cm, int *path, int pathlen) {
  int j, sp_state, prev_sp_state;
  List *pred;

  if (cm->unspooler == NULL) return;

  pred = lst_new_int(cm->unspooler->nstates_spooled);
  prev_sp_state = -1;
  for (j = 0; j < pathlen; j++) {
    if (!(path[j] >= 0 && path[j] <= cm->unspooler->nstates_spooled))
      die("ERROR cm_spooled_to_unspooled: path[%i]=%i, should be in [0, %i]\n",
	  j, path[j], cm->unspooler->nstates_spooled);

    sp_state = path[j];
    path[j] = cm_get_unspooled_state(cm, path[j], pred);

    if (path[j] == -1) 
      die("ERROR: failure mapping to uspooled state at position %d.\n", j);

    if (sp_state != prev_sp_state) {
      /* if the current (spooled) state is not conditioned on any
         other state, then its predecessor cannot matter, so the list
         can be cleared */
      if (lst_size(cm->unspooler->spooled_to_unspooled[sp_state]->children) == 0)
        lst_clear(pred);

      lst_push_int(pred, sp_state);
    }

    prev_sp_state = sp_state;
  }

  lst_free(pred);
}

/* maps a sequence (array) of category numbers from the unspooled space to
   the spooled space, using the current unspooler.  Original
   sequence is overwritten */
void cm_unspooled_to_spooled(CategoryMap *cm, int *path, int pathlen) {
  int j;
  if (cm->unspooler == NULL) return;
  for (j = 0; j < pathlen; j++) {
    if (!(path[j] >= 0 && path[j] <= cm->unspooler->nstates_unspooled))
      die("ERROR cm_unspooled_to_spooled: path[%i]=%i, should be in [0,%i]\n",
	  j, path[j], cm->unspooler->nstates_unspooled);
    path[j] = cm->unspooler->unspooled_to_spooled[path[j]];
  }
}

/* returns the spooled category number corresponding to a given
   unspooled category number; simplifies calling code by taking care
   of cases with and without unspooler */
int cm_unspooled_to_spooled_cat(CategoryMap *cm, int spooled_cat) {
  return (cm->unspooler != NULL ? 
          cm->unspooler->unspooled_to_spooled[spooled_cat] : 
          spooled_cat);
}

/* conditioned_on must be an array of integer lists; specifically, the
   ith element must be the list of state numbers on which the ith
   state is conditioned. */
Unspooler *cm_create_unspooler(int nstates_spooled, List **conditioned_on) {
  UnspoolNode *n;
  int i, j;
  Stack *s;
  Unspooler *unsp;
  int *mark;
  int capacity;

  unsp = (Unspooler*)smalloc(sizeof(Unspooler));
  unsp->nstates_spooled = nstates_spooled;
  unsp->nstates_unspooled = 0;
  unsp->spooled_to_unspooled = 
    (UnspoolNode**)smalloc(nstates_spooled * sizeof(UnspoolNode*));
  capacity = nstates_spooled * nstates_spooled;
  unsp->unspooled_to_spooled = (int*)smalloc(capacity * sizeof(int));

  mark = (int*)smalloc(nstates_spooled * sizeof(int));
  s = stk_new_ptr(nstates_spooled);
  for (i = 0; i < nstates_spooled; i++) {
    /* erase marks (used to detect cycles) */
    for (j = 0; j < nstates_spooled; j++) mark[j] = 0;

    unsp->spooled_to_unspooled[i] = cm_new_unspool_node(i);
    stk_push_ptr(s, unsp->spooled_to_unspooled[i]);
    while ((n = (UnspoolNode*)stk_pop_ptr(s)) != NULL) {
      if (conditioned_on[n->oldstate] == NULL ||
          lst_size(conditioned_on[n->oldstate]) == 0) {
        n->newstate = unsp->nstates_unspooled++;

        /* mapping to spooled space */
        if (n->newstate >= capacity) {
          capacity *= 2;
          unsp->unspooled_to_spooled = 
            (int*)srealloc(unsp->unspooled_to_spooled, 
                          capacity * sizeof(int));          
        }
        unsp->unspooled_to_spooled[n->newstate] = i;
      }
      else {
        for (j = 0; j < lst_size(conditioned_on[n->oldstate]); j++) {
          int oldstate = lst_get_int(conditioned_on[n->oldstate], j);
          UnspoolNode *m;

          if (mark[oldstate] == 1)
            die("ERROR: cycle in 'conditioned_on' dependencies.\n");
          mark[oldstate] = 1;

          m = cm_new_unspool_node(oldstate);
          lst_push_ptr(n->children, m);
          stk_push_ptr(s, m);
        }
      }
    }
  }
  stk_free(s);
  sfree(mark);
  return unsp;
}

UnspoolNode *cm_new_unspool_node(int oldstate) {
  UnspoolNode *uspn = (UnspoolNode*)smalloc(sizeof(UnspoolNode));
  uspn->oldstate = oldstate;
  uspn->newstate = -1;
  uspn->children = lst_new_ptr(3);
  return uspn;
}

void cm_free_unspool_node(UnspoolNode *n) {
  lst_free(n->children);
  sfree(n);
}

void cm_free_unspooler(Unspooler *unsp) {
  int i, j;
  UnspoolNode *n;
  Stack *s = stk_new_ptr(unsp->nstates_unspooled);

  /* free all nodes in tree */
  for (i = 0; i < unsp->nstates_spooled; i++) 
    stk_push_ptr(s, unsp->spooled_to_unspooled[i]);
  while ((n = stk_pop_ptr(s)) != NULL) {
    for (j = 0; j < lst_size(n->children); j++)
      stk_push_ptr(s, lst_get_ptr(n->children, j));
    cm_free_unspool_node(n);
  }

  sfree(unsp->spooled_to_unspooled);
  sfree(unsp->unspooled_to_spooled);
  stk_free(s);
  sfree(unsp);
}

/* last item in predecessors is assumed to be the most recently visited */
int cm_get_unspooled_state(CategoryMap *cm, int spooled_state, 
                           List *predecessors) {
  UnspoolNode *n, *child;
  int p, pred_idx, i;

  pred_idx = lst_size(predecessors) - 1;
  n = cm->unspooler->spooled_to_unspooled[spooled_state];
  while (n->newstate == -1) {
    child = NULL;
    while (n != child && pred_idx >= 0) {
      p = lst_get_int(predecessors, pred_idx--);
      for (i = 0; n != child && i < lst_size(n->children); i++) {
        child = (UnspoolNode*)lst_get_ptr(n->children, i);
        if (child->oldstate == p) 
          n = child;
      }
    }
    if (n != child) {
      fprintf(stderr, "ERROR (cm_get_unspooled_state): no match for state %d preceded by state(s) ", spooled_state);
      for (i = 0; i < lst_size(predecessors); i++) 
        fprintf(stderr, "%d ", lst_get_int(predecessors, i));
      fprintf(stderr, "\n");
      return -1;
    }
    if (n != child)
      die("ERROR cm_get_unspooled_state: n != child\n");
  }
  return n->newstate;
}

/* simple wrapper for cm_read that opens specified filename and reads
   a category map from it; aborts with error message if unable to open
   the file.  Saves typing in mains for command-line programs */
CategoryMap* cm_read_from_fname(char *fname) {
  CategoryMap *cm = NULL;
  FILE *F;
  if ((F = fopen(fname, "r")) == NULL || 
      (cm = cm_read(F)) == NULL) 
    die("ERROR: cannot read category map from %s.\n", fname);
  fclose(F);
  return cm;
}

/* given list of spooled category names/numbers, return a list of
   corresponding unspooled category numbers */
List *cm_get_unspooled_list(CategoryMap *cm, List *spooled) {
  List *spooled_catnos, *unspooled_catnos;
  int mark[cm->ncats+1];
  int i;

  spooled_catnos = cm_get_category_list(cm, spooled, 0);
  if (cm->unspooler == NULL) return spooled_catnos;

  unspooled_catnos = lst_new_int(lst_size(spooled_catnos) * 3);
  for (i = 0; i <= cm->ncats; i++) mark[i] = 0;
  for (i = 0; i < lst_size(spooled_catnos); i++) 
    mark[lst_get_int(spooled_catnos, i)] = 1;

  for (i = 0; i < cm->unspooler->nstates_unspooled; i++) 
    if (mark[cm->unspooler->unspooled_to_spooled[i]])
      lst_push_int(unspooled_catnos, i);

  lst_free(spooled_catnos);
  return unspooled_catnos;
}
