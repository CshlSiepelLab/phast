/* $Id: genepred.c,v 1.4 2004-06-22 22:14:06 acs Exp $
   Written by Adam Siepel, 2004
   Copyright 2004, Adam Siepel, University of California */

/** \file genepred.c
    Reading and writing of genepred files, used for UCSC Genome Browser.
    \ingroup feature
*/

#include <gff.h>
#include <assert.h>
#include <ctype.h>
#include <misc.h>
#include <hashtable.h>

#define GENEPRED_SOURCE "genepred"

/** Fill out a GFF_Set from a genepred file. */
void gff_read_from_genepred(GFF_Set *gff, FILE *F) {
  String *line = str_new(STR_LONG_LEN);
  List *l = lst_new_ptr(12), *tmpl1 = lst_new_ptr(10), 
    *tmpl2 = lst_new_ptr(10);
  int i, lineno = 0, frame = 0;
  Hashtable *hash = hsh_new(10000);

  while (str_readline(line, F) != EOF) {
    int txStart = 0, txEnd = 0, cdsStart = 0, cdsEnd = 0, 
      exonCount = 0, num = 0;
    String *name, *chrom, *tmpstr;
    char group[STR_MED_LEN];
    char strand;

    lineno++;

    str_trim(line);
    if (line->length == 0) continue;

    str_split(line, "\t", l);

    if (lst_size(l) < 10)
      die("ERROR (line %d): >= 10 columns required in genepred file.\n", lineno);

    name = lst_get_ptr(l, 0);
    chrom = lst_get_ptr(l, 1);

    tmpstr = lst_get_ptr(l, 2);
    
    if (tmpstr->length != 1 || 
        (tmpstr->chars[0] != '+' && tmpstr->chars[0] != '-'))
      die("ERROR (line %d): bad strand in genepred file: \"%s\".\n", 
          lineno, tmpstr->chars);

    strand = tmpstr->chars[0];

    if (str_as_int(lst_get_ptr(l, 3), &txStart) != 0 ||
        str_as_int(lst_get_ptr(l, 4), &txEnd) != 0 ||
        str_as_int(lst_get_ptr(l, 5), &cdsStart) != 0 ||
        str_as_int(lst_get_ptr(l, 6), &cdsEnd) != 0)
      die("ERROR (line %d): can't parse txStart, txEnd, cdsStart, or cdsEnd in genepred file.\n", lineno);

    txStart++; cdsStart++;      /* switch to GFF coord convention */

    if (cdsStart < txStart || cdsEnd > txEnd)
      die("ERROR (line %d): cds bounds outside of tx bounds in genepred file.\n", 
          lineno);

    if (str_as_int(lst_get_ptr(l, 7), &exonCount) != 0)
      die("ERROR (line %d): can't parse exonCount in genepred file.\n", 
          lineno);

    str_split(lst_get_ptr(l, 8), ",", tmpl1);
    str_split(lst_get_ptr(l, 9), ",", tmpl2);

    /* make sure group name is unique */
    if ((num = (int)hsh_get(hash, name->chars)) > 0) {
      num++;
      sprintf(group, "transcript_id \"%s.%d\"", name->chars, num);
      hsh_reset(hash, name->chars, (void*)num);
    }
    else {
      sprintf(group, "transcript_id \"%s\"", name->chars);
      hsh_put(hash, name->chars, (void*)1);
    }
    
    if (exonCount != lst_size(tmpl1) || lst_size (tmpl1) != lst_size(tmpl2))
      die("ERROR (line %d): exonStarts or exonEnds don't match exonCount in genepred file.\n", lineno);

    for (i = 0; i < exonCount; i++) {
      int eStart = 0, eEnd = 0;

      if (str_as_int(lst_get_ptr(tmpl1, i), &eStart) != 0 ||
              str_as_int(lst_get_ptr(tmpl2, i), &eEnd) != 0)
        die("ERROR (line %d): can't parse exonStarts or exonEnds in genepred file.\n", lineno);

      eStart++;

      /* create one feature for the whole exon and another for the CDS
         portion (if necessary) */
      lst_push_ptr(gff->features, 
                   gff_new_feature(str_dup(chrom), 
                                   str_new_charstr(GENEPRED_SOURCE), 
                                   str_new_charstr(GFF_EXON_TYPE), 
                                   eStart, eEnd, 0, strand, GFF_NULL_FRAME,
                                   str_new_charstr(group), TRUE));

      if ((eStart >= cdsStart && eStart <= cdsEnd) || /* left end in cds */
          (eEnd >= cdsStart && eEnd <= cdsEnd) || /* right end in cds */
          (cdsStart >= eStart && cdsEnd <= eEnd)) /* cds completely within exon */
        lst_push_ptr(gff->features, 
                     gff_new_feature(str_dup(chrom), 
                                     str_new_charstr(GENEPRED_SOURCE), 
                                     str_new_charstr(GFF_CDS_TYPE), 
                                     max(cdsStart, eStart), min(cdsEnd, eEnd), 
                                     0, strand, frame,
                                     str_new_charstr(group), TRUE));
    }

    /* FIXME: frame! */

    lst_free_strings(tmpl1);
    lst_free_strings(tmpl2);
    lst_free_strings(l);
  }

  str_free(line);
  lst_free(l);
  lst_free(tmpl1);
  lst_free(tmpl2);
  hsh_free(hash);
}

/** Write a GFF_Set in genepred format.  Features must already be
    grouped as desired.  Features will be sorted within groups as a
    side effect. */
void gff_print_genepred(FILE *OUTF,
                                /**< Output stream */
                        GFF_Set *feats
                                /**< Set to write */
                        ) {
  int i, j;
  String *exonStarts = str_new(STR_LONG_LEN);
  String *exonEnds = str_new(STR_LONG_LEN);

  if (feats->groups == NULL)
    die("ERROR: features must be grouped to print as genepred.\n");

  for (i = 0; i < lst_size(feats->groups); i++) {
    int cdsStart = -1, cdsEnd = -1, nexons = 0;
    char strand = '\0';
    String *seqname = NULL;
    GFF_FeatureGroup *g = lst_get_ptr(feats->groups, i);

    str_clear(exonStarts);
    str_clear(exonEnds);

    lst_qsort(g->features, gff_feature_comparator);

    for (j = 0; j < lst_size(g->features); j++) {
      GFF_Feature *f = lst_get_ptr(g->features, j);
      if (str_equals_charstr(f->feature, GFF_EXON_TYPE)) {
        str_append_int(exonStarts, f->start-1);
        str_append_char(exonStarts, ',');
        str_append_int(exonEnds, f->end);
        str_append_char(exonEnds, ',');
        if (nexons == 0) {
          strand = f->strand;
          seqname = f->seqname;
        }
        else if (strand != f->strand || !str_equals(seqname, f->seqname))
          die("ERROR (gff_print_genepred): inconsistent strand or seqname in GFF group \"%s\".\n", g->name->chars);
        nexons++;
      }
      else if (str_equals_charstr(f->feature, GFF_CDS_TYPE)) {
        if (cdsStart == -1) cdsStart = f->start-1;
        if (f->end > cdsEnd) cdsEnd = f->end;
      }
    }

    if (nexons > 0)
      fprintf(OUTF, "%s\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n", 
              g->name->chars, seqname->chars, strand, g->start-1, g->end, 
              cdsStart, cdsEnd, nexons, exonStarts->chars, exonEnds->chars);
  }

  str_free(exonStarts);
  str_free(exonEnds);
}
