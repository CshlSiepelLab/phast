/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: genepred.c,v 1.7 2008-11-12 02:07:59 acs Exp $ */

/** \file genepred.c
    Reading and writing of genepred files, used for UCSC Genome Browser.
    \ingroup feature
*/

#include <gff.h>
#include <ctype.h>
#include <misc.h>
#include <hashtable.h>
#include <string.h>

#define GENEPRED_SOURCE "genepred"

/** Fill out a GFF_Set from a genepred file. */
void gff_read_from_genepred(GFF_Set *gff, FILE *F) {
  String *line = str_new(STR_LONG_LEN);
  List *l = lst_new_ptr(12), *tmpl1 = lst_new_ptr(10), 
    *tmpl2 = lst_new_ptr(10), *framefeats = lst_new_ptr(50);
  int i, lineno = 0, hasbin=-1;
  Hashtable *hash = hsh_new(10000);

  while (str_readline(line, F) != EOF) {
    int txStart = 0, txEnd = 0, cdsStart = 0, cdsEnd = 0, 
      exonCount = 0, num = 0;
    String *name, *chrom, *tmpstr;
    GFF_Feature *f;
    char group[STR_MED_LEN];
    char strand;

    checkInterruptN(lineno, 1000);
    lineno++;
    
    if (line->chars[0] == '#') continue;

    str_trim(line);
    if (line->length == 0) continue;

    str_split(line, "\t", l);

    if (lst_size(l) < 10)
      die("ERROR (line %d): >= 10 columns required in genepred file.\n", lineno);

    if (hasbin == -1) {  //determine if this genepred has a "bin" column
      if (lst_size(l) >= 11 &&
	  str_as_int(lst_get_ptr(l, 0), &i)==0 &&
	  str_as_int(lst_get_ptr(l, 4), &i)==0 &&
	  str_as_int(lst_get_ptr(l, 5), &i)==0 &&
	  str_as_int(lst_get_ptr(l, 6), &i)==0 &&
	  str_as_int(lst_get_ptr(l, 7), &i)==0) {
	tmpstr = lst_get_ptr(l, 3);
	if (str_equals_charstr(tmpstr, "+") ||
	    str_equals_charstr(tmpstr, "-"))
	  hasbin = 1;
      }
      if (hasbin == -1) hasbin=0;
    }

    name = lst_get_ptr(l, 0+hasbin);
    chrom = lst_get_ptr(l, 1+hasbin);

    tmpstr = lst_get_ptr(l, 2+hasbin);
    
    if (tmpstr->length != 1 || 
        (tmpstr->chars[0] != '+' && tmpstr->chars[0] != '-'))
      die("ERROR (line %d): bad strand in genepred file: \"%s\".\n", 
          lineno, tmpstr->chars);

    strand = tmpstr->chars[0];

    if (str_as_int(lst_get_ptr(l, 3+hasbin), &txStart) != 0 ||
        str_as_int(lst_get_ptr(l, 4+hasbin), &txEnd) != 0 ||
        str_as_int(lst_get_ptr(l, 5+hasbin), &cdsStart) != 0 ||
        str_as_int(lst_get_ptr(l, 6+hasbin), &cdsEnd) != 0)
      die("ERROR (line %d): can't parse txStart, txEnd, cdsStart, or cdsEnd in genepred file.\n", lineno);

    txStart++; cdsStart++;      /* switch to GFF coord convention */

    if (cdsStart < txStart || cdsEnd > txEnd)
      die("ERROR (line %d): cds bounds outside of tx bounds in genepred file.\n", 
          lineno);

    if (str_as_int(lst_get_ptr(l, 7+hasbin), &exonCount) != 0)
      die("ERROR (line %d): can't parse exonCount in genepred file.\n", 
          lineno);

    str_split(lst_get_ptr(l, 8+hasbin), ",", tmpl1);
    str_split(lst_get_ptr(l, 9+hasbin), ",", tmpl2);

    /* make sure group name is unique */
    if ((num = hsh_get_int(hash, name->chars)) > 0) {
      num++;
      sprintf(group, "transcript_id \"%s.%d\"", name->chars, num);
      hsh_reset_int(hash, name->chars, num);
    }
    else {
      sprintf(group, "transcript_id \"%s\"", name->chars);
      hsh_put_int(hash, name->chars, 1);
    }
    
    if (exonCount != lst_size(tmpl1) || lst_size (tmpl1) != lst_size(tmpl2))
      die("ERROR (line %d): exonStarts or exonEnds don't match exonCount in genepred file.\n", lineno);

    lst_clear(framefeats);
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
          (cdsStart >= eStart && cdsEnd <= eEnd)) { /* cds completely within exon */
        f = gff_new_feature(str_dup(chrom), 
                            str_new_charstr(GENEPRED_SOURCE), 
                            str_new_charstr(GFF_CDS_TYPE), 
                            max(cdsStart, eStart), min(cdsEnd, eEnd), 
                            0, strand, GFF_NULL_FRAME,
                            str_new_charstr(group), TRUE);
        lst_push_ptr(gff->features, f);
        lst_push_ptr(framefeats, f);
      }
    }

    /* set frame */
    if (lst_size(framefeats) > 0) {
      int frame = 0;            
      if (strand == '-') lst_reverse(framefeats);
                                /* framefeats should now be sorted
                                   5'->3' */
      for (i = 0; i < lst_size(framefeats); i++) {
        f = lst_get_ptr(framefeats, i);
        f->frame = frame;
        frame = ((frame + f->end - f->start + 1) % 3);
      }
    }

    lst_free_strings(tmpl1);
    lst_free_strings(tmpl2);
    lst_free_strings(l);
  }

  str_free(line);
  lst_free(l);
  lst_free(tmpl1);
  lst_free(tmpl2);
  lst_free(framefeats);
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
  String *exonStarts = str_new(STR_LONG_LEN), *exonEnds = str_new(STR_LONG_LEN),
    *cdsStarts = str_new(STR_LONG_LEN), *cdsEnds = str_new(STR_LONG_LEN);

  if (feats->groups == NULL)
    die("ERROR: features must be grouped to print as genepred.\n");

  for (i = 0; i < lst_size(feats->groups); i++) {
    int cdsStart = -1, cdsEnd = -1, nexons = 0, ncds_exons = 0;
    char strand = '\0';
    String *seqname = NULL;
    GFF_FeatureGroup *g = lst_get_ptr(feats->groups, i);

    checkInterruptN(i, 100);
    str_clear(exonStarts);
    str_clear(exonEnds);
    str_clear(cdsStarts);
    str_clear(cdsEnds);

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

        str_append_int(cdsStarts, f->start-1);
        str_append_char(cdsStarts, ',');
        str_append_int(cdsEnds, f->end);
        str_append_char(cdsEnds, ',');
        if (nexons == 0 && ncds_exons == 0) {
          strand = f->strand;
          seqname = f->seqname;
        }
        else if (strand != f->strand || !str_equals(seqname, f->seqname))
          die("ERROR (gff_print_genepred): inconsistent strand or seqname in GFF group \"%s\".\n", g->name->chars);
        ncds_exons++;
      }
    }

    if (nexons > 0)
      fprintf(OUTF, "%s\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n", 
              g->name->chars, seqname->chars, strand, g->start-1, g->end, 
              cdsStart, cdsEnd, nexons, exonStarts->chars, exonEnds->chars);

    else if (ncds_exons > 0)    /* this is in case only CDS is specified */
      fprintf(OUTF, "%s\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n", 
              g->name->chars, seqname->chars, strand, g->start-1, g->end, 
              cdsStart, cdsEnd, ncds_exons, cdsStarts->chars, cdsEnds->chars);

  }

  str_free(exonStarts);
  str_free(exonEnds);
  str_free(cdsStarts);
  str_free(cdsEnds);
}
