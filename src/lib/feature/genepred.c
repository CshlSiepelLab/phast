/* $Id: genepred.c,v 1.1 2004-06-21 19:50:33 acs Exp $
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

#define GENEPRED_SOURCE "genepred"
#define GENEPRED_EXON "exon"
#define GENEPRED_CDS "CDS"
#define GENEPRED_UTR5 "5'UTR"
#define GENEPRED_UTR3 "3'UTR"

/** Fill out a GFF_Set from a genepred file. */
void gff_read_from_genepred(GFF_Set *gff, FILE *F, int do_utr) {
  String *line = str_new(STR_LONG_LEN);
  List *l = lst_new_ptr(12), *tmpl1 = lst_new_ptr(10), 
    *tmpl2 = lst_new_ptr(10);
  int i, lineno = 0, frame = 0;

  while (str_readline(line, F) != EOF) {
    int txStart = 0, txEnd = 0, cdsStart = 0, cdsEnd = 0, 
      exonCount = 0;
    String *name, *chrom, *tmpstr;
    char group[STR_MED_LEN];
    char strand;

    lineno++;

    str_trim(line);
    if (line->length == 0) continue;

    str_split(line, "\t", l);

    if (lst_size(l) != 10)
      die("ERROR (line %d): 10 columns required in genepred file.\n", lineno);

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

    if (str_as_int(lst_get_ptr(l, 7), &exonCount) != 0)
      die("ERROR (line %d): can't parse exonCount in genepred file.\n", 
          lineno);

    str_split(lst_get_ptr(l, 8), ",", tmpl1);
    str_split(lst_get_ptr(l, 9), ",", tmpl2);
    
    if (exonCount != lst_size(tmpl1) || lst_size (tmpl1) != lst_size(tmpl2))
      die("ERROR (line %d): exonStarts or exonEnds don't match exonCount in genepred file.\n", lineno);

    for (i = 0; i < exonCount; i++) {
      int eStart = 0, eEnd = 0;

      if (str_as_int(lst_get_ptr(tmpl1, i), &eStart) != 0 ||
              str_as_int(lst_get_ptr(tmpl2, i), &eEnd) != 0)
        die("ERROR (line %d): can't parse exonStarts or exonEnds in genepred file.\n", lineno);

      eStart++;
      sprintf(group, "transcript_id \"%s\"", name->chars);

      /* create one feature for the whole exon and another for the CDS
         portion (if necessary) */
      lst_push_ptr(gff->features, 
                   gff_new_feature(str_dup(chrom), 
                                   str_new_charstr(GENEPRED_SOURCE), 
                                   str_new_charstr(GENEPRED_EXON), 
                                   eStart, eEnd, 0, strand, GFF_NULL_FRAME,
                                   str_new_charstr(group), TRUE));

      if ((eStart >= cdsStart && eStart <= cdsEnd) ||
          (eEnd >= cdsStart && eEnd <= cdsEnd))
        lst_push_ptr(gff->features, 
                     gff_new_feature(str_dup(chrom), 
                                     str_new_charstr(GENEPRED_SOURCE), 
                                     str_new_charstr(GENEPRED_CDS), 
                                     max(cdsStart, eStart), min(cdsEnd, eEnd), 
                                     0, strand, frame,
                                     str_new_charstr(group), TRUE));

      if (do_utr && eStart < cdsStart)
        lst_push_ptr(gff->features, 
                     gff_new_feature(str_dup(chrom), 
                                     str_new_charstr(GENEPRED_SOURCE), 
                                     strand == '+' ? 
                                     str_new_charstr(GENEPRED_UTR5) : 
                                     str_new_charstr(GENEPRED_UTR3), 
                                     eStart, min(cdsStart-1, eEnd), 
                                     0, strand, GFF_NULL_FRAME,
                                     str_new_charstr(group), TRUE));

      if (do_utr && eEnd > cdsEnd)
        lst_push_ptr(gff->features, 
                     gff_new_feature(str_dup(chrom), 
                                     str_new_charstr(GENEPRED_SOURCE), 
                                     strand == '+' ? 
                                     str_new_charstr(GENEPRED_UTR3) :
                                     str_new_charstr(GENEPRED_UTR5), 
                                     max(cdsEnd+1, eStart), eEnd,
                                     0, strand, GFF_NULL_FRAME,
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
}
