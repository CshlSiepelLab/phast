/* $Id: bed.h,v 1.2 2004-06-14 03:06:21 acs Exp $
   Written by Adam Siepel, 2004
   Copyright 2004, Adam Siepel, University of California */

/* Functions for reading and writing BED files, for use with the UCSC
   genome browser. */

#ifndef BED_H
#define BED_H

#include <gff.h>

void gff_read_from_bed(GFF_Set *gff, FILE *F);

void gff_print_bed(FILE *OUTF, GFF_Set *gff, char *groupby, List *include);

#endif
