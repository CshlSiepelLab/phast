/* $Id: bed.h,v 1.4 2004-06-23 06:03:09 acs Exp $
   Written by Adam Siepel, 2004
   Copyright 2004, Adam Siepel, University of California */

/** \file bed.h
   Reading and writing of BED files.  See
   http://genome.ucsc.edu/goldenPath/help/customTrack.html 
   \ingroup feature
*/ 

#ifndef BED_H
#define BED_H

#include <gff.h>

void gff_read_from_bed(GFF_Set *gff, FILE *F);

void gff_print_bed(FILE *OUTF, GFF_Set *gff, int use_groups);

#endif
