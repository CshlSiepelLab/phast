/* $Id: genepred.h,v 1.3 2004-06-22 22:14:06 acs Exp $
   Written by Adam Siepel, 2004
   Copyright 2004, Adam Siepel, University of California */

/** \file genepred.h
    Reading and writing of genepred files, used for UCSC Genome Browser.
    \ingroup feature
*/

#ifndef GENEPRED_H
#define GENEPRED_H

#include <gff.h>

void gff_read_from_genepred(GFF_Set *gff, FILE *F);

void gff_print_genepred(FILE *OUTF, GFF_Set *feats);

#endif
