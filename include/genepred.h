/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: genepred.h,v 1.4 2008-11-12 02:07:59 acs Exp $ */

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
