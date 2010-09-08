/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: bed.h,v 1.5 2008-11-12 02:07:59 acs Exp $ */

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
