/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/


/** @file bed.h
   Reading and writing of BED files.  See
   http://genome.ucsc.edu/goldenPath/help/customTrack.html 
   @ingroup feature
*/ 

#ifndef BED_H
#define BED_H

#include <gff.h>

/** Create new feature set from file
  @param[in] F File descriptor of BED file
  @param[out] gff Feature set to populate with features from file  */
void gff_read_from_bed(GFF_Set *gff, FILE *F);

/** Save feature set to a file
  @param[in] OUTF File descriptor of output file
  @param[in] gff Feature set to save to file
  @param[in] use_groups If TRUE, all members of a group
                                   are described by a single line in
                                   the BED file (12-column format); if
                                   FALSE, any groups are ignored and
                                   each feature gets its own line.
  @note Features must already be grouped if use_groups == TRUE.
*/
void gff_print_bed(FILE *OUTF, GFF_Set *gff, int use_groups);

#endif
