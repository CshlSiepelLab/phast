/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/


/** @file bed.h
   Reading and writing of features to/from .BED files.

   Bed files are read or written as sets of features.  
   These sets of features are stored in GFF_Set objects
   which contain all information pertinant to a each feature.
   
   See http://genome.ucsc.edu/goldenPath/help/customTrack.html 
   @ingroup feature
*/ 

#ifndef BED_H
#define BED_H

#include <gff.h>

/** Create new feature set from file
  @param[in] F File descriptor of BED file
  @param[out] gff Feature set to populate with features from file  
  @note
  @code 
    #include <stdio.h>
    #include <gff.h> //To use GFF_Set object
    #include <bed.h> //To read from BED file

    int main(int argc, char *argv[]) {

      char *filename = "file.bed";
      GFF_Set *FeatureSet = gff_new_set();

      gff_read_from_bed(FeatureSet, phast_fopen(filename, "r"));
    }
  @endcode
*/
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
 @code 
    #include <stdio.h>
    #include <gff.h>  //To use GFF_Set object
    #include <bed.h>  //To read in from BED file
    #include <misc.h> //To use TRUE/FALSE

    int main(int argc, char *argv[]) {

      char *filename = "file.bed";
      GFF_Set *FeatureSet = gff_new_set();

      gff_read_from_bed(FeatureSet, phast_fopen(filename, "r"));
      gff_print_bed(stdout, FeatureSet, FALSE); //If we had groups we would have used TRUE
    }
    OUTPUT:
	chr22	49300442	49300573	hpmrc.1	0	+
	chr22	49301519	49302857	hpmrc.2	0	+
	chr22	49302939	49303037	hpmrc.3	0	+
	chr22	49303476	49304046	hpmrc.4	0	+
	chr22	49309025	49309154	hpmrc.5	0	+
	chr22	49310089	49310131	hpmrc.6	0	+
	chr22	49310179	49310218	hpmrc.7	0	+
	chr22	49314992	49315043	hpmrc.8	0	+
	chr22	49315445	49315520	hpmrc.9	0	+
	chr22	49315555	49315641	hpmrc.10	0	+

  @endcode

*/
void gff_print_bed(FILE *OUTF, GFF_Set *gff, int use_groups);

#endif
