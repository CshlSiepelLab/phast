/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file gff.h
    Reading and writing of sequence features in wig format. 
    Obeys file specification at
    http://genome.ucsc.edu/goldenPath/help/wiggle.html
    @ingroup feature
*/


#ifndef WIG_H
#define WIG_H

#include <gff.h>


/** Check if a string is a wig file header and parse the arguments.
    @param line[in] A string which may be a wig header file (fixed or variable step)
    @param fixed[out] If not NULL, will be set to 1 if line is a header for a fixedStep wig, 
      and 0 if line is  header for a variableStep wig.  Will not be set if line is not a wig 
      header.  If not NULL, must be allocated large enough to hold the chrom value.
    @param chrom[out] If not NULL, will be set to the chrom argument in a wig header.  Will 
      not be set if line is not a wig header.
    @param start[out] Will be set to the start argument if line is a fixedStep wig header.
    @param step[out] Will be set to the step argument if line is a fixedStep wig header
    @param span[out] Will be set to the span argument if line is a wig header.  If line is a 
      wig header but span is not given, it will be set to the default value of 1.  Will not 
      be set if line is not a wig header.
    @return 1 if line is a wig header, 0 otherwise.  None 
 */
int wig_parse_header(String *line, int *fixed, char *chrom, int *start, int *step, int *span);

/** Read a wig file and store as a GFF.  Can read either fixed or variable wig files.
  @param F wig file to read
  @return A newly allocated GFF set populated with the data from the wig file.
  @note although not technically allowed in the specification, this function will accept a wig
  file with a mix of fixed/variable sections.
 */
GFF_Set *gff_read_wig(FILE *F);

/** Write a GFF object as a fixedStep wig file.
  @param F file stream to write to
  @param set The GFF set to write.
  @note all elements in set must have the same length and have a defined score.
 */
void wig_print(FILE *outfile, GFF_Set *set);

#endif

