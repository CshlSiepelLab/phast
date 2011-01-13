/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file genepred.h
    Reading and writing of genepred files, used for UCSC Genome Browser.
    \ingroup feature
*/

#ifndef GENEPRED_H
#define GENEPRED_H

#include <gff.h>

/** Fill out a GFF_Set from a genepred file.
  @param[in] F Gff file to read from
  @param[out] gff Gff object populated with data from file F
 */
void gff_read_from_genepred(GFF_Set *gff, FILE *F);

/** Write a GFF_Set in genepred format.
    @pre Features must already be grouped as desired.
    @param OUTF Output Stream
    @param feats Set to write to output stream
    @note Features will be sorted within groups as a side effect.
*/
void gff_print_genepred(FILE *OUTF, GFF_Set *feats);

#endif
