/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file multi_msa.h
   Read more than one 'Multiple Sequence Alignment' into a single object
   @ingroup msa
*/

#ifndef MULTI_MSA_H
#define MULTI_MSA_H

#include <msa.h>
#include <indel_history.h>

/** Structure holding multiple individual MSA objects. */
typedef struct {
  int nblocks;                  /**< Number of individual MSA objects */
  MSA **blocks;                 /**< The actual array of MSA objects */
  List *seqnames;               /**< The names of the sequence blocks (file
				  names from file list) */
  IndelHistory **ih;            /**< Array of indel histories (may be null) */
} Multi_MSA;

/** Creates a new multi-fasta object containing alignments read from files
    specified in a config file.
    Config file contains the number of alignments,
    format, alphabet and list of files to be read. If "alphabet" is
    NULL, default alphabet for DNA will be used.  This routine will
    abort if the sequence contains a character not in the alphabet.
   @param F File descriptor of file containing multi MSA
   @param do_ih Whether to create indel history
  */
Multi_MSA *msa_multimsa_new(FILE *F, int do_ih);

#endif
