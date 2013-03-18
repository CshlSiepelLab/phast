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
  int num_msa;                  /**< Number of individual MSA objects */
  MSA **msa;                 /**< The actual array of MSA objects */
  List *seqnames;               /**< The names of the msas (file
				  names from file list) */
  MSA *concat_msa;              /**< If not NULL, contains all the MSAs in 
				     concatenated together */
  int *start_coord;             /**< If not NULL, contains start coordinate
				     (1-based) of each msa from the frame
				     of reference of the 1st species */
  int *end_coord;               /**< If not NULL, contains end coordinate
				     (1-based, inclusive) of each msa from
				     the frame of reference of the 1st 
				     species */
  String **type;                /**< If not NULL, contains a string identifying
				     the type of alignment (ie, from the 
				     feature field of a gff if multi-msa created
				     from gff */
} Multi_MSA;

/** Creates a new multi-msa object containing alignments read from files
    specified in a config file.
    Config file contains the number of alignments,
    format, alphabet and list of files to be read. If "alphabet" is
    NULL, default alphabet for DNA will be used.  This routine will
    abort if the sequence contains a character not in the alphabet.
   @param F File descriptor of file containing multi MSA
  */
Multi_MSA *multimsa_new_from_files(FILE *F);

/** Creates a new multi-msa object containing all alignments passed in.
    @param msas An array of MSA objects.  These msas will NOT be copied, rather the Multi_MSA object will contain pointers to them.
    @param nmsa The length of msas.
    @return A Multi_MSA object containing pointers to all the msas.  Note that the original msas are not copied!
 */
Multi_MSA *multimsa_new_from_msas(MSA **msas, int nmsa);

/** Compute concatenation of MSAs in Multi_MSA object
    @param mmsa A multi_msa object.  If mmsa->concat_msa is not-NULL, assume it is correct and return without doing anything.  Otherwise, set mmsa->concat_msa to the concatenation of all msas in the Multi_MSA (concat_msa will be newly allocated).
 */
void multimsa_make_concat(Multi_MSA *mmsa);

/** Add a MSA to a Multi_MSA object.
    Does not copy the alignment.  If mmsa->concat_msa is not NULL, the new alignment will also be added to the concatenated alignment.
    @param mmsa A Multi_MSA object
    @param msa The alignment to add to mmsa.
 */
void multimsa_add_msa(Multi_MSA *mmsa, MSA *msa);


#endif
