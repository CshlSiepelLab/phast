/* $Id: mutli_msa.h,v 0.1 2008/08/09
   Written by Adam Diehl, 2008
   Copyright 2008, Adam Diehl, Cornell University */

/** \file multi_msa.h
   Multiple sequence alignments.
   \ingroup msa
*/

#ifndef MULTI_MSA_H
#define MULTI_MSA_H

#include <msa.h>
#include <indel_history.h>

/** Structure holding multiple individual MSA objects. **/
typedef struct {
  int nblocks;                  /*< Number of individual MSA objects */
  MSA **blocks;                 /*< The actual array of MSA objects */
  List *seqnames;               /*< The names of the sequence blocks (file
				  names from file list) */
  IndelHistory **ih;            /*< Array of indel histories (may be null) */
} Multi_MSA;

Multi_MSA *msa_multimsa_new(FILE *F, int do_ih);

#endif
