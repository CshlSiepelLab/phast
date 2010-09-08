/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: multi_msa.h,v 1.2 2008-11-12 02:07:59 acs Exp $ */

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
