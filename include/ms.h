/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/


/** @file ms.h
   Multiple Sequences object

   Object holds multiple sequences (like an MSA) but are not necessairly aligned
   @ingroup motif
*/ 


#ifndef MS_H
#define MS_H
#include <lists.h>

typedef struct {
  double rangeLow;
  double rangeHigh;
  int nseqs;                    /**< Number of sequences */
  char *alphabet;               /**< Alphabet (see #DEFAULT_ALPHABET) */
  int inv_alphabet[NCHARS];     /**< Inverse of 'alphabet' (maps characters to their index in 'alphabet') */
  char **names;			/**< Names of each sequence (2D char array) */
  char **seqs;			/**< Sequence data (2D char array) */
  int *idx_offsets;		/**< Index offset of each sequence*/
  char *missing;                /**< Recognized missing data characters */
  int is_missing[NCHARS];       /**< Fast lookup of whether character
                                   is missing data char */
} MS;

#endif
