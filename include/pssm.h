/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: pssm.h,v 1.2 2008-11-12 02:07:59 acs Exp $ */

#ifndef PSSM_H
#define PSSM_H

#include <vector.h>

/* simple motif struct.  Eventually could be extended with background
   model, logodds scores, higher-order models, etc. */
typedef struct {
  int width;
  int alphsize;
  char *alphabet;
  Vector **probs;
} PSSM;


PSSM *mot_new(int width, char *alphabet);
PSSM *mot_read(FILE *F);
void mot_write(FILE *F, PSSM *m);
void mot_free(PSSM *m);

#endif
