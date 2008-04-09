/* $Id: pssm.h,v 1.1 2008-04-09 13:01:58 acs Exp $ */

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
