/* calculation of Robinson Foulds distances */

#ifndef RF_H
#define RF_H

#include <stdio.h>

/* bitset for up to many thousands of leaves */
typedef struct {
  int W;            /* number of 64-bit words */
  uint64_t *w;      /* words */
} BitMask;

/* dynamic array of BitMask* */
typedef struct {
  BitMask **a; int size, cap;
} MaskVec;

double tr_robinson_foulds(TreeNode *t1, TreeNode *t2);

#endif
