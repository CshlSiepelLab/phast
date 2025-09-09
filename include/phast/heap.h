/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* Simple min-heap implementation.  Implemented as a pairing heap,
   with lazy updates for high amortized efficiency (insert,
   decrease-key, delete-min all O(log n) or better). */

#ifndef HEAP_H
#define HEAP_H

#include <stdio.h>
#include <stdlib.h>

typedef struct HeapNode {
  double val;
  void *auxdata;
  struct HeapNode *child;
  struct HeapNode *sibling;
} HeapNode;

HeapNode* hp_new_node(double val, void *auxdata);
HeapNode* hp_meld(HeapNode *a, HeapNode *b);
HeapNode* hp_insert(HeapNode *heap, double val, void *auxdata);
HeapNode* hp_meld_two_pass(HeapNode *node);
HeapNode* hp_delete_min(HeapNode *heap, void **min_auxdata);
void hp_free(HeapNode *heap);
void hp_dump(HeapNode *heap, FILE *F);

#endif
