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

#include <phast/heap.h>

HeapNode* hp_new_node(double val, void *auxdata) {
  HeapNode *node = malloc(sizeof(HeapNode));
  node->val = val;
  node->auxdata = auxdata;
  node->child = NULL;
  node->sibling = NULL;
  return node;
}

HeapNode* hp_meld(HeapNode *a, HeapNode *b) {
  if (a == NULL)
    return b;
  if (b == NULL)
    return a;
  if (a->val <= b->val) {
    b->sibling = a->child;
    a->child = b;
    return a;
  }
  else {
    a->sibling = b->child;
    b->child = a;
    return b;
  }
}

HeapNode* hp_insert(HeapNode *heap, double val, void *auxdata) {
  HeapNode *node = hp_new_node(val, auxdata);
  return hp_meld(heap, node);
}

HeapNode* hp_meld_two_pass(HeapNode *node) {
  if (node == NULL || node->sibling == NULL)
    return node;

  HeapNode *a = node;
  HeapNode *b = node->sibling;
  HeapNode *rest = b->sibling;

  a->sibling = NULL;
  b->sibling = NULL;
  return hp_meld(hp_meld(a, b), hp_meld_two_pass(rest));
}

HeapNode* hp_delete_min(HeapNode *heap, void **min_auxdata) {
  if (heap == NULL) return NULL;
  *min_auxdata = heap->auxdata;
  HeapNode *new_heap = hp_meld_two_pass(heap->child);
  free(heap);
  return new_heap;
}
