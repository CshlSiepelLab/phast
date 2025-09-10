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
#include <phast/stacks.h>

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
    b->sibling = a->child; /* put b in front of a's children */
    a->child = b;
    return a;
  }
  else {
    a->sibling = b->child; /* put a in front of b's children */
    b->child = a;
    return b;
  }
}

HeapNode* hp_insert(HeapNode *heap, double val, void *auxdata) {
  HeapNode *node = hp_new_node(val, auxdata);
  return hp_meld(heap, node);
}

HeapNode* hp_meld_two_pass_recur(HeapNode *node) {
  if (node == NULL || node->sibling == NULL)
    return node;

  HeapNode *a = node;
  HeapNode *b = node->sibling;
  HeapNode *rest = b->sibling;

  a->sibling = NULL;
  b->sibling = NULL;
  return hp_meld(hp_meld(a, b), hp_meld_two_pass(rest));
  /* FIXME: better to do iteratively */
}

/* iterative version of two-pass meld, needed for larger heaps */
/* Phase 1: pairwise meld siblings left->right, pushing results on a stack (via sibling).
   Phase 2: pop stack right->left, melding into a single root. */
HeapNode* hp_meld_two_pass(HeapNode *first) {
  if (first == NULL || first->sibling == NULL)
    return first;

  /* ----- pass 1: pair neighbors and push merged pairs onto stack ----- */
  HeapNode *stack = NULL;
  HeapNode *p = first;

  while (p != NULL) {
    HeapNode *a = p;
    HeapNode *b = a->sibling;
    HeapNode *next = (b != NULL ? b->sibling : NULL);

    /* detach a and (maybe) b from the sibling chain */
    a->sibling = NULL;
    if (b != NULL) b->sibling = NULL;

    /* pairwise meld a and b (if b exists) */
    if (b != NULL) a = hp_meld(a, b);

    /* push result onto stack (reuse sibling as stack next) */
    a->sibling = stack;
    stack = a;

    p = next;
  }

  /* ----- pass 2: meld stack right->left into a single heap ----- */
  HeapNode *acc = NULL;
  while (stack != NULL) {
    HeapNode *t = stack;
    stack = stack->sibling;
    t->sibling = NULL;
    acc = (acc == NULL) ? t : hp_meld(acc, t);
  }

  return acc;
}

HeapNode* hp_delete_min(HeapNode *heap, void **min_auxdata) {
  if (heap == NULL) return NULL;
  if (min_auxdata != NULL) *min_auxdata = heap->auxdata;
  HeapNode *new_heap = hp_meld_two_pass(heap->child);
  free(heap);
  return new_heap;
}

/* warning: does not free auxdata */
void hp_free(HeapNode *heap) {
  if (heap == NULL) 
    return;

  Stack *st = stk_new_ptr(128);   /* initial cap; grows as needed */
  stk_push_ptr(st, heap);

  while (!stk_empty(st)) {
    HeapNode *node = (HeapNode*)stk_pop_ptr(st);
    if (node == NULL) continue;
    if (node->sibling != NULL) stk_push_ptr(st, node->sibling);
    if (node->child != NULL) stk_push_ptr(st, node->child);
    free(node);
  }

  stk_free(st);
}

void hp_dump(HeapNode *heap, FILE *F) {
  if (heap == NULL)
    return;
  
  fprintf(F, "Heap Node:\n\tval = %f\n", heap->val);

  if (heap->sibling == NULL)
    fprintf(F, "\tsibling = NULL\n");
  else {
    fprintf(F, "\tsibling = (%f)\n", heap->sibling->val);
    hp_dump(heap->sibling, F);
  }

  if (heap->child == NULL)
    fprintf(F, "\tchild = NULL\n");
  else {
    fprintf(F, "\tchild = (%f)\n", heap->child->val);
    hp_dump(heap->child, F);
  }
}
