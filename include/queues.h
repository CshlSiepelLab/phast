/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* queues - Simple array-based queues and supporting functions (uses 'lists'). */

/* See documentation for lists2.h.

   $Id: queues.h,v 1.2 2008-11-12 02:07:59 acs Exp $ */

#ifndef QUEUES_H
#define QUEUES_H

#include "lists.h"
#include "external_libs.h"

typedef List Queue;

/** Inline functions: Note: all functions are inline; 
    there is no queues.c **/

/* Test whether queue is empty.
   Returns 1 if empty, 0 otherwise */
static PHAST_INLINE
int que_empty(Queue *q) { return lst_empty(q); }

/* Pop object from queue.
   Object is removed.
   Returns address of object at beginning of queue, or NULL if queue
   is empty.  If object is an int, you will need to access it with
   *((int*)lst_get(...)); if object is a Node* and you seek the
   attribute called "data", you will need to do
   ((Node*)lst_get(...))->data. */
static PHAST_INLINE
void* que_pop(Queue* q) {
  if (que_empty(q)) return NULL;
  return lst_arr_get(q, q->lidx++);
}

/* Peek at object at beginning of queue.
   Object is not removed.
   Returns address of object at beginning of queue, or NULL if queue
   is empty. See que_pop for details.*/
static PHAST_INLINE
void* que_peek(Queue* q) {
  if (que_empty(q)) return NULL;
  return lst_arr_get(q, q->lidx);
}

/* Pop integer from queue.
   Integer is removed.
   Returns integer from beginning of queue or 0 if queue is empty (use
   stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when queue is empty. */
static PHAST_INLINE
int que_pop_int(Queue *q) {
  int *ptr = (int*)que_pop(q);
  return (ptr == NULL ? 0 : *ptr);
}

/* Pop double from queue.
   Double is removed.
   Returns double from beginning of queue or 0 if queue is empty (use
   stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when queue is empty. */
static PHAST_INLINE
double que_pop_dbl(Queue *q) {
  double *ptr = (double*)que_pop(q);
  return (ptr == NULL ? 0 : *ptr);
}

/* Pop pointer from queue.
   Pointer is removed.
   Returns pointer from beginning of queue or NULL if queue is empty. */
static PHAST_INLINE
void* que_pop_ptr(Queue *q) {
  void **ptr = (void**)que_pop(q);
  return (ptr == NULL ? NULL : *ptr);
}

/* Peek at integer at beginning of queue.
   Integer is not removed.
   Returns integer from beginning of queue or 0 if queue is empty (use
   stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when queue is empty. */
static PHAST_INLINE
int que_peek_int(Queue *q) {
  int *ptr = (int*)que_peek(q);
  return (ptr == NULL ? 0 : *ptr);
}

/* Peek at double at beginning of queue.
   Double is not removed.
   Returns double from beginning of queue or 0 if queue is empty (use
   stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when queue is empty. */
static PHAST_INLINE
double que_peek_dbl(Queue *q) {
  double *ptr = (double*)que_peek(q);
  return (ptr == NULL ? 0 : *ptr);
}

/* Peek at pointer at beginning of queue.
   Pointer is not removed.
   Returns pointer from beginning of queue or NULL if queue is empty. */
static PHAST_INLINE
void* que_peek_ptr(Queue *q) {
  void **ptr = (void**)que_peek(q);
  return (ptr == NULL ? NULL : *ptr);
}

/* Push object on queue */
static PHAST_INLINE
void que_push(Queue *q, void *o) { lst_push(q, o); }

/* Push integer on queue */
static PHAST_INLINE
void que_push_int(Queue *q, int i) { lst_push_int(q, i); }

/* Push double on queue */
static PHAST_INLINE
void que_push_dbl(Queue *q, double d) { lst_push_dbl(q, d); }

/* Push pointer on queue */
static PHAST_INLINE
void que_push_ptr(Queue *q, void *ptr) { lst_push_ptr(q, ptr); }

/* Size of queue.
   Returns number of elements */
static PHAST_INLINE
int que_size(Queue *q) { return lst_size(q); }

/* Create new queue of objects.
   Returns newly allocated Queue object, with
   specified starting size. */
static PHAST_INLINE
Queue* que_new(int nelements,	/* Starting number of elements */
	      int elementsz)	/* Size of each element (bytes) */
{ return lst_new(nelements, elementsz); }


/* Create new queue of integers. 
   Like que_new, but with element size fixed at sizeof(int).
   Returns newly allocated Queue object, with specified starting size. */
static PHAST_INLINE
Queue* que_new_int(int nelements) /* Starting number of elements */
{ return lst_new(nelements, sizeof(int)); }

/* Create new queue of doubles. 
   Like que_new, but with element size fixed at sizeof(double).
   Returns newly allocated Queue object, with specified starting size. */
static PHAST_INLINE
Queue* que_new_dbl(int nelements) /* Starting number of elements */
{ return lst_new(nelements, sizeof(double)); }

/* Create new queue of pointers.
   Like que_new, but with element size fixed at sizeof(void*).
   Returns newly allocated Queue object, with specified starting size. */
static PHAST_INLINE
Queue* que_new_ptr(int nelements)
{ return lst_new(nelements, sizeof(void*)); }

/* Clear contents of queue.
   Contents will be cleared but memory will remain allocated. */   
static PHAST_INLINE
void que_clear(Queue *q) { lst_clear(q); }

/* Free memory associated with queue.
   Queue object itself will be freed also. */
static PHAST_INLINE
void que_free(Queue *q) { lst_free(q); }

#endif

