/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/
/**
   @file queues.h
   queues - Simple array-based queues and supporting functions (uses 'lists'). 

   Supports storage of objects of arbitrary size in queue format. 
   Convenience functions are available for common data types such as ints,
   doubles, and pointers.  
   @note Inline functions: Note: all functions are inline; there is no queues.c
   \ingroup base
   @see lists.h
*/

#ifndef QUEUES_H
#define QUEUES_H

#include "lists.h"
#include "external_libs.h"

typedef List Queue;

/** \name Queue allocation functions 
 \{ */

/** Create new queue of objects.
   @param nelements Starting number of elements
   @param elementsz Size of each element (bytes)
   @result newly allocated Queue object, with specified starting size.
   @note this is the more generic version of functions like que_new_int or que_new_dbl
   @see que_new_int 
   @see que_new_dbl
  */
static PHAST_INLINE
Queue* que_new(int nelements,   /* Starting number of elements */
              int elementsz)    /* Size of each element (bytes) */
{ return lst_new(nelements, elementsz); }


/** Create new queue of integers. 
   @param nelements Starting number of integers for the queue to hold (will expand dynamically)
   @result newly allocated Queue object, with specified starting number of ints 
   @note Like que_new, but with element size fixed at sizeof(int). 
   @see que_new */
static PHAST_INLINE
Queue* que_new_int(int nelements) /* Starting number of elements */
{ return lst_new(nelements, sizeof(int)); }

/** Create new queue of doubles. 
   @param nelements Starting number of doubles for the queue to hold (will expand dynamically)
   @result newly allocated Queue object, with specified starting number of doubles
   @note Like que_new, but with element size fixed at sizeof(double). 
   @see que_new */
static PHAST_INLINE
Queue* que_new_dbl(int nelements) /* Starting number of elements */
{ return lst_new(nelements, sizeof(double)); }

/** Create new queue of pointers.
   @param nelements Starting number of pointers for the queue to hold (will expand dynamically)
   @result newly allocated Queue object, with specified starting number of pointers
   @note Like que_new, but with element size fixed at sizeof(void*). */
static PHAST_INLINE
Queue* que_new_ptr(int nelements)
{ return lst_new(nelements, sizeof(void*)); }


/** \} \name Queue misc functions 
 \{ */


/** Test whether queue is empty.
   @result 1 if empty, 0 otherwise. */
static PHAST_INLINE
int que_empty(Queue *q) { return lst_empty(q); }

/** Size of queue.
   @result number of elements currently in queue */
static PHAST_INLINE
int que_size(Queue *q) { return lst_size(q); }

/** \} \name Queue push functions 
 \{ */

/** Push pointer to object on queue
   @note same as que_push_ptr
   @see que_push_ptr
 */
static PHAST_INLINE
void que_push(Queue *q, void *o) { lst_push(q, o); }

/** Push integer on queue */
static PHAST_INLINE
void que_push_int(Queue *q, int i) { lst_push_int(q, i); }

/** Push double on queue */
static PHAST_INLINE
void que_push_dbl(Queue *q, double d) { lst_push_dbl(q, d); }

/** Push pointer on queue 
  @note same as que_push
  @see que_push
*/
static PHAST_INLINE
void que_push_ptr(Queue *q, void *ptr) { lst_push_ptr(q, ptr); }

/** \} \name Queue peek functions 
 \{ */


/** Peek at object at beginning of queue.
   @result address of object at beginning of queue, or NULL if stack is empty. 
   @note Object is not removed. 
   @see que_pop */
static PHAST_INLINE
void* que_peek(Queue* q) {
  if (que_empty(q)) return NULL;
  return lst_arr_get(q, q->lidx);
}

/** Peek at integer at beginning of stack.
   @result integer from beginning of queue or 0 if queue is empty (use que_empty).
   @warning return value will be ambiguous when using numeric data
   containing zeroes.  Use que_empty to tell when stack is empty. 
   @note Integer is not removed. 
   @see que_empty */
static PHAST_INLINE
int que_peek_int(Queue *q) {
  int *ptr = (int*)que_peek(q);
  return (ptr == NULL ? 0 : *ptr);
}

/** Peek at double at beginning of queue.
   @result double from beginning of queue or 0 if queue is empty (use que_empty).
   @warning return value will be ambiguous when using numeric data
   containing zeroes.  Use que_empty to tell when queue is empty. 
   @note Double is not removed
   @see que_empty*/
static PHAST_INLINE
double que_peek_dbl(Queue *q) {
  double *ptr = (double*)que_peek(q);
  return (ptr == NULL ? 0 : *ptr);
}

/** Peek at pointer at beginning of queue.
   @result pointer from beginning of queue or NULL if queue is empty. 
   @note pointer is not removed */
static PHAST_INLINE
void* que_peek_ptr(Queue *q) {
  void **ptr = (void**)que_peek(q);
  return (ptr == NULL ? NULL : *ptr);
}

/** \} \name Queue pop functions 
 \{ */


/** Pop object from queue.
   @result address of object at beginning, or NULL if queue is empty. 
   @warning Object is removed. 
   @note If object is an int, you will need to access it with
   *((int*)lst_get(...)); if object is a Node* and you seek the
   attribute called "data", you will need to do
   ((Node*)lst_get(...))->data. */
static PHAST_INLINE
void* que_pop(Queue* q) {
  if (que_empty(q)) return NULL;
  return lst_arr_get(q, q->lidx++);
}

/** Pop integer from queue.
   @result integer from beginning of queue, or 0 if queue is empty (use stk_empty).
   @warning Integer is removed.
   @warning return value will be ambiguous when using numeric data
   containing zeroes.  Use que_empty to tell when queue is empty. 
   @see que_empty */
static PHAST_INLINE
int que_pop_int(Queue *q) {
  int *ptr = (int*)que_pop(q);
  return (ptr == NULL ? 0 : *ptr);
}

/** Pop double from queue.
   @result double from beginning of queue or 0 if queue is empty (use que_empty).
   @warning Double is removed
   @warning return value will be ambiguous when using numeric data
   containing zeroes.  Use que_empty to tell when stack is empty. 
   @see que_empty */
static PHAST_INLINE
double que_pop_dbl(Queue *q) {
  double *ptr = (double*)que_pop(q);
  return (ptr == NULL ? 0 : *ptr);
}

/** Pop pointer from queue.
   @result pointer from beginning of queue or NULL if queue is empty. 
   @warning Pointer is removed */
static PHAST_INLINE
void* que_pop_ptr(Queue *q) {
  void **ptr = (void**)que_pop(q);
  return (ptr == NULL ? NULL : *ptr);
}

/** \} \name Queue cleanup functions 
 \{ */

/** Clear contents of queue.
   @warning Contents will be cleared but memory will remain allocated. To free elements in the queue use que_free
   @see que_free
 */   
static PHAST_INLINE
void que_clear(Queue *q) { lst_clear(q); }

/** Free memory associated with queue.
   @warning pointers in the queue will be cleared, but not the contents they point to.
   @note Queue object itself will be freed also. */
static PHAST_INLINE
void que_free(Queue *q) { lst_free(q); }

/** \} */
#endif

