/* queues - Simple array-based queues and supporting functions (uses 'lists'). */

/* See documentation for lists2.h.

   $Id: queues.h,v 1.1.1.1 2004-06-03 22:43:11 acs Exp $
   Written by Adam Siepel, Spring 2001, Summer 2002
   Copyright 2001, 2002, Adam Siepel, University of California
*/

#ifndef QUEUES_H
#define QUEUES_H

#include "lists.h"

typedef List Queue;

/* Test whether queue is empty.
   Returns 1 if empty, 0 otherwise */
int que_empty(Queue *q);

/* Pop object from queue.
   Object is removed.
   Returns address of object at beginning of queue, or NULL if queue
   is empty.  If object is an int, you will need to access it with
   *((int*)lst_get(...)); if object is a Node* and you seek the
   attribute called "data", you will need to do
   ((Node*)lst_get(...))->data. */
void* que_pop(Queue* q); 

/* Peek at object at beginning of queue.
   Object is not removed.
   Returns address of object at beginning of queue, or NULL if queue
   is empty. See que_pop for details.*/
void* que_peek(Queue* q);

/* Pop integer from queue.
   Integer is removed.
   Returns integer from beginning of queue or 0 if queue is empty (use
   stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when queue is empty. */
int que_pop_int(Queue *q);

/* Pop double from queue.
   Double is removed.
   Returns double from beginning of queue or 0 if queue is empty (use
   stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when queue is empty. */
double que_pop_dbl(Queue *q);

/* Pop pointer from queue.
   Pointer is removed.
   Returns pointer from beginning of queue or NULL if queue is empty. */
void* que_pop_ptr(Queue *q);

/* Peek at integer at beginning of queue.
   Integer is not removed.
   Returns integer from beginning of queue or 0 if queue is empty (use
   stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when queue is empty. */
int que_peek_int(Queue *q);

/* Peek at double at beginning of queue.
   Double is not removed.
   Returns double from beginning of queue or 0 if queue is empty (use
   stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when queue is empty. */
double que_peek_dbl(Queue *q);

/* Peek at pointer at beginning of queue.
   Pointer is not removed.
   Returns pointer from beginning of queue or NULL if queue is empty. */
void* que_pop_ptr(Queue *q);

/* Push object on queue */
void que_push(Queue *q, void *o);

/* Push integer on queue */
void que_push_int(Queue *q, int i);

/* Push double on queue */
void que_push_dbl(Queue *q, double d);

/* Push pointer on queue */
void que_push_ptr(Queue *q, void *ptr);

/* Size of queue.
   Returns number of elements */
int que_size(Queue *q);

/* Create new queue of objects.
   Returns newly allocated Queue object, with
   specified starting size. */
Queue* que_new(int nelements,	/* Starting number of elements */
	      int elementsz);	/* Size of each element (bytes) */

/* Create new queue of integers. 
   Like que_new, but with element size fixed at sizeof(int).
   Returns newly allocated Queue object, with specified starting size. */
Queue* que_new_int(int nelements); /* Starting number of elements */

/* Create new queue of doubles. 
   Like que_new, but with element size fixed at sizeof(double).
   Returns newly allocated Queue object, with specified starting size. */
Queue* que_new_dbl(int nelements); /* Starting number of elements */

/* Create new queue of pointers.
   Like que_new, but with element size fixed at sizeof(void*).
   Returns newly allocated Queue object, with specified starting size. */
Queue* que_new_ptr(int nelements);

/* Clear contents of queue.
   Contents will be cleared but memory will remain allocated. */   
void que_clear(Queue *q);

/* Free memory associated with queue.
   Queue object itself will be freed also. */
void que_free(Queue *q);

/***************************************************************************/
/* Inline functions; note: these are redefined in queues.c for             */
/* situations in which inlining is not available                           */
/***************************************************************************/

extern inline
int que_empty(Queue *q) { return lst_empty(q); }

extern inline
void* que_pop(Queue* q) {
  if (que_empty(q)) return NULL;
  return lst_arr_get(q, q->lidx++);
}

extern inline
void* que_peek(Queue* q) {
  if (que_empty(q)) return NULL;
  return lst_arr_get(q, q->lidx);
}

extern inline
int que_pop_int(Queue *q) {
  int *ptr = (int*)que_pop(q);
  return (ptr == NULL ? 0 : *ptr);
}

extern inline
double que_pop_dbl(Queue *q) {
  double *ptr = (double*)que_pop(q);
  return (ptr == NULL ? 0 : *ptr);
}

extern inline
void* que_pop_ptr(Queue *q) {
  void **ptr = (void**)que_pop(q);
  return (ptr == NULL ? NULL : *ptr);
}

extern inline
int que_peek_int(Queue *q) {
  int *ptr = (int*)que_peek(q);
  return (ptr == NULL ? 0 : *ptr);
}

extern inline
double que_peek_dbl(Queue *q) {
  double *ptr = (double*)que_peek(q);
  return (ptr == NULL ? 0 : *ptr);
}

extern inline
void* que_peek_ptr(Queue *q) {
  void **ptr = (void**)que_peek(q);
  return (ptr == NULL ? NULL : *ptr);
}

extern inline
void que_push(Queue *q, void *o) { lst_push(q, o); }

extern inline
void que_push_int(Queue *q, int i) { lst_push_int(q, i); }

extern inline
void que_push_dbl(Queue *q, double d) { lst_push_dbl(q, d); }

extern inline
void que_push_ptr(Queue *q, void *ptr) { lst_push_ptr(q, ptr); }

extern inline
int que_size(Queue *q) { return lst_size(q); }

extern inline
Queue* que_new(int nelements, int elementsz)
{ return lst_new(nelements, elementsz); }

extern inline
Queue* que_new_int(int nelements) 
{ return lst_new(nelements, sizeof(int)); }

extern inline
Queue* que_new_dbl(int nelements) 
{ return lst_new(nelements, sizeof(double)); }

extern inline
Queue* que_new_ptr(int nelements) 
{ return lst_new(nelements, sizeof(void*)); }

extern inline
void que_clear(Queue *q) { lst_clear(q); }

extern inline
void que_free(Queue *q) { lst_free(q); }

#endif

