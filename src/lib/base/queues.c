/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* Simple array-based queues and supporting functions. 

   See documentation for lists.h.

   $Id: queues.c,v 1.2 2008-11-12 02:07:59 acs Exp $
*/

#include <stdlib.h>
#include "queues.h"

int que_empty(Queue *q) { return lst_empty(q); }

void* que_pop(Queue* q) {
  if (que_empty(q)) return NULL;
  return lst_arr_get(q, q->lidx++);
}

void* que_peek(Queue* q) {
  if (que_empty(q)) return NULL;
  return lst_arr_get(q, q->lidx);
}

int que_pop_int(Queue *q) {
  int *ptr = (int*)que_pop(q);
  return (ptr == NULL ? 0 : *ptr);
}

double que_pop_dbl(Queue *q) {
  double *ptr = (double*)que_pop(q);
  return (ptr == NULL ? 0 : *ptr);
}

void* que_pop_ptr(Queue *q) {
  void **ptr = (void**)que_pop(q);
  return (ptr == NULL ? NULL : *ptr);
}

int que_peek_int(Queue *q) {
  int *ptr = (int*)que_peek(q);
  return (ptr == NULL ? 0 : *ptr);
}

double que_peek_dbl(Queue *q) {
  double *ptr = (double*)que_peek(q);
  return (ptr == NULL ? 0 : *ptr);
}

void* que_peek_ptr(Queue *q) {
  void **ptr = (void**)que_peek(q);
  return (ptr == NULL ? NULL : *ptr);
}

void que_push(Queue *q, void *o) { lst_push(q, o); }

void que_push_int(Queue *q, int i) { lst_push_int(q, i); }

void que_push_dbl(Queue *q, double d) { lst_push_dbl(q, d); }

void que_push_ptr(Queue *q, void *ptr) { lst_push_ptr(q, ptr); }

int que_size(Queue *q) { return lst_size(q); }

Queue* que_new(int nelements, int elementsz)
{ return lst_new(nelements, elementsz); }

Queue* que_new_int(int nelements) 
{ return lst_new(nelements, sizeof(int)); }

Queue* que_new_dbl(int nelements) 
{ return lst_new(nelements, sizeof(double)); }

Queue* que_new_ptr(int nelements) 
{ return lst_new(nelements, sizeof(void*)); }

void que_clear(Queue *q) { lst_clear(q); }

void que_free(Queue *q) { lst_free(q); }
