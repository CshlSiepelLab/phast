/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* stacks - Simple array-based stacks and supporting functions (uses 'lists'). */

/* See documentation for lists2.h.

   $Id: stacks.h,v 1.2 2008-11-12 02:07:59 acs Exp $ */

#ifndef STACKS_H
#define STACKS_H

#include "lists.h"
#include "external_libs.h"

typedef List Stack;

/* Test whether stack is empty.
   Returns 1 if empty, 0 otherwise */
static PHAST_INLINE
int stk_empty(Stack *s) { return lst_empty(s); }

/* Pop object from stack.
   Object is removed.
   Returns address of object at top of stack, or NULL if stack is
   empty.  If object is an int, you will need to access it with
   *((int*)lst_get(...)); if object is a Node* and you seek the
   attribute called "data", you will need to do
   ((Node*)lst_get(...))->data. */
static PHAST_INLINE
void* stk_pop(Stack* s) {
  if (stk_empty(s)) return NULL;
  return lst_arr_get(s, --s->ridx);
}

/* Peek at object at top of stack.
   Object is not removed.
   Returns address of object at top of stack, or NULL if stack is
   empty.  See stk_pop for details. */
static PHAST_INLINE
void* stk_peek(Stack* s) {
  if (stk_empty(s)) return NULL;
  return lst_arr_get(s, s->ridx-1);
}

/* Pop integer from stack.
   Integer is removed.
   Returns integer from top of stack or 0 if stack is empty (use stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when stack is empty. */
static PHAST_INLINE
int stk_pop_int(Stack *s) {
  int *ptr = (int*)stk_pop(s);
  return (ptr == NULL ? 0 : *ptr);
}

/* Pop double from stack.
   Double is removed.
   Returns double from top of stack or 0 if stack is empty (use stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when stack is empty. */
static PHAST_INLINE
double stk_pop_dbl(Stack *s) {
  double *ptr = (double*)stk_pop(s);
  return (ptr == NULL ? 0 : *ptr);
}

/* Pop pointer from stack.
   Pointer is removed.
   Returns pointer from top of stack or NULL if stack is empty. */
static PHAST_INLINE
void* stk_pop_ptr(Stack *s) {
  void **ptr = (void**)stk_pop(s);
  return (ptr == NULL ? NULL : *ptr);
}

/* Peek at integer at top of stack.
   Integer is not removed.
   Returns integer from top of stack or 0 if stack is empty (use stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when stack is empty. */
static PHAST_INLINE
int stk_peek_int(Stack *s) {
  int *ptr = (int*)stk_peek(s);
  return (ptr == NULL ? 0 : *ptr);
}

/* Peek at double at top of stack.
   Double is not removed.
   Returns double from top of stack or 0 if stack is empty (use stk_empty).

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when stack is empty. */
static PHAST_INLINE
double stk_peek_dbl(Stack *s) {
  double *ptr = (double*)stk_peek(s);
  return (ptr == NULL ? 0 : *ptr);
}

/* Peek at pointer at top of stack.
   Pointer is not removed.
   Returns pointer from top of stack or NULL if stack is empty. */
static PHAST_INLINE
void* stk_peek_ptr(Stack *s) {
  void **ptr = (void**)stk_peek(s);
  return (ptr == NULL ? NULL : *ptr);
}

/* Push object on stack */
static PHAST_INLINE
void stk_push(Stack *s, void *o) { lst_push(s, o); }

/* Push integer on stack */
static PHAST_INLINE
void stk_push_int(Stack *s, int i) { lst_push_int(s, i); }

/* Push double on stack */
static PHAST_INLINE
void stk_push_dbl(Stack *s, double d) { lst_push_dbl(s, d); }

/* Push pointer on stack */
static PHAST_INLINE
void stk_push_ptr(Stack *s, void *ptr) { lst_push_ptr(s, ptr); }

/* Obtain size of stack (number of elements).
   Returns number of elements */
static PHAST_INLINE
int stk_size(Stack *s) { return lst_size(s); }

/* Create new stack of objects.
   Returns newly allocated Stack object, with specified starting
   size. */
static PHAST_INLINE
Stack* stk_new(int nelements,	/* Starting number of elements */
	      int elementsz)	/* Size of each element (bytes) */
{ return lst_new(nelements, elementsz); }


/* Create new stack of integers. 
   Like stk_new, but with element size fixed at sizeof(int).
   Returns newly allocated Stack object, with specified starting size. */
static PHAST_INLINE
Stack* stk_new_int(int nelements) /* Starting number of elements */
{ return lst_new(nelements, sizeof(int)); }

/* Create new stack of doubles. 
   Like stk_new, but with element size fixed at sizeof(double).
   Returns newly allocated Stack object, with specified starting size. */
static PHAST_INLINE
Stack* stk_new_dbl(int nelements) /* Starting number of elements */
{ return lst_new(nelements, sizeof(double)); }


/* Create new stack of pointers.
   Like stk_new, but with element size fixed at sizeof(void*).
   Returns newly allocated Stack object, with specified starting size. */
static PHAST_INLINE
Stack* stk_new_ptr(int nelements)
{ return lst_new(nelements, sizeof(void*)); }

/* Clear contents of stack.
   Contents will be cleared but memory will remain allocated. */   
static PHAST_INLINE
void stk_clear(Stack *s) { lst_clear(s); }

/* Free memory associated with stack.
   Stack object itself will be freed also. */
static PHAST_INLINE
void stk_free(Stack *s) { lst_free(s); }

#endif

