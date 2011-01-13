/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file stacks.h
   Simple array-based stacks and supporting functions.

   Supports storage of objects of arbitrary size in a stack format.  
   Convenience functions are available available for common data 
   types such as ints, doubles, and pointers. 
   @ingroup base
   @see lists.h
*/

#ifndef STACKS_H
#define STACKS_H

#include "lists.h"
#include "external_libs.h"

/** Basic Stack object. 
    @note Stack object has same format as list objects.
     The only difference between the two is that a stack object
     is required to use stack functions.
    @see list.h
*/
typedef List Stack;


/** \name Stack allocation functions 
 \{ */

/** Create new stack with specified starting size.
   @param nelements Initial capacity of new stack
   @param elementz Size of each element (bytes)
   @result newly allocated stack
   @note stack will grow as needed, but is most efficient when nelements=expected stack size
*/
static PHAST_INLINE
Stack* stk_new(int nelements,	/* Starting number of elements */
	      int elementsz)	/* Size of each element (bytes) */
{ return lst_new(nelements, elementsz); }


/** Create new stack of integers. 
   @param nelements Initial capacity of stack of integers the stack
   @result newly allocated Stack object that stores integers. 
   @note stack will grow as needed, but is most efficient when nelements=expected stack size
   @see stk_new*/
static PHAST_INLINE
Stack* stk_new_int(int nelements) /* Starting number of elements */
{ return lst_new(nelements, sizeof(int)); }

/** Create new stack of doubles. 
   @param nelements Initial capacity of new stack
   @result newly allocated Stack object that stores doubles. 
   @note stack will grow as needed, but is most efficient when nelements=expected stack size
   @see stk_new*/
static PHAST_INLINE
Stack* stk_new_dbl(int nelements) /* Starting number of elements */
{ return lst_new(nelements, sizeof(double)); }


/** Create new stack of pointers.
   @param nelements Initial capacity of new stack 
   @result newly allocated Stack object, that stores pointers 
   @note stack will grow as needed, but is most efficient when nelements=expected stack size
   @see stk_new*/
static PHAST_INLINE
Stack* stk_new_ptr(int nelements)
{ return lst_new(nelements, sizeof(void*)); }


/** \} \name Stack misc functions 
 \{ */

/** Test whether stack is empty.
   @result 1 if empty, 0 otherwise. */
static PHAST_INLINE
int stk_empty(Stack *s) { return lst_empty(s); }

/** Obtain size of stack.
   @result number of elements */
static PHAST_INLINE
int stk_size(Stack *s) { return lst_size(s); }


/** \} \name Stack push functions 
 \{ */

/** Push object on stack. */
static PHAST_INLINE
void stk_push(Stack *s, void *o) { lst_push(s, o); }

/** Push integer on stack. */
static PHAST_INLINE
void stk_push_int(Stack *s, int i) { lst_push_int(s, i); }

/** Push double on stack */
static PHAST_INLINE
void stk_push_dbl(Stack *s, double d) { lst_push_dbl(s, d); }

/** Push pointer on stack */
static PHAST_INLINE
void stk_push_ptr(Stack *s, void *ptr) { lst_push_ptr(s, ptr); }


/** \} \name Stack peek functions 
 \{ */

/** Peek at object at top of stack.
   @result address of object at top of stack, or NULL if stack is empty. 
   @note Object is not removed. 
   @see stk_pop */
static PHAST_INLINE
void* stk_peek(Stack* s) {
  if (stk_empty(s)) return NULL;
  return lst_arr_get(s, s->ridx-1);
}

/** Peek at integer at top of stack.
   @result integer from top of stack or 0 if stack is empty (use stk_empty).

   @warning return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when stack is empty. 
   @note Integer is not removed. 
   @see stk_empty */
static PHAST_INLINE
int stk_peek_int(Stack *s) {
  int *ptr = (int*)stk_peek(s);
  return (ptr == NULL ? 0 : *ptr);
}

/** Peek at double at top of stack.
   @result double from top of stack or 0 if stack is empty (use stk_empty).
   @warning return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when stack is empty. 
   @note Double is not removed
   @see stk_empty*/
static PHAST_INLINE
double stk_peek_dbl(Stack *s) {
  double *ptr = (double*)stk_peek(s);
  return (ptr == NULL ? 0 : *ptr);
}

/** Peek at pointer at top of stack.
   @result pointer from top of stack or NULL if stack is empty. 
   @note pointer is not removed */
static PHAST_INLINE
void* stk_peek_ptr(Stack *s) {
  void **ptr = (void**)stk_peek(s);
  return (ptr == NULL ? NULL : *ptr);
}



/** \} \name Stack pop functions 
 \{ */

/** Pop object from stack.
   @result address of object at top of stack, or NULL if stack is
   empty. 
   @warning Object is removed. 
   @note If object is an int, you will need to access it with
   *((int*)lst_get(...)); 
   if object is a Node* and you seek the
   attribute called "data", you will need to do
   ((Node*)lst_get(...))->data. */
static PHAST_INLINE
void* stk_pop(Stack* s) {
  if (stk_empty(s)) return NULL;
  return lst_arr_get(s, --s->ridx);
}

/** Pop integer from stack.
   @result integer from top of stack, or 0 if stack is empty (use stk_empty).
   @warning Integer is removed.
   @warning return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when stack is empty. 
   @see stk_empty */
static PHAST_INLINE
int stk_pop_int(Stack *s) {
  int *ptr = (int*)stk_pop(s);
  return (ptr == NULL ? 0 : *ptr);
}

/** Pop double from stack.
   @result double from top of stack or 0 if stack is empty (use stk_empty).
   @warning Double is removed
   @warning return value will be ambiguous when using numeric data
   containing zeroes.  Use stk_empty to tell when stack is empty. 
   @see stk_empty */
static PHAST_INLINE
double stk_pop_dbl(Stack *s) {
  double *ptr = (double*)stk_pop(s);
  return (ptr == NULL ? 0 : *ptr);
}

/** Pop pointer from stack.
   @result pointer from top of stack or NULL if stack is empty. 
   @warning Pointer is removed */
static PHAST_INLINE
void* stk_pop_ptr(Stack *s) {
  void **ptr = (void**)stk_pop(s);
  return (ptr == NULL ? NULL : *ptr);
}

/** \} \name Stack cleanup functions 
 \{ */

/** Clear contents of stack.
   @note Contents will be cleared but memory will remain allocated unlike stk_free
   @see stk_free */   
static PHAST_INLINE
void stk_clear(Stack *s) { lst_clear(s); }

/** Free memory associated with stack.
   @note Contents and Stack object itself will be freed unlike stk_clear. 
   @warning if the stack contains pointers, the objects pointed to are not freed
   @see stk_clear*/
static PHAST_INLINE
void stk_free(Stack *s) { lst_free(s); }
/** \} */

#endif

