/* Simple array-based stacks and supporting functions. 

   See documentation for lists.h.

   $Id: stacks.c,v 1.1.1.1 2004-06-03 22:43:11 acs Exp $
   Written by Adam Siepel, Spring 2001, Summer 2002
   Copyright 2001, 2002, Adam Siepel, University of California
*/

#include <stdlib.h>
#include "stacks.h"

int stk_empty(Stack *s) { return lst_empty(s); }

void* stk_pop(Stack* s) {
  if (stk_empty(s)) return NULL;
  return lst_arr_get(s, --s->ridx);
}

void* stk_peek(Stack* s) {
  if (stk_empty(s)) return NULL;
  return lst_arr_get(s, s->ridx-1);
}

int stk_pop_int(Stack *s) {
  int *ptr = (int*)stk_pop(s);
  return (ptr == NULL ? 0 : *ptr);
}

double stk_pop_dbl(Stack *s) {
  double *ptr = (double*)stk_pop(s);
  return (ptr == NULL ? 0 : *ptr);
}

void* stk_pop_ptr(Stack *s) {
  void **ptr = (void**)stk_pop(s);
  return (ptr == NULL ? NULL : *ptr);
}

int stk_peek_int(Stack *s) {
  int *ptr = (int*)stk_peek(s);
  return (ptr == NULL ? 0 : *ptr);
}

double stk_peek_dbl(Stack *s) {
  double *ptr = (double*)stk_peek(s);
  return (ptr == NULL ? 0 : *ptr);
}

void* stk_peek_ptr(Stack *s) {
  void **ptr = (void**)stk_peek(s);
  return (ptr == NULL ? NULL : *ptr);
}

void stk_push(Stack *s, void *o) { lst_push(s, o); }

void stk_push_int(Stack *s, int i) { lst_push_int(s, i); }

void stk_push_dbl(Stack *s, double d) { lst_push_dbl(s, d); }

void stk_push_ptr(Stack *s, void *ptr) { lst_push_ptr(s, ptr); }

int stk_size(Stack *s) { return lst_size(s); }

Stack* stk_new(int nelements, /* Starting number of elements */
	      int elementsz)	    /* Size of each element (bytes) */
{ return lst_new(nelements, elementsz); }

Stack* stk_new_int(int nelements) /* Starting number of elements */
{ return lst_new(nelements, sizeof(int)); }

Stack* stk_new_dbl(int nelements) /* Starting number of elements */
{ return lst_new(nelements, sizeof(double)); }

Stack* stk_new_ptr(int nelements) 
{ return lst_new(nelements, sizeof(void*)); }

void stk_clear(Stack *s) { lst_clear(s); }

void stk_free(Stack *s) { lst_free(s); }
