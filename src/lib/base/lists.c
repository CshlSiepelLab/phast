/* $Id: lists.c,v 1.2 2004-06-09 17:10:29 acs Exp $
   Written by Adam Siepel, Spring 2001 and Summer 2002
   Copyright 2001, 2002, Adam Siepel, University of California */

/** \file lists.c 
   Simple array-based lists and supporting functions.  
   Elements stored in contiguous memory locations (for efficient
   caching performance).  Supports storage of objects of arbitrary
   size.  Convenience functions are available for common data types
   such as ints, doubles, and pointers.  The "stacks" and "queues"
   libraries are layered on top of this one.  Use these lists when
   memory locality is important or when you need to access elements by
   index; use linked-lists instead when (mid-list) insertions and
   deletions are important. 
   \ingroup base
*/

/*
   TODO
   - generalize lst_bsearch_int to work with any kind of object and with
   specified comparison function
   - add insertion functions to complement deletion functions
   - create a linked-list alternative, either behind the same
   interface (implementation strategy could be selected as parameter
   to lst_new()) or as a separate library
*/


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "lists.h"
#include "stringsplus.h"
#include <misc.h>

List* lst_new(int nelements, int elementsz) {
  List *l = (List*)smalloc(sizeof(List));
  l->ridx = l->lidx = 0;
  l->CAPACITY = nelements;
  l->elementsz = elementsz;
  l->array = (void**)smalloc(nelements * elementsz);
  assert(l->array != NULL);
  l->step = ceil(l->elementsz * 1.0/sizeof(void*));
  return l;
}

void lst_free(List* l) {
  free(l->array);
  free(l);
}

void lst_push(List *l, void* o) {
  int i;
  if (l->ridx >= l->CAPACITY) {
    if (l->lidx > 0) {
      for (i = l->lidx; i < l->ridx; i++) 
        lst_arr_set(l, i - l->lidx, lst_arr_get(l, i));
      l->ridx -= l->lidx;
      l->lidx = 0;
    }

    else {
      l->CAPACITY *= 2;
      l->array = (void**)srealloc(l->array, l->CAPACITY * l->elementsz);
      assert(l->array != NULL);
    }
  }
  lst_arr_set(l, l->ridx++, o);
}

void* lst_get(List *l, int i) {
  if (i >= lst_size(l)) return NULL;
  return lst_arr_get(l, l->lidx + i);
}

void lst_set(List *l, int i, void *o) {
  assert(i >= 0 && i < lst_size(l));
  lst_arr_set(l, l->lidx + i, o);
}


/* low-level array access routine, encapsulates memcpy, helps simplify
   indexing.  For use internally only. */ 
void lst_arr_set(List *l, int i, void *o) {
  if (l->elementsz <= sizeof(void*))
    l->array[i] = *((void**)o);	
  else
    memcpy(&l->array[i * l->step], o, l->elementsz);
}

/* low-level array access routine, helps simplify indexing.  For use
 * internally only */
void* lst_arr_get(List *l, int i) {
  return (&l->array[i * l->step]);
}

void lst_cpy(List* dest, List* src) {
  int i;
  lst_clear(dest);
  for (i = 0; i < lst_size(src); i++)
    lst_push(dest, lst_get(src, i));  
}

int lst_delete_idx(List *l, int idx) {
  int i;
  if (idx >= lst_size(l))
    return 1;
  for (i = idx + 1; i < lst_size(l); i++)
    lst_arr_set(l, l->lidx + i - 1, lst_arr_get(l, l->lidx + i));
  l->ridx--;
  return 0;
}

int lst_delete_obj(List *l, void *o) {
  int idx;
  if ((idx = lst_find(l, o)) != -1)
    return lst_delete_idx(l, idx);
  return 1;
}

int lst_delete_obj_compare(List *l, void *o, int (*compare)(void*, void*)) {
  int idx;
  if ((idx = lst_find_compare(l, o, compare)) != -1)
    return lst_delete_idx(l, idx);
  return 1;
}

/* insert *following* specified index; use -1 to ins at front of list */
int lst_insert_idx(List *l, int idx, void *o) {
  int i;
  if (idx < -1 || idx >= lst_size(l))
    return 1;
  if (idx == lst_size(l) - 1) 
    lst_push(l, o);
  else {
    lst_push(l, lst_get(l, lst_size(l)-1));
    for (i = lst_size(l) - 2; i >= idx+2; i--)
      lst_arr_set(l, l->lidx + i, lst_arr_get(l, l->lidx + i - 1));
    lst_arr_set(l, l->lidx + idx + 1, o);
  }
  return 0;
}

void lst_insert_idx_int(List *l, int idx, int i) 
{  lst_insert_idx(l, idx, &i); }

void lst_insert_idx_dbl(List *l, int idx, double d) 
{  lst_insert_idx(l, idx, &d); }

void lst_insert_idx_ptr(List *l, int idx, void *ptr) 
{  lst_insert_idx(l, idx, &ptr); }

int lst_find(List *l, void* ptr) {
  int i;
  if (l->elementsz <= sizeof(void*)) {
    for (i = 0; i < lst_size(l); i++)
      if (*((void**)lst_get(l, i)) == *((void**)ptr))
	return i;
  }
  else {
    for (i = 0; i < lst_size(l); i++)
      if (memcmp(lst_get(l, i), ptr, l->elementsz) == 0)
	return i;
  }
  return -1;
}

int lst_find_compare(List *l, void* ptr, 
                     int (*compare)(void*, void*)) {
  int i;
  for (i = 0; i < lst_size(l); i++)
    if (compare(lst_get(l, i), ptr))
      return i;
  return -1;  
}


/* Expects list to be in ascending order.  When query falls between
   two values, the index of the smaller one is returned.  Returns -1
   if query is smaller than all values in the list. */

/* todo: support comparison functino, ASCENDING or DESCENDING order,
   option to report not found or report next lowest index (can always
   obtain next highest from this */
int lst_bsearch_int(List *lst, int val) {
  int l = 0;
  int r = lst_size(lst) - 1;
  if (r < 0 || val < lst_get_int(lst, l)) return -1;
  if (val > lst_get_int(lst, r)) return r;      
  while (l <= r) {
    int m = (l + r)/2;
    int candidate = lst_get_int(lst, m);
    if (val == candidate) 
      return m;
    else if (val < candidate) 
      r = m - 1;
    else                        /* val > candidate */
      l = m + 1;
  }

 /* the value's proper place must either be between l-1 and l or
  * between r and r+1. */
  if (val > lst_get_int(lst, r))
    return r;
  return(l - 1);
}

void lst_qsort(List *l, int (*compare)(const void *, const void *)) {
  /* will we need to use max(elementsz, sizeof(void*))? */
  qsort(&l->array[l->lidx], lst_size(l), l->elementsz, compare);
}

/* free the elements of a list of strings */
void lst_free_strings(List *l) {
  int i;
  for (i = 0; i < lst_size(l); i++) 
    if (lst_get_ptr(l, i) != NULL) str_free(lst_get_ptr(l, i));
}

/****************************************************************************/
/*          Comparison functions for use with numeric lists                 */
/****************************************************************************/

int lst_int_compare_asc(const void* ptr1, const void* ptr2) {
  int val1 = *((int*)ptr1);
  int val2 = *((int*)ptr2);
  return (val1 - val2);
}

int lst_int_compare_desc(const void* ptr1, const void* ptr2) {
  int val1 = *((int*)ptr1);
  int val2 = *((int*)ptr2);
  return (val2 - val1);
}

int lst_dbl_compare_asc(const void* ptr1, const void* ptr2) {
  double val1 = *((double*)ptr1);
  double val2 = *((double*)ptr2);
  if (val1 == val2) return 0;
  else if (val1 < val2) return -1;
  return 1;                    
}

int lst_dbl_compare_desc(const void* ptr1, const void* ptr2) {
  double val1 = *((double*)ptr1);
  double val2 = *((double*)ptr2);
  if (val1 == val2) return 0;
  else if (val1 < val2) return 1;
  return -1;                    
}


void lst_clear(List* l) 
{  l->ridx = l->lidx = 0; }

int lst_size(List *l) 
{  return (l->ridx - l->lidx); }

int lst_empty(List *l) 
{  return (l->lidx >= l->ridx); }

void lst_push_int(List *l, int i) 
{  lst_push(l, &i); }

void lst_push_dbl(List *l, double d) 
{  lst_push(l, &d); }

void lst_push_ptr(List *l, void *ptr) 
{  lst_push(l, &ptr); }

int lst_get_int(List* l, int i) {  
  int *ptr = (int*)lst_get(l, i);
  return (ptr == NULL ? 0 : *ptr);
}

double lst_get_dbl(List* l, int i) { 
  double *ptr =  (double*)lst_get(l, i); 
  return (ptr == NULL ? 0 : *ptr);
}

void* lst_get_ptr(List *l, int i) {
  void **ptr =  (void**)lst_get(l, i); 
  return (ptr == NULL ? NULL : *ptr);
}

void lst_set_int(List *l, int idx, int i) 
{  lst_set(l, idx, &i); }

void lst_set_dbl(List *l, int idx, double d) 
{  lst_set(l, idx, &d); }

void lst_set_ptr(List *l, int idx, void *ptr) 
{  lst_set(l, idx, &ptr); }

int lst_find_int(List *l, int val) 
{  return lst_find(l, (void*)(&val)); }

int lst_find_dbl(List *l, double val) 
{  return lst_find(l, (void*)(&val)); } 

int lst_find_ptr(List *l, void *ptr) 
{ return lst_find(l, &ptr); }

/** Sort list of integers using qsort */
void lst_qsort_int(List *l, order_t ord) 
{ lst_qsort(l, ord == ASCENDING ? lst_int_compare_asc : 
            lst_int_compare_desc); }

void lst_qsort_dbl(List *l, order_t ord) 
{ lst_qsort(l, ord == ASCENDING ? lst_dbl_compare_asc : 
            lst_dbl_compare_desc); }

List* lst_new_int(int nelements) 
{ return lst_new(nelements, sizeof(int)); }

List* lst_new_dbl(int nelements) 
{ return lst_new(nelements, sizeof(double)); }

List* lst_new_ptr(int nelements) 
{ return lst_new(nelements, sizeof(void*)); }

