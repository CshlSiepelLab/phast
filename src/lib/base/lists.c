/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: lists.c,v 1.7 2008-11-12 02:07:59 acs Exp $ */

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
  if (l->array == NULL)
    die("ERROR lst_new l->array is NULL\n");
  l->step = ceil(l->elementsz * 1.0/sizeof(void*));
  return l;
}

void lst_free(List* l) {
  free(l->array);
  free(l);
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
  if ((idx = lst_find_compare(l, o, compare)) == -1)
    return 0;
  lst_delete_idx(l, idx);
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
    void *lastelement;

  /* push copy of last element onto end of list.  We want to use
     lst_push for automatic realloc (if necessary) but have to be
     careful to pass a *copy* of the element, in case there is a
     realloc (pointer will be stale) */
    if (l->elementsz <= sizeof(void*)) {
      lastelement = *((void**)lst_get(l, lst_size(l)-1));
      lst_push(l, &lastelement); 
    }
    else {
      lastelement = smalloc(l->elementsz);
      memcpy(lastelement, lst_get(l, lst_size(l)-1), l->elementsz);
      lst_push(l, lastelement);
      free(lastelement);
    }

    /* now shift each downstream element right by one */
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
  lst_clear(l);
}

/* reverse the order of a list */
void lst_reverse(List *l) {
  void *tmp = smalloc(l->elementsz);
  int i, j;
  for (i = 0, j = lst_size(l) - 1; i < j; i++, j--) {
    memcpy(tmp, lst_get(l, j), l->elementsz);
    lst_set(l, j, lst_get(l, i));
    lst_set(l, i, tmp);
  }
  free(tmp);
}

/* some simple statistical functions for lists of doubles */
double lst_dbl_mean(List *l) {
  int i;
  double sum = 0;
  for (i = 0; i < lst_size(l); i++)
    sum += lst_get_dbl(l, i);
  return sum / lst_size(l);
}

double lst_dbl_stdev(List *l) {
  int i;
  double sum = 0, mean = lst_dbl_mean(l);  
  for (i = 0; i < lst_size(l); i++)
    sum += pow(lst_get_dbl(l, i) - mean, 2);
  return sqrt(sum / (lst_size(l) - 1));
}

/* assumes list is sorted in ascending order */
void lst_dbl_quantiles(List *l, double *quantiles, int nquantiles, 
                       double *quantile_vals) {
  int i, max_idx = lst_size(l) - 1;
  for (i = 0; i < nquantiles; i++) {
    double lower_idx, upper_idx, lower_val, upper_val;
    if (!(quantiles[i] >= 0 && quantiles[i] <= 1))
      die("ERROR lst_dbl_quantiles: quantiles[%i]=%f, should be in [0,1]\n",
	  i, quantiles[i]);
    lower_idx = floor(quantiles[i] * max_idx);
    upper_idx = ceil(quantiles[i] * max_idx);
    lower_val = lst_get_dbl(l, lower_idx);
    upper_val = lst_get_dbl(l, upper_idx);

    if (lower_idx == upper_idx) 
      quantile_vals[i] = lower_val;
    else {                      /* linearly interpolate */
      double lower_q = (double)lower_idx / max_idx;
      double upper_q = (double)upper_idx / max_idx;
      double lambda = (quantiles[i] - lower_q) / (upper_q - lower_q);
      quantile_vals[i] = lower_val + lambda * (upper_val - lower_val);
    }
  }
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

/** Sort list of integers using qsort */
void lst_qsort_int(List *l, order_t ord) 
{ lst_qsort(l, ord == ASCENDING ? lst_int_compare_asc : 
            lst_int_compare_desc); }

void lst_qsort_dbl(List *l, order_t ord) 
{ lst_qsort(l, ord == ASCENDING ? lst_dbl_compare_asc : 
            lst_dbl_compare_desc); }

List* lst_new_int(int nelements) 
{ return lst_new(nelements, max(sizeof(int), sizeof(void*))); }
/* have to be sure size isn't smaller than size of void*; can be an
   issue with 64-bit architectures */

List* lst_new_dbl(int nelements) 
{ return lst_new(nelements, max(sizeof(double), sizeof(void*))); }

List* lst_new_ptr(int nelements) 
{ return lst_new(nelements, sizeof(void*)); }

