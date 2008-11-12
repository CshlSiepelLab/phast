/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: lists.h,v 1.6 2008-11-12 02:07:59 acs Exp $ */

/** \file lists.h
   Simple array-based lists and supporting functions. 

   Supports storage of objects of arbitrary size.  Convenience
   functions are available for common data types such as ints,
   doubles, and pointers.  The "stacks" and "queues" libraries are
   layered on top of this one.  Use these lists when memory locality
   is important or when you need to access elements by index; use
   linked-lists instead when (mid-list) insertions and deletions are
   important. 
   \ingroup base
*/

/*
   TODO
   * generalize lst_bsearch to work with any kind of object and with
   specified comparison function
   * add insertion functions to complement deletion functions
   * create a linked-list alternative, either behind the same
   interface (implementation strategy could be selected as parameter
   to lst_new()) or as a separate library
*/

#ifndef LISTS_H
#define LISTS_H

#include <stdlib.h> 
#include <string.h>

/** Basic List object */
struct lst_struct {
  void** array;                 /**< storage array for list elements */
  int lidx;                     /**< leftmost index of active list
                                   (inclusive) */
  int ridx;                     /**< righmost index of active list
                                   (exclusive) */
  int CAPACITY;                 /**< current storage capacity (number
                                   of elements) */
  int elementsz;                /**< number of bytes for each list element */
  int step;                     /**< number of array elements occupied
                                   by each list element (multiple are
                                   possible)  */ 
};
typedef struct lst_struct List;

/** List sorting */
typedef enum {
  ASCENDING,			/**< sort in ascending order */
  DESCENDING			/**< sort in descending order*/
} order_t;


void lst_arr_set(List *l, int i, void *o);
void* lst_arr_get(List *l, int i);
int lst_int_compare_asc(const void* ptr1, const void* ptr2);
int lst_int_compare_desc(const void* ptr1, const void* ptr2);
int lst_dbl_compare_asc(const void* ptr1, const void* ptr2);
int lst_dbl_compare_desc(const void* ptr1, const void* ptr2);


/* Create new list.
   Returns newly allocated List object, with specified starting size. */
List* lst_new(int nelements,	/* Starting number of elements */
	      int elementsz);	/* Size of each element (bytes) */

/* Free memory associated with list.
   List object itself will be freed also. */
void lst_free(List* l);		

/* Push object onto end of list */
void lst_push(List *l, 
              void* o);		/* Pointer to object to be copied into
                               list (elementsz bytes) -- e.g., if
                               object is an int i, o must equal
                               (void*)(&i); if object is a Node* nptr,
                               o must equal (void*)(&nptr); if object
                               is a Node n, o must equal
                               (void*)(&n). */

/* Retrieve ith object in list. 
   Returns address of element in list, or NULL if index is out of
   bounds.  If object is an int, you will need to access it with
   *((int*)lst_get(...)); if object is a Node* and you seek the
   attribute called "data", you will need to do
   (*((Node**)lst_get(...)))->data. */
void* lst_get(List *l, int i);

/* Set value of ith object in list. 
   elementsz bytes will be copied. */
void lst_set(List *l, int i, 
	     void *o);		/* Pointer to object to be copied into
				   list (see details under lst_push) */

/* Copy entire contents of list.
   elementsz bytes will be copied for each element. Destination list
   must already be initialized with proper elementsz. */
void lst_cpy(List* dest, List* src);

/* Delete ith object in list. 
   Returns 0 on successful deletion, 1 otherwise. 

   Warning: Deletion is highly inefficient. */
int lst_delete_idx(List *l, int i);

/* Delete first matching object in list.
   Performs direct comparison of list elements (perhaps not desirable if
   pointers; see below).
   Returns 0 on successful deletion, 1 otherwise. 

   Warning: Deletion is highly inefficient. */
int lst_delete_obj(List *l, 
		   void *o);	/* Pointer to object to be compared
				   (see details under lst_push) */

/* Delete first matching object in list.
   Compares objects using specified comparison function.
   Returns 0 on successful deletion, 1 otherwise. 

   Warning: Deletion is highly inefficient. */
int lst_delete_obj_compare(List *l, 
			   void *o, /* Pointer to object to be
				       compared (see details under
				       lst_push) */
			   int (*compare)(void*, void*)
                                /* Comparison function: arguments are
                                   pointers to objects in the list.
                                   For example, if the list contains
                                   ints, then the function should
                                   compare "*((int*)arg1)" and
                                   "*((int*)arg2)". If it contains
                                   pointers to some object Node, then
                                   the function might want to compare
                                   "**((Node**)arg1)" and
                                   "**((Node**)arg2)".  */
			   );

int lst_insert_idx(List *l, int idx, void *o);
void lst_insert_idx_int(List *l, int idx, int i);
void lst_insert_idx_dbl(List *l, int idx, double d);
void lst_insert_idx_ptr(List *l, int idx, void *ptr);


/* Search for object in list (linear).
   Performs direct comparison of list elements (perhaps not desirable if
   pointers; see below).
   Returns index of match or -1 if no match */
int lst_find(List *l, 
	     void* ptr);	/* Pointer to object to be compared
				   (see details under lst_push) */

/* Search for object in list (linear).
   Compares objects using specified comparison function.
   Returns index of match or -1 if no match */
int lst_find_compare(List *l, 
		     void* ptr, /* Pointer to object to be compared
				   (see details under lst_push) */
		     int (*compare)(void*, void*)
				/* Comparison function: arguments are
				   pointers to objects in the list
				   (see details under
				   lst_delete_obj_compare). */
		     );

/* Sort list using qsort.
   Compares objects using specified comparison function. */
void lst_qsort(List *l, 
	       int (*compare)(const void *, const void *)
				/* Comparison function: arguments are
				   pointers to objects in the list
				   (see details under
				   lst_delete_obj_compare). */
	       );


/* Clear contents of list.
   Contents will be cleared but memory will remain allocated. */   
void lst_clear(List* l);

/* Obtain size of list (number of elements).
   Returns number of elements */
int lst_size(List *l);

/* Test whether list is empty.
   Returns 1 if empty, 0 otherwise */
int lst_empty(List *l); 

/* Push integer onto end of list */
void lst_push_int(List *l, int i);

/* Push double onto end of list */
void lst_push_dbl(List *l, double d);

/* Push pointer onto end of list */
void lst_push_ptr(List *l, void *ptr);

/* Retrieve ith integer in list .
   Returns integer at ith position in list or 0 if i is out of bounds.

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Make sure index is within bounds */
int lst_get_int(List* l, int i);

/* Retrieve ith double in list .
   Returns double at ith position in list or 0 if i is out of bounds.

   Warning: return value will be ambiguous when using numeric data
   containing zeroes.  Make sure index is within bounds */
double lst_get_dbl(List* l, int i);

/* Retrieve ith pointer in list .
   Returns pointer at ith position in list or NULL if i is out of bounds. */
void* lst_get_ptr(List* l, int i);

/* Set value of ith integer in list. */
void lst_set_int(List *l, int idx, int i);

/* Set value of ith double in list. */
void lst_set_dbl(List *l, int idx, double d);

/* Set value of ith pointer in list. */
void lst_set_ptr(List *l, int idx, void *ptr);

int lst_bsearch_int(List *lst, int val);

/* Sort list of integers using qsort. */
void lst_qsort_int(List *l, 
		   order_t ord); /* Sorting order */

/* Sort list of doubles using qsort. */
void lst_qsort_dbl(List *l, 
		   order_t ord); /* Sorting order */



/* Create new list of integers. 
   Like lst_new, but with element size fixed at sizeof(int).
   Returns newly allocated List object, with specified starting size. */
List* lst_new_int(int nelements); /* Starting number of elements */

/* Create new list of doubles.
   Like lst_new, but with element size fixed at sizeof(double).
   Returns newly allocated List object, with specified starting size. */
List* lst_new_dbl(int nelements); /* Starting number of elements */

/* Create new list of pointers.
   Like lst_new, but with element size fixed at sizeof(void*).
   Returns newly allocated List object, with specified starting size. */
List* lst_new_ptr(int nelements); /* Starting number of elements */

void lst_free_strings(List *l);

void lst_reverse(List *l);

double lst_dbl_mean(List *l);
double lst_dbl_stdev(List *l);
void lst_dbl_quantiles(List *l, double *quantiles, int nquantiles, 
                       double *quantile_vals);


/***************************************************************************/
/* Inline functions; note: these are redefined in lists.c for              */
/* situations in which inlining is not available                           */
/***************************************************************************/

extern inline
void lst_arr_set(List *l, int i, void *o) {
  if (l->elementsz <= sizeof(void*))
    l->array[i] = *((void**)o);	/* ?? */
  else
    memcpy(&l->array[i * l->step], o, l->elementsz);
}

extern inline
void* lst_arr_get(List *l, int i) {
  return (&l->array[i * l->step]);
}

extern inline
int lst_size(List *l) 
{  return (l->ridx - l->lidx); }

extern inline
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
      l->array = (void**)realloc(l->array, l->CAPACITY * l->elementsz);
    }
  }
  lst_arr_set(l, l->ridx++, o);
}

extern inline
void* lst_get(List *l, int i) {
  if (i >= lst_size(l)) return NULL;
  return lst_arr_get(l, l->lidx + i);
}

extern inline
void lst_set(List *l, int i, void *o) {
  lst_arr_set(l, l->lidx + i, o);
}

extern inline
void lst_clear(List* l) 
{  l->ridx = l->lidx = 0; }

extern inline
int lst_empty(List *l) 
{  return (l->lidx >= l->ridx); }

extern inline
void lst_push_int(List *l, int i) 
{  lst_push(l, &i); }

extern inline
void lst_push_dbl(List *l, double d) 
{  lst_push(l, &d); }

extern inline
void lst_push_ptr(List *l, void *ptr) 
{  lst_push(l, &ptr); }

extern inline
int lst_get_int(List* l, int i) {  
  int *ptr = (int*)lst_get(l, i);
  return (ptr == NULL ? 0 : *ptr);
}

extern inline
double lst_get_dbl(List* l, int i) { 
  double *ptr =  (double*)lst_get(l, i); 
  return (ptr == NULL ? 0 : *ptr);
}

extern inline
void* lst_get_ptr(List *l, int i) {
  void **ptr =  (void**)lst_get(l, i); 
  return (ptr == NULL ? NULL : *ptr);
}

extern inline
void lst_set_int(List *l, int idx, int i) 
{  lst_set(l, idx, &i); }

extern inline
void lst_set_dbl(List *l, int idx, double d) 
{  lst_set(l, idx, &d); }

extern inline
void lst_set_ptr(List *l, int idx, void *ptr) 
{  lst_set(l, idx, &ptr); }

extern inline
int lst_find_int(List *l, int i) 
{ return lst_find(l, &i); }

extern inline
int lst_find_dbl(List *l, double d) 
{ return lst_find(l, &d); }

extern inline
int lst_find_ptr(List *l, void *ptr) 
{ return lst_find(l, &ptr); }

#endif
