/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
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

#include <string.h>
#include <stdlib.h> 
#include <external_libs.h>

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
/** List type. */
typedef struct lst_struct List;

/** List sorting */
typedef enum {
  ASCENDING,			/**< sort in ascending order */
  DESCENDING			/**< sort in descending order*/
} order_t;


int lst_int_compare_asc(const void* ptr1, const void* ptr2);
int lst_int_compare_desc(const void* ptr1, const void* ptr2);
int lst_dbl_compare_asc(const void* ptr1, const void* ptr2);
int lst_dbl_compare_desc(const void* ptr1, const void* ptr2);

void *srealloc(void *ptr, size_t size);



static PHAST_INLINE
void lst_arr_set(List *l, int i, void *o) {
  if (l->elementsz <= sizeof(void*))
    l->array[i] = *((void**)o);	/* ?? */
  else
    memcpy(&l->array[i * l->step], o, l->elementsz);
}

static PHAST_INLINE
void* lst_arr_get(List *l, int i) {
  return (&l->array[i * l->step]);
}


/** Create new list.
  Returns newly allocated List object, with specified starting size. 

  @param nelements Starting number of elements.
  @param Size of each element (bytes).
  
  \sa lst_new_int, lst_new_dbl, lst_new_ptr, lst_push.
*/
List* lst_new(int nelements,	/* Starting number of elements */
	      int elementsz);	/* Size of each element (bytes) */

/** Create new list of integers. 
  Returns newly allocated List object, with starting size fixed at sizeof(int).

  \code 
  // same as
  lst_new(nelements, sizeof(int))
  \endcode
   
  @param nelements Starting number of elements.
  
  \sa lst_new, lst_new_dbl, lst_new_ptr.
*/
List* lst_new_int(int nelements); /* Starting number of elements */

/** Create new list of doubles.
  Returns newly allocated List object, with starting size fixed at sizeof(double).

  \code 
  // same as
  lst_new(nelements, sizeof(double))
  \endcode

  @param nelements Starting number of elements.
  \sa lst_new, lst_new_int, lst_new_ptr.  
*/
List* lst_new_dbl(int nelements); /* Starting number of elements */

/** Create new list of pointers.
  Returns newly allocated List object, with starting size fixed at sizeof(void*).

  \code 
  // same as
  lst_new(nelements, sizeof(void*))
  \endcode

  @param nelements Starting number of elements.
  \sa lst_new, lst_new_int, lst_new_dbl.
*/
List* lst_new_ptr(int nelements); /* Starting number of elements */

/** Free memory associated with list.
   List object itself will be freed also. 
   
   \note If the list contains pointers to other objects, as in the case of a list
   created with lst_new_ptr(), the objects targeted by the pointers will not be
    freed by this function.
   \endnote
   
   @param l List to be freed.
   
   \sa lst_push.
*/
void lst_free(List* l);		

/** Copy entire contents of list.
   elementsz bytes will be copied for each element. Destination list
   must already be initialized with proper elementsz. 

  @param dest Destination list.
  @param src Source list.
*/
void lst_cpy(List* dest, List* src);


/** Obtain size of list (number of elements).
   Returns number of elements 

  @param l Target list.
*/
static PHAST_INLINE int lst_size(List *l)
{ return l->ridx - l->lidx;}

/** Test whether list is empty.
   Returns 1 if empty, 0 otherwise 

  @param l Target list.

  \sa lst_clear.
*/
static PHAST_INLINE
int lst_empty(List *l) 
{ return (l->lidx >= l->ridx);}


/** Free the elements of a list of strings.

  \warning This will not free the list! Use with care!
  
  @param l Target list containing strings.
  \sa lst_free
*/
void lst_free_strings(List *l);

/** Reverse the order of a list.

  @param l Target list.  
*/
void lst_reverse(List *l);

/** \name List get & set operations. */

/** \{ */

/** Retrieve ith object in list. 
   Returns address of element in list, or NULL if index is out of
   bounds.  
   
   If object is an int, you will need to access it with *((int*)lst_get(...)); 
   \code
   // adding an int to a list of ints
   List * lst = lst_new_int(10);
   int a;
   lst_push_int(lst, 5);
   
   // retrieving the element
   a = *((int*) lst_get(lst, 0));

   // alternative way using helper functions
   a = lst_get_int(lst, 0);
   \endcode
   
      
  if object is a Node pointer and you seek the attribute called "data", you will need
    to do (*((Node**)lst_get(...)))->data;
  \code
  // adding a Node pointer to a list of pointers
  List * lst = lst_new_ptr(10);
  Node * nptr = ...; // some initialization
  lst_push_ptr(lst, nptr);

  // accessing the data element
  data = (*((Node**)lst_get(lst, 0)))->data;
  
  // alternative way using helper functions
  data = ((Node*)lst_get_ptr(lst, 0))->data;
  \endcode
  
  if object is a Node and you seek the attribute called "data", you will need
   to do ((Node*)lst_get(...))->data;
  \code
  // adding a copy of a Node to a list of Nodes
  List * lst = lst_new(10, sizeof(Node));
  Node n = ...; // some initialization
  lst_push(lst, (void*)(&n));

  data = ((Node*)lst_get(lst, 0))->data;
  // there is no alternative way to do this
  \endcode

  @param l List containing the desired object.
  @param i Index of the object to be retrieved.   
  
  \sa lst_new, lst_push, lst_get_int, lst_get_dbl, lst_get_ptr.
*/
static PHAST_INLINE
void* lst_get(List *l, int i) {
  if (i >= lst_size(l)) return NULL;
  return lst_arr_get(l, l->lidx + i);
}

/** Retrieve ith integer in list .
   Returns integer at ith position in list or 0 if i is out of bounds.

  \warning Return value will be ambiguous when using numeric data
   containing zeroes.  Make sure index is within bounds
   
  @param l List containing the desired object.
  @param i Index of the object to be retrieved.   

  \sa lst_get, lst_get_dbl, lst_get_ptr.
*/
static PHAST_INLINE
int lst_get_int(List* l, int i) {
  int *ptr = (int*)lst_get(l, i);
  return (ptr == NULL ? 0 : *ptr);
}

/** Retrieve ith double in list .
   Returns double at ith position in list or 0 if i is out of bounds.

  \warning Return value will be ambiguous when using numeric data
   containing zeroes.  Make sure index is within bounds 

  @param l List containing the desired object.
  @param i Index of the object to be retrieved.   

  \sa lst_get, lst_get_int, lst_get_ptr.
*/
static PHAST_INLINE
double lst_get_dbl(List* l, int i) {
  double *ptr =  (double*)lst_get(l, i); 
  return (ptr == NULL ? 0 : *ptr);
}

/** Retrieve ith pointer in list .
   Returns pointer at ith position in list or NULL if i is out of bounds. 

  @param l List containing the desired object.
  @param i Index of the object to be retrieved.   

  \sa lst_get, lst_get_int, lst_get_ptr.
*/
static PHAST_INLINE
void* lst_get_ptr(List* l, int i) {
  void **ptr =  (void**)lst_get(l, i); 
  return (ptr == NULL ? NULL : *ptr);
}

/** Set value of ith object in list. 

  Replace ith object (elementsz bytes will be copied).
  
  See lst_push() or lst_get() for details object argument.
  
  @param l List where the element is to be copied into.
  @param i Index of the object to be replaced.
  @param o Pointer to the object to be copied into the list.

  \sa lst_set_int, lst_set_dbl, lst_set_ptr, lst_get.
*/
static PHAST_INLINE
void lst_set(List *l, int i, 
	     void *o) {		/* Pointer to object to be copied into
				   list (see details under lst_push) */
  lst_arr_set(l, l->lidx + i, o);
}

/** Set value of ith integer in list. 

  @param l List where the element is to be copied into.
  @param idx Index of the object to be replaced.
  @param i Integer value.

  \sa lst_set, lst_set_dbl, lst_set_ptr.
*/
static PHAST_INLINE
void lst_set_int(List *l, int idx, int i)
{  lst_set(l, idx, &i); }


/** Set value of ith double in list. 

  @param l List where the element is to be copied into.
  @param idx Index of the object to be replaced.
  @param d Double value.

  \sa lst_set, lst_set_dbl, lst_set_ptr.
*/
static PHAST_INLINE
void lst_set_dbl(List *l, int idx, double d)
{  lst_set(l, idx, &d); }


/** Set value of ith pointer in list. 

  @param l List where the element is to be copied into.
  @param idx Index of the object to be replaced.
  @param ptr Pointer value.

  \sa lst_set, lst_set_int, lst_set_ptr.
*/
static PHAST_INLINE
void lst_set_ptr(List *l, int idx, void *ptr)
{  lst_set(l, idx, &ptr); }



/** \} */

/** \name List insert & delete operations. */

/** \{ */

/** Push object onto end of list.

Pointer to object to be copied into list (elementsz bytes).

e.g., if object is an int i, o must equal (void*)(&i); 
\code
// adding an int to a list of ints
List * lst = lst_new(10, sizeof(int));
int i = 2;
lst_push(lst, (void*)(&i));

// alternative way using helper functions
List * lst = lst_new_int(10);
lst_push_int(lst, 2);
\endcode

e.g., if object is a Node* nptr, o must equal (void*)(&nptr);
\code
// adding a Node pointer to a list of pointers
List * lst = lst_new(10, sizeof(void*));
Node * nptr = ...; // some initialization
lst_push(lst, (void*)(&nptr));

// alternative way using helper functions
List * lst = lst_new_ptr(10);
Node * nptr = ...; // some initialization
lst_push_ptr(lst, nptr);
\endcode

e.g., if object is a Node n, o must equal (void*)(&n).
\code
// adding a copy of a Node to a list of Nodes
List * lst = lst_new(10, sizeof(Node));
Node n = ...; // some initialization
lst_push(lst, (void*)(&n));

// there is no alternative way to do this
\endcode

\note The last example shows a list that contains the structures within the
list memory. There is no need to free the node elements when freeing the list.
Moreover, changing the contents of the variable #n will not affect the value
of the node within the list.
\endnote

  @param l List where object is to be appended.
  @param o Object pointer.
  
  \sa lst_push_int, lst_push_dbl, lst_push_ptr.
 */
static PHAST_INLINE 
void lst_push(List *l, 
              void* o) {	/* Pointer to object to be copied into
                               list (elementsz bytes) -- e.g., if
                               object is an int i, o must equal
                               (void*)(&i); if object is a Node* nptr,
                               o must equal (void*)(&nptr); if object
                               is a Node n, o must equal
                               (void*)(&n). */

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
    }
  }
  lst_arr_set(l, l->ridx++, o);
}
                               
/** Push integer onto end of list.

  @param l List where object is to be appended.
  @param i Integer value.
  
  \sa lst_push, lst_push_dbl, lst_push_ptr.
*/
static PHAST_INLINE
void lst_push_int(List *l, int i) 
{  lst_push(l, &i); }


/** Push double onto end of list .
  @param l List where object is to be appended.
  @param d Double value.
  
  \sa lst_push, lst_push_int, lst_push_ptr.
*/
static PHAST_INLINE
void lst_push_dbl(List *l, double d) 
{  lst_push(l, &d); }

/** Push pointer onto end of list.

  @param l List where object is to be appended.
  @param ptr Pointer value.
  
  \sa lst_push, lst_push_int, lst_push_dbl.
*/
static PHAST_INLINE
void lst_push_ptr(List *l, void *ptr)
{  lst_push(l, &ptr); }


/** Insert object after given index.

  Insert *following* specified index; use -1 to insert at front of list.
  
  @param l Target list.
  @param idx Index after which the object is to be inserted.
  @param o Pointer to object (see lst_push()).
  
  \sa lst_push, lst_insert_idx_int, lst_insert_idx_dbl, lst_insert_idx_ptr.
*/
int lst_insert_idx(List *l, int idx, void *o);

/** Insert integer after given index.

  Insert *following* specified index; use -1 to insert at front of list.
  
  @param l Target list.
  @param idx Index after which the object is to be inserted.
  @param i Integer value.
  
  \sa lst_push, lst_insert_idx, lst_insert_idx_dbl, lst_insert_idx_ptr.
*/
void lst_insert_idx_int(List *l, int idx, int i);
/** Insert double after given index.

  Insert *following* specified index; use -1 to insert at front of list.
  
  @param l Target list.
  @param idx Index after which the object is to be inserted.
  @param d Double value.
  
  \sa lst_push, lst_insert_idx, lst_insert_idx_int, lst_insert_idx_ptr.
*/
void lst_insert_idx_dbl(List *l, int idx, double d);

/** Insert integer after given index.

  Insert *following* specified index; use -1 to insert at front of list.
  
  @param l Target list.
  @param idx Index after which the object is to be inserted.
  @param ptr Pointer value.
  
  \sa lst_push, lst_insert_idx, lst_insert_idx_int, lst_insert_idx_dbl.
*/
void lst_insert_idx_ptr(List *l, int idx, void *ptr);

/** Delete ith object in list. 
   Returns 0 on successful deletion, 1 otherwise. 

   \warning Deletion is highly inefficient. 

  @param l Target list.
  @param i Index of object to be deleted.
  
  \sa lst_delete_obj, lst_delete_obj_compare.
*/
int lst_delete_idx(List *l, int i);

/** Delete first matching object in list.
   Performs direct comparison of list elements (perhaps not desirable if
   pointers; see below).
   Returns 0 on successful deletion, 1 otherwise. 

   \warning Deletion is highly inefficient. 

  @param l Target list.
  @param o Pointer to object used for comparison (see lst_push() for details).
  
  \sa lst_delete_idx, lst_delete_obj_compare.
*/
int lst_delete_obj(List *l, 
		   void *o);	/* Pointer to object to be compared
				   (see details under lst_push) */

/** Delete first matching object in list.
   Compares objects using specified comparison function.
   Returns 0 on successful deletion, 1 otherwise. 

   \warning Deletion is highly inefficient. 

  @param l Target list.
  @param o Pointer to object used for comparison (see lst_push() for details).
  @param compare Comparision function: arguments are pointers to objects in the list.
                                   For example, if the list contains
                                   ints, then the function should
                                   compare "*((int*)arg1)" and
                                   "*((int*)arg2)". If it contains
                                   pointers to some object Node, then
                                   the function might want to compare
                                   "**((Node**)arg1)" and
                                   "**((Node**)arg2)".  
  
  \sa lst_delete_idx, lst_delete_obj.
*/
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


/** Clear contents of list.
   Contents will be cleared but memory will remain allocated.

  @param l Target list.
*/
static PHAST_INLINE
void lst_clear(List* l) 
{  l->ridx = l->lidx = 0; }


/** \} */

/** \name List sort & search operations. */

/** \{ */

/** Sort list using qsort.
   Compares objects using specified comparison function. 

  @param l Target list.
  @param compare Comparison function: arguments are pointers to objects in the list
				   (see details under lst_delete_obj_compare()). 

  \sa lst_qsort_int, lst_qsort_dbl.
*/
void lst_qsort(List *l, 
	       int (*compare)(const void *, const void *)
				/* Comparison function: arguments are
				   pointers to objects in the list
				   (see details under
				   lst_delete_obj_compare). */
	       );

/** Sort list of integers using qsort. 

  @param l Target list.
  @param ord Sorting order.

  \sa lst_qsort_int, lst_qsort_dbl.
*/
void lst_qsort_int(List *l, 
		   order_t ord); /* Sorting order */

/** Sort list of doubles using qsort. 
  @param l Target list.
  @param ord Sorting order.

  \sa lst_qsort, lst_qsort_int.
*/
void lst_qsort_dbl(List *l, 
		   order_t ord); /* Sorting order */

/** Search for object in list (linear).
   Performs direct comparison of list elements (perhaps not desirable if
   pointers; see below).
   Returns index of match or -1 if no match 

  @param l Target list.
  @param ptr Pointer to object to be compared (see details under lst_push()).

  \sa lst_find_compare, lst_bsearch_int.
*/
int lst_find(List *l, 
	     void* ptr);	/* Pointer to object to be compared
				   (see details under lst_push) */

/** Search for object in list (linear).
   Compares objects using specified comparison function.
   Returns index of match or -1 if no match 

  @param l Target list.
  @param ptr Pointer to object to be compared (see details under lst_push()).
  @param compare Comparison function: arguments are pointers to objects in the list
				   (see details under lst_delete_obj_compare()).

  \sa lst_find, lst_bsearch_int.
*/
int lst_find_compare(List *l, 
		     void* ptr, /* Pointer to object to be compared
				   (see details under lst_push) */
		     int (*compare)(void*, void*)
				/* Comparison function: arguments are
				   pointers to objects in the list
				   (see details under
				   lst_delete_obj_compare). */
		     );

/** Binary search of a list of integers.

   Expects list to be in ascending order.  When query falls between
   two values, the index of the smaller one is returned.  Returns -1
   if query is smaller than all values in the list. 
   
  @param lst Target list.
  @param val Integer value to be searched.   
*/
int lst_bsearch_int(List *lst, int val);


/** \} */

/** \name List math operations. */

/** \{ */

/** Compute the mean of a list of doubles.

  @param l List of doubles.
*/
double lst_dbl_mean(List *l);
/** Compute the standard deviation of a list of doubles.

  Computes the standard deviation of a sample \f$ \frac{\sum_{i=1}^N x_i}{N - 1} \f$.
  
  @param l List of doubles.
*/
double lst_dbl_stdev(List *l);

/** Compute the quantiles of a list of doubles.

  \warning Assumes list is sorted in ascending order.
  
  @param l List of doubles.
  @param quantiles Quantiles to be computed; values must be between 0 and 1.
  @param nquantiles Number of quantiles to be computed.
  @param quantile_vals Pre-allocated array where quantile values are to be stored.
*/
void lst_dbl_quantiles(List *l, double *quantiles, int nquantiles, 
                       double *quantile_vals);

/** \} */

static PHAST_INLINE
int lst_find_int(List *l, int i) 
{ return lst_find(l, &i); }

static PHAST_INLINE
int lst_find_dbl(List *l, double d) 
{ return lst_find(l, &d); }

static PHAST_INLINE
int lst_find_ptr(List *l, void *ptr) 
{ return lst_find(l, &ptr); }

#endif
