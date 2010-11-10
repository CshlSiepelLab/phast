/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: hashtable.c,v 1.8 2009/02/19 19:44:15 agd27 Exp $*/

/* hashtable - Fast, simple array-based hash table, optimized for
   'put' and 'get' ('delete' somewhat inefficient).  Stores copies of
   keys but not of data objects, which are managed as void*s (memory
   management expected to be done externally) */

#include <stdlib.h>
#include <lists.h>
#include <hashtable.h>
#include <math.h>
#include <misc.h>

/* Create new hashtable with initial capacity as specified (in number
   of items).  
   Returns new hashtable with initial capacity as specified. */
Hashtable* hsh_new(int est_capacity) {
  Hashtable* ht;
  int i;
  ht = (Hashtable*)smalloc(sizeof(Hashtable));
  ht->nbuckets = ceil(est_capacity*1.0/LOADING_FACTOR);
  if (ht->nbuckets < 10) ht->nbuckets = 10;
  ht->keys = (List**)smalloc(ht->nbuckets*sizeof(List*));
  ht->vals = (List**)smalloc(ht->nbuckets*sizeof(List*));
  for (i = 0; i < ht->nbuckets; i++) 
    ht->keys[i] = ht->vals[i] = NULL;    
  return ht;
}

/* makes copy of hashtable.  Warning: if vals are pointers,
   only copies pointers.  Does copy keys. */
Hashtable *hsh_copy(Hashtable *src) {
  Hashtable *ht;
  int i, j, len;
  char *key, *keycpy;
  ht = (Hashtable*)smalloc(sizeof(Hashtable));
  ht->nbuckets = src->nbuckets;
  ht->keys = (List**)smalloc(ht->nbuckets*sizeof(List*));
  ht->vals = (List**)smalloc(ht->nbuckets*sizeof(List*));
  for (i=0; i<ht->nbuckets; i++) {
    if (src->keys[i] == NULL) {
      ht->keys[i] = NULL;
      ht->vals[i] = NULL;
    }
    else {
      len = lst_size(src->keys[i]);
      ht->keys[i] = lst_new_ptr(len);
      ht->vals[i] = lst_new_ptr(len);
      for (j=0; j<len; j++) {
	key = (char*)lst_get_ptr(src->keys[i], j);
	keycpy = smalloc(sizeof(char)*(strlen(key)+1));
	strcpy(keycpy, key);
	lst_push_ptr(ht->keys[i], keycpy);
	lst_push_ptr(ht->vals[i], lst_get_ptr(src->vals[i], j));
      }
    }
  }
  return ht;
}

void hsh_put_int(Hashtable *ht, const char *key, int val) {
  hsh_put(ht, key, int_to_ptr(val));
}

/* Equality function to pass to list_find_equivalent; simply tests for
   exact string match */
int equal(void *key1ptr, void *key2) {
  return !strcmp(*((char**)key1ptr), (char*)key2);
}

/* Retrieve object associated with specified key.
   Returns pointer to object or -1 if not found. 

   Warning: Convention of returning -1 when object is not found is
   inappropriate when objects are integers (needs to be fixed).*/
void* hsh_get(Hashtable* ht, const char *key) {
  unsigned int bucket;
  int idx;
  bucket = hsh_hash_func(ht, key);
  if (ht->keys[bucket] == NULL || 
      (idx = lst_find_compare(ht->keys[bucket], (void*)key, equal)) == -1) 
    return (void*)-1;
  return lst_get_ptr(ht->vals[bucket], idx);
}

int hsh_get_int(Hashtable *ht, const char *key) {
  return ptr_to_int(hsh_get(ht, key));
}


/* Delete entry with specified key.  
   Returns 1 if item found and deleted, 0 if item not found */
int hsh_delete(Hashtable* ht, const char *key) {
  unsigned int bucket;
  int idx;
  bucket = hsh_hash_func(ht, key);
  if (ht->keys[bucket] == NULL || 
      (idx = lst_find_compare(ht->keys[bucket], (void*)key, equal)) == -1) 
    return 0;
  lst_delete_idx(ht->keys[bucket], idx);
  lst_delete_idx(ht->vals[bucket], idx);
  return 1;
}

/* reset value for given key; returns 0 on success, 1 if item isn't found */
int hsh_reset(Hashtable *ht, const char* key, void* val) {
  unsigned int bucket;
  int idx;
  bucket = hsh_hash_func(ht, key);
  if (ht->keys[bucket] == NULL || 
      (idx = lst_find_compare(ht->keys[bucket], (void*)key, equal)) == -1) 
    return 1;
  
  lst_set_ptr(ht->vals[bucket], idx, val);
  return 0;
}

int hsh_reset_int(Hashtable *ht, const char *key, int val) {
  return hsh_reset(ht, key, int_to_ptr(val));
}


/* Free all resources; does *not* free memory associated with values */
void hsh_free(Hashtable *ht) {
  int i, j;
  for (i = 0; i < ht->nbuckets; i++) {
    if (ht->keys[i] != NULL) {
      for (j = 0; j < lst_size(ht->keys[i]); j++)
        sfree(lst_get_ptr(ht->keys[i], j));
      lst_free(ht->keys[i]);
      lst_free(ht->vals[i]);      
    }
  }
  sfree(ht->keys);
  sfree(ht->vals);
  sfree(ht);
}

/* Free all resources; *does* free memory associated with values */
void hsh_free_with_vals(Hashtable *ht) {
  int i, j;
  for (i = 0; i < ht->nbuckets; i++) {
    if (ht->keys[i] != NULL) {
      for (j = 0; j < lst_size(ht->keys[i]); j++)
        sfree(lst_get_ptr(ht->keys[i], j));
      for (j = 0; j < lst_size(ht->vals[i]); j++)
        sfree(lst_get_ptr(ht->vals[i], j));
      lst_free(ht->keys[i]);
      lst_free(ht->vals[i]);      
    }
  }
  sfree(ht->keys);
  sfree(ht->vals);
  sfree(ht);
}

List *hsh_keys(Hashtable *ht) {
  int i, j;
  List *retval = lst_new_ptr(ht->nbuckets); /* often under-loaded */
  for (i = 0; i < ht->nbuckets; i++) {
    if (ht->keys[i] == NULL) continue;
    for (j = 0; j < lst_size(ht->keys[i]); j++)
      lst_push_ptr(retval, lst_get_ptr(ht->keys[i], j));
  }
  return retval;
}

/* Clear keys and values in a hashtable without freeing the hashtable. The end
   result is equivaslent to a newly-allocated hashtable. */
void hsh_clear_with_vals(Hashtable *ht) {
  int i, j;
  for (i = 0; i < ht->nbuckets; i++) {
    if (ht->keys[i] != NULL) {
      for (j = 0; j < lst_size(ht->keys[i]); j++)
        sfree(lst_get_ptr(ht->keys[i], j));
      for (j = 0; j < lst_size(ht->vals[i]); j++)
        sfree(lst_get_ptr(ht->vals[i], j));
      lst_free(ht->keys[i]);
      lst_free(ht->vals[i]);      
      ht->keys[i] = ht->vals[i] = NULL;
    }
  }
}

/* Clear keys in a hashtable without freeing the hashtable or values. The end
   result is equivaslent to a newly-allocated hashtable, but objects pointed
   to by the hash are left intact. */
void hsh_clear(Hashtable *ht) {
  int i, j;
  for (i = 0; i < ht->nbuckets; i++) {
    if (ht->keys[i] != NULL) {
      for (j = 0; j < lst_size(ht->keys[i]); j++)
	sfree(lst_get_ptr(ht->keys[i], j));
      lst_free(ht->keys[i]);
      lst_free(ht->vals[i]);
      ht->keys[i] = ht->vals[i] = NULL;
    }
  }
}
