/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file hashtable.h
 Fast, simple array-based hash table, optimized for 'put' and 'get' ('delete' somewhat inefficient). 
  @ingroup base
*/


#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <lists.h>
#include <misc.h>
#include <external_libs.h>

/** Recommended by Kernighan and Pike */
#define MULTIPLIER 31
#define LOADING_FACTOR 5

typedef struct hash_table Hashtable;
/** Hash table struct  */
struct hash_table {
  int nbuckets;                 /**< Number of 'buckets' */
  List **keys,                  /**< List of char* keys */
    **vals;                     /**< Corresponding list of void* values */
};

/** \name HashTable allocation functions 
 \{ */

/** Create new hashtable.  
   @param est_capacity Estimated needed capacity (in number of items)
   @return New hashtable with initial capacity as specified.
 */
Hashtable* hsh_new(int est_capacity);

/**  Copies a hash table.
  @param ht Hash table to copy from
  @result Clone of original hash table
  @warning If values are pointers, only pointers will be copied, not values being pointed to.
  @note Does copy keys
*/

Hashtable* hsh_copy(Hashtable *ht);

/** \} \name HashTable misc. functions
 \{ */

/* we'll only inline the functions likely to be used heavily in inner
   loops */  

/** Hashing function mapping key to index of array holding values 
   @param ht Hash Table to calculate mapping for 
   @param key Key to get index in hash table for
   @result Index in array holding values associated with key
*/
static PHAST_INLINE
unsigned int hsh_hash_func(Hashtable *ht, const char* key) {
  unsigned int h = 0;
  int i = 0;
  for (i = 0; key[i] != '\0'; i++)
    h = MULTIPLIER * h + key[i];
  return h % ht->nbuckets;
}

/** Make a list of all the keys in the hash table.
  @param ht Hash Table to list keys for
  @result List of all keys in hash table ht
*/
List *hsh_keys(Hashtable *ht);


/** \} \name HashTable put functions 
 \{ */


/** Put a new value into hash table referred to by key.
   @param ht Hash Table to add entry to
   @param key Key associated with value so we can retrieve/modify it later
   @param val Value associated with key that we wish to store
*/
static PHAST_INLINE
void hsh_put(Hashtable *ht, const char* key, void* val) {
  unsigned int bucket = hsh_hash_func(ht, key);
  char *keycpy;
  if (ht->keys[bucket] == NULL) {
    ht->keys[bucket] = lst_new_ptr(LOADING_FACTOR);
    ht->vals[bucket] = lst_new_ptr(LOADING_FACTOR);
  }
  keycpy = (char*)smalloc(sizeof(char) * (strlen(key) + 1));
  strcpy(keycpy, key);

  lst_push_ptr(ht->keys[bucket], keycpy);
  lst_push_ptr(ht->vals[bucket], val);
}

/** Add an integer to the hash table 
  @param ht Hash table to add integer to
  @param key Key to refer to value by
  @param val Data that will be added to hash table
*/
void hsh_put_int(Hashtable *ht, const char *key, int val);

/** \name HashTable get functions 
 \{ */

/* hsh_get involves call to lst_find_compare with function pointer;
   not sure how best to inline, if at all */ 
/** Retrieve object associated with specified key.
  @param ht Hash Table to retrieve value from 
  @param Key key associated with the value to retrieve
  @result Object associated with key, if key is not found -1 returned
  
*/
void* hsh_get(Hashtable* ht, const char *key);

/** Retrieve integer associated with specified key.
  @param ht Hash Table to retrieve integer value from 
  @param key Key associated with the integer value to retrieve
  @result Integer associated with key, if key is not found -1 returned
*/
int hsh_get_int(Hashtable *ht, const char *key);

/** \name HashTable remove functions 
 \{ */

/** Delete entry from hash table with specified key.
  @param Hash Table to delete entry from
  @param key Key associated with entry in hash table to delete
  @result 1 if successful, 0 if item not found
*/
int hsh_delete(Hashtable* ht, const char *key);

/** Reset value for given key.
  @param ht Hash Table containing element to modify
  @param key Key associated with entry
  @param val Value to be associated with key
  @result 0 if successful, 1 if item not found
*/
int hsh_reset(Hashtable *ht, const char* key, void* val);

/** Reset integer value for given key.
  @param ht Hash Table containing integer to modify
  @param key Key associated with integer value
  @param val Integer value to be associated with key
  @result 0 if successful, 1 if item not found
*/
int hsh_reset_int(Hashtable *ht, const char *key, int val);

/** \name HashTable cleanup functions 
 \{ */

/** Free hash table but not values
  @param ht Hash Table to free
  @warning Does *NOT* free memory associated with hash table values
*/
void hsh_free(Hashtable *ht);

/** Free hash table and values
  @param ht Hash Table to free
  @warning Also frees memory associated with hash table values
*/
void hsh_free_with_vals(Hashtable *ht);


/** Clear keys but not values in a hash table.
  @param ht Hash Table to clear out
  @note The end result is equivalent to a newly-allocated hash table, but objects pointed to by the hash are left intact
*/
void hsh_clear(Hashtable *ht);

/** Clear keys and values in a hash table without freeing the hash table.
  @param ht Hash Table to clear out
  @note The end result is equivalent to a newly-allocated hash table
*/
void hsh_clear_with_vals(Hashtable *ht);

/** \} */

#endif
