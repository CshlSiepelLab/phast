/* hashtable - Fast, simple array-based hash table, optimized for 'put' and 'get' ('delete' somewhat inefficient). */

/* $Id: hashtable.h,v 1.1.1.1 2004-06-03 22:43:11 acs Exp $ 
   Written by Adam Siepel, 2002.
   Copyright 2002, Adam Siepel, University of California.
*/


#ifndef HASHTABLE_H
#define HASHTABLE_H

#include "lists.h"

typedef struct hash_table Hashtable;
struct hash_table {
  int nbuckets;                 /* Number of 'buckets' */
  List **keys,                  /* List of char* keys */
    **vals;                     /* Corresponding list of void* values */
};

/* Create new hashtable with initial capacity as specified (in number
   of items).  
   Returns new hashtable with initial capacity as specified. */
Hashtable* hsh_new(int est_capacity);

/* Insert object with specified key and value. */
void hsh_put(Hashtable *ht, char* key, void* val);

/* Retrieve object associated with specified key.
   Returns pointer to object or -1 if not found. 

   Warning: Convention of returning -1 when object is not found is
   inappropriate when objects are integers (needs to be fixed).*/
void* hsh_get(Hashtable* ht, char *key);

/* Delete entry with specified key.  
   Returns 1 if item found and deleted, 0 if item not found */
int hsh_delete(Hashtable* ht, char *key);

/* Free all resources; does *not* free memory associated with values */
void hsh_free(Hashtable *ht);

/* Free all resources; *does* free memory associated with values */
void hsh_free_with_vals(Hashtable *ht);

unsigned int hsh_hash_func(Hashtable *ht, char* key);

#endif
