/* bitset.h — generic N-bit sets + hash map */

#ifndef BITSET_H
#define BITSET_H

#include <stdint.h>
#include "phast/lists.h"

/* ===== BSet ===== */
typedef struct {
  int nbits;        /* total bits */
  int nwords;       /* number of uint64_t words */
  uint64_t *w;      /* words[0..nwords-1] */
} BSet;

/* Construction / basics */
BSet *bs_new(int nbits);
BSet *bs_clone(const BSet *src);
void    bs_copy(BSet *dst, const BSet *src);   /* same nbits */
void    bs_free(BSet *bs);
void    bs_zero(BSet *bs);
int     bs_equals(const BSet *a, const BSet *b);

/* Bit operations */
void bs_set_bit(BSet *bs, int idx);
void bs_or(BSet *dst, const BSet *a, const BSet *b);   /* dst = a | b */
void bs_and(BSet *dst, const BSet *a, const BSet *b);  /* dst = a & b */
void bs_flip(BSet *bs);                                    /* invert all bits within nbits */

/* List helpers (using your List* with int* elements) */
void  bs_set_bits_from_list(BSet *bs, List *idxs);   /* set bits in idxs */
List *bs_to_index_list(const BSet *bs);                    /* returns List* of int* */

/* Hash of raw words (FNV-1a 64-bit) */
uint64_t bs_hash64(const BSet *bs);

/* ===== Hash map BSet -> void* (open addressing) ===== */
typedef struct {
  int used;        /* 0 empty, 1 used */
  uint64_t h;      /* cached hash */
  BSet *key;     /* owned copy of key */
  void   *val;     /* user value */
} BSHashEntry;

typedef struct {
  int nslots;      /* power-of-two capacity */
  int nused;       /* number of used entries */
  BSHashEntry *a;  /* slots */
} BSHash;

BSHash *bs_hash_new(int init_capacity);  /* cap will round up to a power of two */
void    bs_hash_free(BSHash *H, int free_vals);
void   *bs_hash_get(const BSHash *H, const BSet *key);
int     bs_hash_put(BSHash *H, const BSet *key, void *value); /* 1 insert, 0 replace */

#endif /* BITSET_H */
