/* bitset.c */

#include <string.h>
#include <assert.h>
#include "phast/misc.h"   
#include "phast/lists.h"  
#include "phast/bitset.h"

/* ===== internals ===== */
static inline int words_for_bits(int nbits) { return (nbits + 63) >> 6; }
static inline uint64_t last_mask(int nbits) {
  int r = (nbits & 63);
  return (r == 0) ? ~0ULL : ((1ULL << r) - 1ULL);
}

/* ===== BSet impl ===== */
BSet *bs_new(int nbits) {
  assert(nbits > 0);
  BSet *bs = (BSet*)smalloc(sizeof(BSet));
  bs->nbits  = nbits;
  bs->nwords = words_for_bits(nbits);
  bs->w      = (uint64_t*)scalloc(bs->nwords, sizeof(uint64_t));
  return bs;
}

BSet *bs_clone(const BSet *src) {
  BSet *bs = bs_new(src->nbits);
  memcpy(bs->w, src->w, src->nwords * sizeof(uint64_t));
  return bs;
}

void bs_copy(BSet *dst, const BSet *src) {
  assert(dst->nbits == src->nbits);
  memcpy(dst->w, src->w, dst->nwords * sizeof(uint64_t));
}

void bs_free(BSet *bs) {
  if (!bs) return;
  sfree(bs->w);
  sfree(bs);
}

void bs_zero(BSet *bs) {
  memset(bs->w, 0, bs->nwords * sizeof(uint64_t));
}

int bs_equals(const BSet *a, const BSet *b) {
  if (a->nbits != b->nbits) return 0;
  for (int i = 0; i < a->nwords; ++i)
    if (a->w[i] != b->w[i]) return 0;
  return 1;
}

void bs_set_bit(BSet *bs, int idx) {
  assert(idx >= 0 && idx < bs->nbits);
  bs->w[idx >> 6] |= (1ULL << (idx & 63));
}

void bs_or(BSet *dst, const BSet *a, const BSet *b) {
  assert(dst->nbits == a->nbits && a->nbits == b->nbits);
  for (int i = 0; i < dst->nwords; ++i) dst->w[i] = a->w[i] | b->w[i];
}

void bs_and(BSet *dst, const BSet *a, const BSet *b) {
  assert(dst->nbits == a->nbits && a->nbits == b->nbits);
  for (int i = 0; i < dst->nwords; ++i) dst->w[i] = a->w[i] & b->w[i];
}

void bs_flip(BSet *bs) {
  for (int i = 0; i < bs->nwords; ++i) bs->w[i] = ~bs->w[i];
  /* clear padding bits beyond nbits in the last word */
  bs->w[bs->nwords - 1] &= last_mask(bs->nbits);
}

/* ===== List helpers ===== */
/* Assumes 'idxs' is a List* of int* (heap-allocated ints).
   Bits out of range [0..nbits-1] are ignored with an assert. */
void bs_set_bits_from_list(BSet *bs, List *idxs) {
  int n = lst_size(idxs);
  for (int i = 0; i < n; ++i) {
    int idx = lst_get_int((List*)idxs, i);
    assert(idx >= 0 && idx < bs->nbits);
    bs_set_bit(bs, idx);
  }
}

/* Returns a List* of int* (indices of 1-bits), caller owns ints & list */
List *bs_to_index_list(const BSet *bs) {
  List *out = lst_new_ptr(16);
  for (int w = 0; w < bs->nwords; ++w) {
    uint64_t word = bs->w[w];
    while (word) {
      int t = __builtin_ctzll(word);  /* index of lowest 1 bit */
      int idx = (w << 6) + t;
      if (idx < bs->nbits) 
        lst_push_int(out, idx);
      word &= (word - 1);             /* clear lowest 1 bit */
    }
  }
  return out;
}

/* ===== Hash (FNV-1a) ===== */
uint64_t bs_hash64(const BSet *bs) {
  uint64_t h = 1469598103934665603ULL;
  const uint8_t *bytes = (const uint8_t*)bs->w;
  size_t nbytes = (size_t)bs->nwords * sizeof(uint64_t);
  for (size_t i = 0; i < nbytes; ++i) {
    h ^= bytes[i];
    h *= 1099511628211ULL;
  }
  return h;
}

/* ===== BSet -> void* hash map (open addressing) ===== */

static inline int is_power_of_two(int x) { return (x & (x-1)) == 0; }
static int round_up_pow2(int x) {
  int p = 1; while (p < x) p <<= 1; return p;
}
static inline int key_eq(const BSet *a, const BSet *b) {
  return bs_equals(a, b);
}

BSHash *bs_hash_new(int init_capacity) {
  if (init_capacity < 16) init_capacity = 16;
  if (!is_power_of_two(init_capacity)) init_capacity = round_up_pow2(init_capacity);
  BSHash *H = (BSHash*)smalloc(sizeof(BSHash));
  H->nslots = init_capacity;
  H->nused  = 0;
  H->a      = (BSHashEntry*)scalloc(H->nslots, sizeof(BSHashEntry));
  return H;
}

static void bs_hash_resize(BSHash *H, int new_cap) {
  BSHashEntry *old = H->a; int oldN = H->nslots;
  H->nslots = round_up_pow2(new_cap);
  H->a = (BSHashEntry*)scalloc(H->nslots, sizeof(BSHashEntry));
  H->nused = 0;
  for (int i = 0; i < oldN; ++i) if (old[i].used) {
    /* reinsert */
    uint64_t h = old[i].h;
    int mask = H->nslots - 1;
    int j = (int)(h & mask);
    while (H->a[j].used) j = (j + 1) & mask;
    H->a[j] = old[i];
    H->nused++;
  }
  sfree(old);
}

void bs_hash_free(BSHash *H, int free_vals) {
  if (!H) return;
  for (int i = 0; i < H->nslots; ++i) {
    if (H->a[i].used) {
      bs_free(H->a[i].key);
      if (free_vals && H->a[i].val) sfree(H->a[i].val);
    }
  }
  sfree(H->a);
  sfree(H);
}

void *bs_hash_get(const BSHash *H, const BSet *key) {
  if (!H) return NULL;
  uint64_t h = bs_hash64(key);
  int mask = H->nslots - 1;
  int j = (int)(h & mask);
  while (H->a[j].used) {
    if (H->a[j].h == h && key_eq(H->a[j].key, key))
      return H->a[j].val;
    j = (j + 1) & mask;
  }
  return NULL;
}

int bs_hash_put(BSHash *H, const BSet *key, void *value) {
  /* resize if load > 0.7 */
  if ((H->nused + 1) * 10 > H->nslots * 7) bs_hash_resize(H, H->nslots << 1);

  uint64_t h = bs_hash64(key);
  int mask = H->nslots - 1;
  int j = (int)(h & mask);
  while (H->a[j].used) {
    if (H->a[j].h == h && key_eq(H->a[j].key, key)) {
      /* replace */
      H->a[j].val = value;
      return 0; /* replaced */
    }
    j = (j + 1) & mask;
  }
  /* insert fresh: make an owned copy of key */
  BSet *copy = bs_clone(key);
  H->a[j].used = 1;
  H->a[j].h    = h;
  H->a[j].key  = copy;
  H->a[j].val  = value;
  H->nused++;
  return 1; /* inserted */
}
