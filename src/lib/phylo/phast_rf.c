/* calculation of robinson foulds distances */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include "phast/stacks.h"
#include "phast/trees.h"
#include "phast/misc.h"
#include "phast/stringsplus.h"
#include "phast/hashtable.h"
#include "phast/rf.h"

static inline int popcount64(uint64_t x) { return __builtin_popcountll(x); }

/* allocate zeroed mask with W words */
static BitMask *bm_new(int W) {
  BitMask *m = smalloc(sizeof(BitMask));
  m->W = W;
  m->w = calloc(W, sizeof(uint64_t));
  return m;
}

static void bm_free(BitMask *m) { if (!m) return; free(m->w); sfree(m); }

/* set a single bit (0-based global bit index) */
static inline void bm_set(BitMask *m, int bit) {
  int wi = bit >> 6, bi = bit & 63;
  m->w[wi] |= (uint64_t)1ULL << bi;
}

/* dst = a OR b (dst may alias a or b) */
/* static inline void bm_or(BitMask *dst, const BitMask *a, const BitMask *b) { */
/*   for (int i = 0; i < dst->W; ++i) dst->w[i] = a->w[i] | b->w[i]; */
/* } */

/* dst = ~src; then clear high unused bits above nbits */
static inline void bm_not(BitMask *dst, const BitMask *src, int nbits) {
  for (int i = 0; i < dst->W; ++i) dst->w[i] = ~src->w[i];
  int rem = nbits & 63;
  if (rem) {
    uint64_t keep = (rem == 64) ? ~0ULL : ((1ULL << rem) - 1ULL);
    dst->w[dst->W - 1] &= keep;
  }
}

/* count set bits */
static inline int bm_popcount(const BitMask *m) {
  int s = 0; for (int i = 0; i < m->W; ++i) s += popcount64(m->w[i]); return s;
}

/* lexicographic compare (lowest word first is fine as long as consistent) */
static int bm_cmp_words(const void *pa, const void *pb) {
  const BitMask *a = *(BitMask* const*)pa;
  const BitMask *b = *(BitMask* const*)pb;
  /* compare from high word to low word for nice order */
  for (int i = a->W - 1; i >= 0; --i) {
    if (a->w[i] < b->w[i]) return -1;
    if (a->w[i] > b->w[i]) return 1;
  }
  return 0;
}

/* deep copy */
static BitMask *bm_clone(const BitMask *m) {
  BitMask *c = bm_new(m->W);
  memcpy(c->w, m->w, sizeof(uint64_t)*m->W);
  return c;
}

/* canonicalize split (make it the smaller side). Returns NEW mask on heap. */
static BitMask *bm_canonical(const BitMask *m, int nbits) {
  int sz = bm_popcount(m);
  int other = nbits - sz;
  if (sz == 0 || other == 0) return NULL;            /* trivial */
  if (sz == 1 || other == 1) return NULL;            /* leaf edge; ignore */

  if (sz <= other) {
    return bm_clone(m);
  } else {
    BitMask *c = bm_new(m->W);
    bm_not(c, m, nbits);
    return c;
  }
}

static void mv_init(MaskVec *mv, int cap) {
  mv->cap = cap > 0 ? cap : 16; mv->size = 0;
  mv->a = smalloc(sizeof(BitMask*) * mv->cap);
}
static void mv_push(MaskVec *mv, BitMask *m) {
  if (m == NULL) return; /* trivial split filtered */
  if (mv->size == mv->cap) {
    mv->cap *= 2;
    mv->a = srealloc(mv->a, sizeof(BitMask*) * mv->cap);
  }
  mv->a[mv->size++] = m;
}
static void mv_free(MaskVec *mv) {
  for (int i = 0; i < mv->size; ++i) bm_free(mv->a[i]);
  sfree(mv->a); mv->a = NULL; mv->size = mv->cap = 0;
}

/* name lookup (sorted common list) */
/* static int name_to_index_in_sorted(List *sorted_names, char *name) { */
/*   int lo = 0, hi = lst_size(sorted_names) - 1; */
/*   while (lo <= hi) { */
/*     int mid = (lo + hi) >> 1; */
/*     String *midS = lst_get_ptr(sorted_names, mid); */
/*     int cmp = strcmp(name, midS->chars); */
/*     if (cmp == 0) return mid; */
/*     if (cmp < 0) hi = mid - 1; else lo = mid + 1; */
/*   } */
/*   return -1; */
/* } */

static inline int name_to_index(Hashtable *ht, const char *name) {
  void *vp = hsh_get(ht, name);
  if (vp == (void*)-1) return -1;
  return ptr_to_int(vp);
}

/* DFS over TreeNode to collect edge splits */
/* Returns heap mask for subtree of u. Parent is used to avoid back-edge in unrooted reps. */
static BitMask *dfs_collect_splits(TreeNode *u, TreeNode *parent,
                                   Hashtable *name2idx, int n, int W,
                                   MaskVec *splits) {
  if (u->lchild == NULL && u->rchild == NULL) {
    /* leaf */
    int idx = name_to_index(name2idx, u->name);
    if (idx < 0) die("Leaf '%s' not in name list.\n", u->name);
    BitMask *m = bm_new(W);
    bm_set(m, idx);
    return m;
  }

  BitMask *mask_u = bm_new(W);

  /* left child */
  if (u->lchild && u->lchild != parent) {
    BitMask *ml = dfs_collect_splits(u->lchild, u, name2idx, n, W, splits);
    /* edge split defined by (u—lchild) => use ml as one side */
    BitMask *canL = bm_canonical(ml, n);
    mv_push(splits, canL);
    /* accumulate */
    for (int i = 0; i < W; ++i) mask_u->w[i] |= ml->w[i];
    bm_free(ml);
  }
  /* right child */
  if (u->rchild && u->rchild != parent) {
    BitMask *mr = dfs_collect_splits(u->rchild, u, name2idx, n, W, splits);
    BitMask *canR = bm_canonical(mr, n);
    mv_push(splits, canR);
    for (int i = 0; i < W; ++i) mask_u->w[i] |= mr->w[i];
    bm_free(mr);
  }

  /* If your TreeNode can have >2 children (multifurcating),
     iterate a generic adjacency list here; the logic is the same:
     for each child c!=parent, collect mc, push canonical(mc), OR into mask_u. */

  return mask_u;
}

/* calculate the symmetric Robinson Foulds distance between two
   trees. Considers topology only.  Leaf names must match exactly.
   This is an O(n log n) implementation  */
double tr_robinson_foulds(TreeNode *t1, TreeNode *t2) {
  /* build sorted lists of leaf names */
  List *names1 = tr_leaf_names(t1);
  List *names2 = tr_leaf_names(t2);
  lst_qsort_str(names1, ASCENDING);
  lst_qsort_str(names2, ASCENDING);
  if (str_list_equal(names1, names2) == FALSE)
    die("ERROR in tr_robinson_foulds: trees do not have matching leaf names.\n");

  int n = lst_size(names1);
  if (n < 3) { lst_free_strings(names1); lst_free_strings(names2); return 0; }

  /* build name2idx map from names1 */
  Hashtable *name2idx = hsh_new(2 * lst_size(names1)); 
  for (int i = 0; i < lst_size(names1); i++) {
    String *s = lst_get_ptr(names1, i);
    hsh_put_int(name2idx, s->chars, i);
  }
  
  const int W = (n + 63) >> 6; /* words needed */

  /* collect canonical splits from each tree */
  MaskVec S1, S2;
  mv_init(&S1, n); mv_init(&S2, n);

  BitMask *rootmask1 = dfs_collect_splits(t1, NULL, name2idx, n, W, &S1);
  BitMask *rootmask2 = dfs_collect_splits(t2, NULL, name2idx, n, W, &S2);
  bm_free(rootmask1); bm_free(rootmask2);

  /* sort both split multisets (lexicographic) */
  qsort(S1.a, S1.size, sizeof(BitMask*), bm_cmp_words);
  qsort(S2.a, S2.size, sizeof(BitMask*), bm_cmp_words);

  /* count common splits (intersection) with two pointers */
  int i = 0, j = 0, common = 0;
  while (i < S1.size && j < S2.size) {
    int cmp = bm_cmp_words(&S1.a[i], &S2.a[j]);
    if (cmp == 0) { common++; i++; j++; }
    else if (cmp < 0) i++;
    else j++;
  }

  /* symmetric RF distance = |S1| + |S2| - 2|common| */
  int RF = (S1.size + S2.size) - 2*common;

  /* cleanup */
  mv_free(&S1); mv_free(&S2);
  lst_free_strings(names2);
  lst_free_strings(names1);
  hsh_free(name2idx);

  return RF;
}

