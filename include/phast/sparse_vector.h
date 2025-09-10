/* Simple implementation of a sparse vector, stored as a list of pairs
   (idx, val), sorted by idx.  Optimized for set and get rather than
   multiplication or addition. */

#ifndef SPARSE_VECTOR_H
#define SPARSE_VECTOR_H

#include <stdlib.h>

typedef struct {
  List *elementlist;
  unsigned int sorted;
  int nnonzero;
  int dim;
  int refcnt; /* reference count for shallow copy */
} SparseVector;

typedef struct {
  int idx;
  double val;
} SparseVectorElement;

SparseVector *spvec_new(int dim, int starting_size);

void spvec_free(SparseVector *svec);

void spvec_zero(SparseVector *svec);

void spvec_copy(SparseVector *dest, SparseVector *src);

void spvec_set(SparseVector *svec, int idx, double val);

void spvec_set_sorted(SparseVector *svec, int idx, double val);

void spvec_set_lazy(SparseVector *svec, int idx, double val);

int spvec_compare_asc(const void* ptr1, const void* ptr2);

void spvec_sort_by_idx(SparseVector *svec);

double spvec_get(SparseVector *svec, int idx);

SparseVectorElement *spvec_bsearch(SparseVector *svec, int idx);

SparseVectorElement* spvec_linsearch(SparseVector *svec, int idx);

unsigned int spvec_bsearch_idx(SparseVector *svec, int idx, int *lidx);

/* these are to support fast copy of sparse matrices */
static inline void spvec_retain(SparseVector *v) { v->refcnt++; }
static inline void spvec_release(SparseVector *v){
  if (--v->refcnt == 0) spvec_free(v);
}

/* fast copy using memcpy */
static inline void spvec_copy_fast(SparseVector *dest, const SparseVector *src) {
  assert(dest->dim == src->dim);
  lst_cpy_fast(dest->elementlist, src->elementlist);  
  dest->nnonzero = src->nnonzero;
  dest->sorted   = src->sorted;  
}

/* fast clone */
static inline SparseVector* spvec_clone_fast(const SparseVector *src) {
  SparseVector *v = spvec_new(src->dim, lst_size(src->elementlist));
  spvec_copy_fast(v, src);
  return v;
}

#endif
