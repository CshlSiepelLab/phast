/* Simple implementation of a sparse vector, stored as a list of pairs
   (idx, val), sorted by idx.  Optimized for set and get rather than
   multiplication or addition. */

/* FIXME: store the element structs contiguously rather than via
   pointers for better cache behavior */

#include <stdlib.h>
#include <assert.h>
#include <phast/lists.h>
#include <phast/sparse_vector.h>
#include <phast/sparse_matrix.h>

SparseVector *spvec_new(int dim, int starting_size) {
  SparseVector *svec = malloc(sizeof(SparseVector));
  svec->elementlist = lst_new_ptr(starting_size);
  svec->sorted = TRUE;
  svec->nnonzero = 0;
  svec->dim = dim;
  return svec;
}

void spvec_free(SparseVector *svec) {        
  spvec_zero(svec);
  lst_free(svec->elementlist);
  free(svec);
}                                            

SparseVectorElement *spvec_new_el(int idx, double val) {
  SparseVectorElement *el = malloc(sizeof(SparseVectorElement));
  el->idx = idx;
  el->val = val;
  return el;
}

void spvec_zero(SparseVector *svec) {
  for (int i = 0; i < lst_size(svec->elementlist); i++)
    free(lst_get_ptr(svec->elementlist, i));
  lst_clear(svec->elementlist);
  svec->nnonzero = 0;
}

void spvec_copy(SparseVector *dest, SparseVector *src) {
  assert(dest->dim == src->dim);
  if (dest->nnonzero == 0 && src->nnonzero == 0)
    return; /* shortcut if both empty */
  dest->sorted = src->sorted;
  dest->nnonzero = src->nnonzero;
  spvec_zero(dest);
  for (int i = 0; i < lst_size(src->elementlist); i++) {
    SparseVectorElement *oldel = lst_get_ptr(src->elementlist, i);
    SparseVectorElement *newel = spvec_new_el(oldel->idx, oldel->val);
    lst_push_ptr(dest->elementlist, newel);
  }
}

/* general set; maintains ascending order by idx but takes O(log n),
   where n is svec->nnonzero */
void spvec_set(SparseVector *svec, int idx, double val) {
  int lidx = -1;
  SparseVectorElement *el = NULL;
  assert(idx >= 0 && idx < svec->dim);
  if (val == 0)
    return;
  
  /* find the insertion point */
  if (spvec_bsearch_idx(svec, idx, &lidx) == TRUE) {
    el = lst_get_ptr(svec->elementlist, lidx); /* found; update val */
    el->val = val;
  }
  else {   /* otherwise insert it at the correct place */
    el = spvec_new_el(idx, val);
    lst_insert_idx(svec->elementlist, lidx, el);
  }
  svec->nnonzero++;
}

/* O(1) set for use when elements can be added in sorted order by
   idx; often possible when working from a loop */
void spvec_set_sorted(SparseVector *svec, int idx, double val) {
  assert(idx >= 0 && idx < svec->dim);
  if (val == 0)
    return;
  
  lst_push_ptr(svec->elementlist, spvec_new_el(idx, val));
  svec->nnonzero++;
}

/* O(1) set that allows sorting to be relaxed until retrieval.
   Calling code must avoid setting multiple values with the same
   index */
void spvec_set_lazy(SparseVector *svec, int idx, double val) {
  assert(idx >= 0 && idx < svec->dim);
  if (val == 0)
    return;
  lst_push_ptr(svec->elementlist, spvec_new_el(idx, val));
  /* note: there may be dups; have to resolve later */
  svec->nnonzero++;
  svec->sorted = FALSE;
}

int spvec_compare_asc(const void* ptr1, const void* ptr2) {
  SparseVectorElement *el1 = (SparseVectorElement*)ptr1;
  SparseVectorElement *el2 = (SparseVectorElement*)ptr2;
  if (el1->idx == el2->idx) return 0;  /* shouldn't happen */
  else if (el1->idx < el2->idx) return -1;
  return 1;                    
}

/* sort elements by index */
void spvec_sort_by_idx(SparseVector *svec) {
  lst_qsort(svec->elementlist, spvec_compare_asc);
}

/* returns val if found or zero otherwise.  Takes O(log n) */
double spvec_get(SparseVector *svec, int idx) {
  SparseVectorElement *el;
  assert(idx >= 0 && idx < svec->dim);
  if (svec->sorted == FALSE)
    spvec_sort_by_idx(svec);  
  el = spvec_bsearch(svec, idx);
  if (el == NULL)
    return 0;
  return el->val;
}

/* Binary search for an index, returns element if found or NULL otherwise */
SparseVectorElement *spvec_bsearch(SparseVector *svec, int idx) {
  SparseVectorElement *match = NULL;
  int lidx;
  
  if (spvec_bsearch_idx(svec, idx, &lidx))
    match = lst_get_ptr(svec->elementlist, lidx);

  return match;
}

/* Binary search for a given index, finds its list index if
   it exists, or if it does not, the list index of the largest val less
   than the target index (the one after which it should fall if inserted).  Return
   value is TRUE if the index does exist in the list and FALSE
   otherwise.  If the query is smaller than all values in the list
   *lidx will be set to -1 and FALSE will be returned */
unsigned int spvec_bsearch_idx(SparseVector *svec, int idx, int *lidx) {
  SparseVectorElement *candidate;
  int l = 0;
  int r = lst_size(svec->elementlist) - 1;

  if (r < 0 || idx < ((SparseVectorElement*)lst_get_ptr(svec->elementlist, 0))->idx) {
    *lidx = -1;
    return FALSE;
  }
  if (idx > ((SparseVectorElement*)lst_get_ptr(svec->elementlist, r))->idx) {
    *lidx = r;
    return FALSE;
  }
  
  while (l <= r) {
    int m = (l + r)/2;
    candidate = lst_get_ptr(svec->elementlist, m);
    if (idx == candidate->idx) {
      *lidx = m;
      return TRUE;
    }
    else if (idx < candidate->idx) 
      r = m - 1;
    else                        /* idx > candidate->idx */
      l = m + 1;
  }

  assert(l == r+1);  /* has to be true on exit from loop */
  
  /* the index must fall between r and r+1 (or equivalently, l-1 and l) */
  *lidx = r;
  return FALSE;
}

