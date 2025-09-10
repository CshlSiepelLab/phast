/* simple implementation of a sparse matrix, stored as an array of
   sparse vectors, each of which represents a row of the matrix.
   Optimized for set and get rather than multiplication or
   addition. */

#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <stdlib.h>
#include <phast/sparse_vector.h>
#include <phast/lists.h>

typedef struct {
  int nrows;
  int ncols;
  SparseVector **rows;
} SparseMatrix;

/* these are for fast copying */
static inline void spmat_copy_fast(SparseMatrix *dest, const SparseMatrix *src) {
  assert(dest->nrows == src->nrows && dest->ncols == src->ncols);
  for (int r = 0; r < src->nrows; r++)
    spvec_copy_fast(dest->rows[r], src->rows[r]);
}

/* ensure row is unique before modifying */
static inline void spmat_ensure_unique_row(SparseMatrix *dest, int r) {
  SparseVector *row = dest->rows[r];

  if (dest->rows[r]->refcnt > 1) {
    SparseVector *clone = spvec_clone_fast(row);
    spvec_release(row);
    dest->rows[r] = clone; 
  }
}

static inline void spmat_replace_row_empty(SparseMatrix *M, int r, int capacity_hint) {
  SparseVector *old = M->rows[r];
  SparseVector *new = spvec_new(M->ncols, capacity_hint > 0 ? capacity_hint : 8);
  M->rows[r] = new;
  spvec_release(old);
}

SparseMatrix *spmat_new(int nrows, int ncols, int colsize);

void spmat_free(SparseMatrix *sm);

void spmat_zero(SparseMatrix *sm);

void spmat_copy(SparseMatrix *dest, SparseMatrix *src);

void spmat_set(SparseMatrix *sm, int row, int col, double val);

void spmat_set_lazy(SparseMatrix *sm, int row, int col, double val);

void spmat_set_sorted(SparseMatrix *sm, int row, int col, double val);

double spmat_get(SparseMatrix *sm, int row, int col);

void spmat_copy_shallow(SparseMatrix *dest, const SparseMatrix *src);

#endif
