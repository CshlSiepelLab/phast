/* simple implementation of a sparse matrix, stored as an array of
   sparse vectors, each of which represents a row of the matrix.
   Optimized for set and get rather than multiplication or
   addition. */


#include <stdlib.h>
#include <assert.h>
#include <phast/lists.h>
#include <phast/sparse_matrix.h>
#include <phast/sparse_vector.h>

/* set up a sparse matrix with nrows rows and ncols cols.  The
   parameter colsize is the expected number of nonzero entries per row
   (will be resized as needed) */
SparseMatrix *spmat_new(int nrows, int ncols, int colsize) {
  SparseMatrix *sm = malloc(sizeof(SparseMatrix));
  sm->nrows = nrows;
  sm->ncols = ncols;
  sm->rows = malloc(nrows * sizeof(SparseVector*));
  for (int i = 0; i < nrows; i++)
    sm->rows[i] = spvec_new(ncols, colsize);
  return sm;
}

void spmat_free(SparseMatrix *sm) {
  for (int i = 0; i < sm->nrows; i++)
    spvec_free(sm->rows[i]);
  free(sm->rows);
  free(sm);
}


void spmat_zero(SparseMatrix *sm) {
  for (int i = 0; i < sm->nrows; i++)
    spvec_zero(sm->rows[i]);
}

void spmat_copy(SparseMatrix *dest, SparseMatrix *src) {
  assert(dest->nrows == src->nrows && dest->ncols == src->ncols);
  for (int i = 0; i < src->nrows; i++) 
    spvec_copy(dest->rows[i], src->rows[i]);
}

void spmat_set(SparseMatrix *sm, int row, int col, double val) {
  assert(row >= 0 && row < sm->nrows && col >= 0 && col < sm->ncols);
  spvec_set(sm->rows[row], col, val);
}

void spmat_set_lazy(SparseMatrix *sm, int row, int col, double val) {
  assert(row >= 0 && row < sm->nrows && col >= 0 && col < sm->ncols);
  spvec_set_lazy(sm->rows[row], col, val);
}

void spmat_set_sorted(SparseMatrix *sm, int row, int col, double val) {
  assert(row >= 0 && row < sm->nrows && col >= 0 && col < sm->ncols);
  spvec_set_sorted(sm->rows[row], col, val);
}
  
double spmat_get(SparseMatrix *sm, int row, int col) {
  assert(row >= 0 && row < sm->nrows && col >= 0 && col < sm->ncols);
  return spvec_get(sm->rows[row], col);
}

