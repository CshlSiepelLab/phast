/* simple implementation of a sparse matrix, stored as an array of
   sparse vectors, each of which represents a row of the matrix.
   Optimized for set and get rather than multiplication or
   addition. */

#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <stdlib.h>
#include <phast/sparse_vector.h>

typedef struct {
  int nrows;
  int ncols;
  SparseVector **rows;
} SparseMatrix;

SparseMatrix *spmat_new(int nrows, int ncols, int colsize);

void spmat_free(SparseMatrix *sm);

void spmat_zero(SparseMatrix *sm);

void spmat_copy(SparseMatrix *dest, SparseMatrix *src);

void spmat_set(SparseMatrix *sm, int row, int col, double val);

void spmat_set_lazy(SparseMatrix *sm, int row, int col, double val);

void spmat_set_sorted(SparseMatrix *sm, int row, int col, double val);

double spmat_get(SparseMatrix *sm, int row, int col);

#endif
