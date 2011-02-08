/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** \file matrix.h
   Structure to hold Matrices of real numbers (doubles) and functions to perform basic matrix operations on them
   \ingroup base
*/

#ifndef MAT_H
#define MAT_H

#include <vector.h>
#include <external_libs.h>

/** Equality threshold -- consider equal if this close */
#define EQ_THRESHOLD 1e-10

/** Matrix structure -- just a 2d array of doubles and its dimensions */
struct matrix_struct {
  double **data;			/**< underlying array of doubles */
  int nrows;			/**< number of rows */
  int ncols;			/**< number of columns */
};
/** Matrix type */
typedef struct matrix_struct Matrix;

/** \name Matrix create/free operations. */

/** \{ */

/** Allocate new matrix.

  @param nrows Number of rows.
  @param ncols Number of columns.

  @see mat_new_from_array.
*/
Matrix *mat_new(int nrows, int ncols);
/** Allocate new matrix using supplied array as storage.

  @param array Pre allocated storage array.
  @param nrows Number of rows.
  @param ncols Number of columns.
  
  @see mat_new, mat_free.
*/
Matrix *mat_new_from_array(double **array, int nrows, int ncols);

/** Allocate matrix, reading entries from file.

  Entries are read row by row.

  @param F Input file stream.
  @param nrows Number of rows.
  @param ncols Number of columns.
  
  @see mat_read.
*/
Matrix *mat_new_from_file(FILE *F, int nrows, int ncols);

/** Allocate a new matrix identical to given input matrix.

  @param src Input matrix.
  @result Newly allocated copy of src
  @see mat_copy.
*/
Matrix *mat_create_copy(Matrix *src);

/** Free matrix.

  @note It will release storage array even if matrix was created using
        mat_new_from_array.

  @param m Input matrix.
*/
void mat_free(Matrix *m);

/** Resize matrix.

  Resize will preserve values that fall within the new dimensions.
  Value that fall outside the dimensions of the source matrix are undefined.
  
  @param m Matrix to be resized.
  @param nrows New row count.
  @param ncols New column count.
*/
void mat_resize(Matrix *m, int nrows, int ncols);

/** \} */

/** \name Matrix initialization operations. */

/** \{ */

/** Set matrix entries to match identity matrix.

  @param m Input matrix.
  
  @see mat_set, mat_set_zero, mat_set_all.
*/
void mat_set_identity(Matrix *m);
/** Set all matrix entries to zero.

  @param m Input matrix.
  @see mat_set, mat_set_identity, mat_set_all.
*/
void mat_zero(Matrix *m);
/** Set all matrix entries to the specified value.

  @param m Input matrix.
  @param val Value to be stored.
  @see mat_set, mat_set_identity, mat_set_all.
*/
void mat_set_all(Matrix *m, double val);
/** Copy all matrix entries from source to destination matrix.

  @note Matrices must have the same dimensions.

  @param src Source matrix.
  @param dest Destination matrix.
  
  @see mat_create_copy.
*/
void mat_copy(Matrix *dest, Matrix *src);

/** Read matrix entries from file.

  Entries are read row by row.

  @param m Pre allocated destination matrix.
  @param F Input file stream.
  
  @see mat_new_from_file, mat_print.
*/
void mat_read(Matrix *m, FILE *F);

/** \} */

/** Retrieve value.

  @param m Input matrix.
  @param row Row index.
  @param col Column index.
  @result Element m[row][col]
  @see mat_get_row, mat_get_col
*/
static PHAST_INLINE 
double mat_get(Matrix *m, int row, int col) {
  return m->data[row][col];
}
/** Retrieve matrix row.

  Allocates a new vector and copies row values from input matrix.

  @param m Input matrix.
  @param row Row index.
  @result Row 'row' or m, m[row]
  @see mat_get, mat_get_col.
*/
Vector *mat_get_row(Matrix *m, int row);
/** Retrieve matrix column.

  Allocates a new vector and copies column values from input matrix.
  
  @param m Input matrix.
  @param col Column index.
  @result Column 'col' of m, m[*][col]
  @see mat_get_row, mat_get.
*/
Vector *mat_get_col(Matrix *m, int col);
/** Store value in matrix.

  @param m Input matrix.
  @param row Row index.
  @param col Column index.
  @param val Value to be stored.
  
  @see mat_set_identity,mat_zero, mat_set_all.
*/
static PHAST_INLINE
void mat_set(Matrix *m, int row, int col, double val) {
  m->data[row][col] = val;
}

/** Print matrix to file.

  Entries are separated by spaces. If the minimum value of the matrix
  is less than 1E-3, matrix is printed in exponential notation.

  @param m Input matrix.
  @param F Output file stream.
  
  @see mat_read.
*/
void mat_print(Matrix *m, FILE *F);

/** \name Matrix math operations. */

/** \{ */

/** Compare two matrices.
  @param A 1st matrix.
  @param B 2nd matrix.
  @result Returns 1 if A and B are same dimensions and data in every 
  row and column are the same.  Otherwise returns 0.
 */
int mat_equal(Matrix *A, Matrix *B);

/** Create transposed version of input matrix.

  @param src Input matrix.
  @result Transpose of src matrix
*/
Matrix *mat_transpose(Matrix *src);
/** Rescale matrix entries by given value.

  @param m Input matrix.
  @param scale_factor Value to use as scale factor.
*/
void mat_scale(Matrix *m, double scale_factor);

/** Multiply two matrices.

  \code 
  // equivalent to this, but with matrices
  prod = m1 * m2;
  \endcode

  \note Destination matrix must be different from source matrices.

  @param prod Pre allocated result matrix.
  @param m1 Input matrix one.
  @param m2 Input matrix two.
*/
void mat_mult(Matrix *prod, Matrix *m1, Matrix *m2);
/** Multiply matrix with vector. 

  \code 
  // equivalent to this, but with matrices
  prod = m * v;
  \endcode
  
  @param prod pre allocated result vector.
  @param m Input matrix.
  @param v Input vector.
*/
void mat_vec_mult(Vector *prod, Matrix *m, Vector *v);
/** Matrix increment.

  \code 
  // equivalent to this, but with matrices
  thism += addm;
  \endcode
  
  @param thism Result matrix.
  @param addm Input matrix.

  \sa mat_minus_eq.
*/
void mat_plus_eq(Matrix *thism, Matrix *addm);
/** Matrix decrement.

  \code 
  // equivalent to this, but with matrices
  thism -= subm;
  \endcode
  
  @param thism Result matrix.
  @param subm Input matrix.

  \sa mat_plus_eq.
*/
void mat_minus_eq(Matrix *thism, Matrix *subm);
/** Linear combination of matrices.

  \code 
  // equivalent to this, but with matrices
  dest = coef1 * src1 + coef2 * src2;
  \endcode
  
  @param dest Result matrix.
  @param src1 Input matrix one.
  @param coef1 Coefficient of input matrix one.
  @param src2 Input matrix two.
  @param coef2 Coefficient of input matrix two.

*/
void mat_linear_comb(Matrix *dest, Matrix *src1, double coef1, 
                     Matrix *src2, double coef2);


/** Linear combination of an arbitrary number of matrices.

@code
// equivalent to this, but with matrices
dest=0;
for (i=0; i < n; i++)
  dest += coef[i]*src[i]
@endcode

@param dest Result matrix.
@param n Number of src matrices.
@param src An array of n input matrices.
@param coef A vector of n coefficients.

*/
void mat_linear_comb_many(Matrix *dest, int n, Matrix **src,
			  double *coef);

/** Invert square, real, non-symmetric matrix

  Uses LU decomposition (LAPACK routines dgetrf and dgetri).
  
  @param M_inv Pre allocated destination matrix.
  @param M Input matrix.
  @result 0 on success, 1 on failure
  @warning Calling this is function will terminate program if program was not compiled
  with LAPACK support.
*/
int mat_invert(Matrix *M_inv, Matrix *M);

/** Pre and post multiply a diagonal matrix.

   Compute A = B * C * D where A, B, C, D are square matrices of the
   same dimension, and C is diagonal.  C is described by a vector
   representing its diagonal elements.
   
   @param A Pre allocated result matrix.
   @param B Input matrix B.
   @param C Input vector representing diagonal of matrix C.
   @param D Input matrix C.
*/
void mat_mult_diag(Matrix *A, Matrix *B, Vector *C, Matrix *D);

#ifndef SKIP_LAPACK

/** Convert a matrix into Lapack's column-major format.
    @param m (input) A matrix.
    @param arr (output) A representation of m in Lapack's columns-major 
               format.
*/
void mat_to_lapack(Matrix *m, LAPACK_DOUBLE *arr);

/** Convert back from Lapack's column-major format to a Matrix.
    @param m (output) A matrix
    @param arr (input) A representation of a matrix in Lapack's
               column-major format.
 */
void mat_from_lapack(Matrix *m, LAPACK_DOUBLE *arr);
#endif

/* \} */

#endif
