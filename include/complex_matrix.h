/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file complex_matrix.h
    Matrices of complex numbers
    @ingroup base
*/

#ifndef ZMAT_H
#define ZMAT_H

#include <complex.h>
#include <complex_vector.h>
#include <external_libs.h>

/** Structure for matrix of complex numbers -- 2d array of Complex
    objects and its dimensions */
typedef struct {
  Complex **data;  /**< Contains matrix data as complex numbers*/
  int nrows;    /**< Number of rows */
  int ncols;   /**< Number of columns */
} Zmatrix;

/** \name Complex Matrix allocation functions 
 \{ */

/** Create a new matrix of nrows * ncols containing complex numbers
  @param nrows Number of rows in new matrix
  @param ncols Number of columns in new matrix
  @result Newly created Zmatrix of size nrows * ncols
*/
Zmatrix *zmat_new(int nrows, int ncols);

/** \name Complex Matrix cleanup functions 
 \{ */


/** Free a Zmatrix
   @param m Zmatrix to free  */
void zmat_free(Zmatrix *m);

/** \name Complex Matrix data access functions 
 \{ */

/* we'll only inline the functions likely to be used heavily in inner
   loops */  

/** Get a single complex number from matrix.
    @param m Matrix containing element to get
    @param row Row within Matrix containing element to get
    @param col Column within Matrix containing element to get
    @result Complex number specified by row & column within matrix
    */
static PHAST_INLINE
Complex zmat_get(Zmatrix *m, int row, int col) {
  return m->data[row][col];
}

/** Set a single complex number in matrix.
    @param m Matrix containing element to set
    @param row Row within matrix containing element to set
    @param col Column within matrix containing element to set
    @param val Value to set element of matrix to (identified by row & column)
 */
static PHAST_INLINE
void zmat_set(Zmatrix *m, int row, int col, Complex val) {
  m->data[row][col] = val;
}

/** Get a single row of complex numbers from a matrix.
  @param m Matrix to get data from
  @param row The row in matrix m to get data from 
  @result Vector containing data for specified row of matrix
*/
Zvector *zmat_get_row(Zmatrix *m, int row);

/** Get a single column of complex numbers from a matrix.
  @param m Matrix to get data from
  @param col Column in matrix m to get data from
  @result Vector containing data from specified column of matrix
*/
Zvector *zmat_get_col(Zmatrix *m, int col);

/** \name Complex Matrix modification functions 
 \{ */

/** Transform specified matrix into identity matrix.
    @param m Matrix to transform into identity matrix
 */
void zmat_set_identity(Zmatrix *m);

/** Set all elements in a matrix to zero.
    @param m Matrix to zero out
 */
void zmat_zero(Zmatrix *m);

/** Set all elements in a matrix to a single value
    @param m Matrix to modify
    @param val New value for every element in matrix
 */
void zmat_set_all(Zmatrix *m, Complex val);


/** Multiply all elements in matrix m by a complex scale_factor.
    @param m Matrix to multiply
    @param scale_factor Amount each element in matrix m should be multiplied by
    @note For real number scale version see zmat_scale()
    @see zmat_scale
 */
void zmat_scale_complex(Zmatrix *m, Complex scale_factor);

/** Multiply all elements in matrix m by a scale_factor of type double. 
    @param m Matrix to multiply
    @param scale_factor Amount each element in matrix m should be multiplied by
    @note For complex number scale version see zmat_scale()
    @see zmat_scale_complex
*/
void zmat_scale(Zmatrix *m, double scale_factor);

/** Matrix multiplication.
  @param[out] prod Result of matrix multiplication
  @param[in] m1 Matrix to multiply by m2
  @param[in] m2 Matrix to be multiplied by m1
  @warning Matrix prod will be overwritten with result
*/
void zmat_mult(Zmatrix *prod, Zmatrix *m1, Zmatrix *m2);

/** Matrix vector multiplication.
  @param[out] prod Result of matrix vector multiplication
  @param[in] m Matrix to be multiplied by v
  @param[in] v Vector to multiply matrix m by
  @warning Matrix prod will be overwritten with result
*/
void zmat_vec_mult(Zvector *prod, Zmatrix *m, Zvector *v);

/** Matrix Multiply complex matrices where result are expected to be real
  @param[out] prod result of matrix multiplication
  @param m1 complex matrix to multiply by m2
  @param m2 complex matrix to be multiplied by m1
  @warning the matrix prod will be overwritten
  @note the results matrix is expected to have real number results
*/
void zmat_mult_real(Matrix *prod, Zmatrix *m1, Zmatrix *m2);

/** Matrix Add.
  For each element [x,y] in thism add element [x,y] from addm. 
  param thism[in,out] Matrix that contains the element to add to
  param addm[in] Matrix that specifies how much to add to each element of thism
*/
void zmat_plus_eq(Zmatrix *thism, Zmatrix *addm);

/** Matrix Subtract.
  For each element[x,y] in thism, subtract element from subm at [x,y].
  param thism[in,out] Matrix that contains the element to subtract from
  param subm[in] Matrix that specifies how much to subtract from each element of thism
*/
void zmat_minus_eq(Zmatrix *thism, Zmatrix *subm);
/** Compute A = B * C * D where A, B, C, D are square matrices of the
   same dimension, and C is diagonal.  Allow B, C, D to be complex
   valued but assume their product is real valued (as when B,C,D
   represent diagonalization of A).  C is described by a vector
   representing its diagonal elements.  A temporary matrix can optionally
   be passed in to improve efficiency. 
   @param A Result matrix
   @param B First matrix in multiplication
   @param C Diagonal matrix represented by a vector
   @param D Last matrix in multiplication
   @param scratch (Optional) Improve efficiency by allocating scratch memory beforehand
 */
void zmat_mult_real_diag(Matrix *A, Zmatrix *B, Zvector *C, Zmatrix *D,
                         Zmatrix *scratch);

/** \name Complex Matrix copy functions 
 \{ */

/** Copy all data from one matrix 'src' to another 'dest'.
    @param src Matrix with data to copy from
    @param dest Matrix to copy data to
 */
void zmat_copy(Zmatrix *dest, Zmatrix *src);

/** Create a copy of all data from matrix 'src'.
    @param src Matrix to clone
    @result Clone of matrix specified in src
 */
Zmatrix *zmat_create_copy(Zmatrix *src);

/** \name Complex Matrix read/save as file functions 
 \{ */

/** Save a Matrix to file.
    @param m Matrix to save to file
    @param F File to save matrix in
 */
void zmat_print(Zmatrix *m, FILE *F);

/** Read a Matrix from a file into an existing matrix object. 
  @pre Matrix m must have defined the number of rows and columns
  @param Allocated Matrix object
  @param F file descriptor of file containing matrix data
*/
void zmat_read(Zmatrix *m, FILE *F);

/** Create a Matrix from matrix data in a file 
  @param F File descriptor of file containing matrix data
  @param nrows Number of rows allocated for new matrix
  @param ncols Number of columns allocated for new matrix
  @result Newly created z-matrix populated by data from file
*/
Zmatrix *zmat_new_from_file(FILE *F, int nrows, int ncols);

/** \name Complex Matrix cast functions 
 \{ */

/** "Cast" complex matrix as real, (extract real component of each element).
   @param[out] dest Result Matrix of real numbers
   @param[in] src Complex number matrix containing data we will "cast"
   @param[in] strict Ensure imaginary components are very close to zero, otherwise die
   @result 0 == Success, 1 == Failure (at least one element had an imaginary component that was not very close to zero)
*/
int zmat_as_real(Matrix *dest, Zmatrix *src, int strict);

/** \} */

#endif
