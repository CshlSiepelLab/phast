/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/**
  @file markov_matrix.h
  Manipulating continuous and discrete Markov matrices. 
  @ingroup base
*/
#ifndef MARKOVMAT_H
#define MARKOVMAT_H

#include <matrix.h>
#include <vector.h>
#include <complex_matrix.h>
#include <complex_vector.h>
#include <external_libs.h>

/** Size of invariant states char array. */
#define NCHARS 256

/** Type of Markov Matrix. */
typedef enum {DISCRETE, /**< Discrete Markov Matrix */
	     CONTINUOUS /**< Continuous Markov Matrix */
	     } mm_type;
/** Number Types used in Markov Matrix. */
typedef enum {REAL_NUM, /**< Real numbers used in matrix */
	     COMPLEX_NUM /**< Complex numbers used in matrix */
	     } number_type;

/** Markov Matrix object */
typedef struct {
  Matrix *matrix; /**< Matrix holding real or complex values */

  number_type eigentype; /**< Whether matrix values are complex or real */ 

  /* these are used if eigen_type == COMPLEX_NUM (general case) */
  Zmatrix *evec_matrix_z; /**< Eigen vector matrix (if eigen_type == COMPLEX_NUM) */
  Zmatrix *evec_matrix_inv_z; /**< Eigen vector inverted matrix (if eigen_type == COMPLEX_NUM) */
  Zvector *evals_z;  /**< Eigen Values (if eigen_type == COMPLEX_NUM) */

  /* these are used if eigen_type == REAL_NUM (always true if reversible,
     results in significant savings in computation) */
  Matrix *evec_matrix_r; /**< Eigen vector matrix (if eigen_type == REAL_NUM) */
  Matrix *evec_matrix_inv_r; /**< Eigen vector inverted matrix (if eigen_type == REAL_NUM) */
  Vector *evals_r;  /**< Eigen values (if eigen_type == REAL_NUM) */

  int diagonalize_error;  /**< Status of diagonalization: -1=diagonalization has not been
                              attempted, 0=diagonalization has succeeded,
                              1=diagonalization has failed */
  int size; /**< Size of matrix */
  char *states; /**< Lookup of state character from state number */
  int inv_states[NCHARS]; /**< Inverse table, for lookup of state number from state character  */
  mm_type type; /**< Whether matrix is Discrete or Continuous */
} MarkovMatrix;

/** \name Markov Matrix allocation functions.
  \{ */

/** Create a new empty Markov Matrix.
    @param size Number of rows in Markov Matrix (#rows == #col)
    @param[in] states  String containing possible states, one for each character
    @param[in] type Whether the matrix is Continuous or Discrete
    @result Newly allocated Markov Matrix with parameters as specified.
*/
MarkovMatrix* mm_new(int size, const char *states, mm_type type);

/** Create a new Markov Matrix from an already populated matrix.
    @param[in] A Matrix containing data to put into Markov Matrix
    @param[in] states String containing possible states, one for each character
    @param[in] type Whether the matrix is Continuous or Discrete
    @result Newly allocated Markov Matrix with data from matrix A
*/
MarkovMatrix* mm_new_from_matrix(Matrix *A, const char *states, mm_type type); 

/** Copy a Markov Matrix creating a new Markov Matrix object
    @param src Markov Matrix to copy
*/
MarkovMatrix *mm_create_copy(MarkovMatrix *src);

/** Create a new Discrete Markov Matrix using a Matrix of counts and a list of states.
  @param counts Matrix of counts corresponding to states
  @param states List of states in the form of a string (one per char)
  @result Newly created Markov Matrix with data
  @note Discrete Markov Matrices only
*/
MarkovMatrix* mm_new_from_counts(Matrix *counts, const char *states);

/** \} \name Markov Matrix cleanup functions. 
  \{ */

/** Free a Markov matrix. 
    @param[in,out] M Markov Matrix to free
*/
void mm_free(MarkovMatrix *M);

/** Free eigenvector and eigenvalue matrices/vectors, if allocated 
   @param[in,out] M Markov Matrix containing eigenvectors and eigenvalues to free
*/
void mm_free_eigen(MarkovMatrix *M);

/** \} \name Markov Matrix sample functions.
   \{ */

/** Given a state, draw the next state.
   @param M Markov Matrix containing next state
   @param state State used to determine next state
   @result Next state ID
*/
int mm_sample_state(MarkovMatrix *M, int state);

/** Given a character label, draw the next state (as a character label).
    @param M Markov Matrix containing next state
    @param c Character Label used to determine next state
    @result Character Label associated with the next state
*/
char mm_sample_char(MarkovMatrix *M, char c);
//int mm_sample_vector(Vector *v);
/**  Sample by background
   @param labels a character string indicating the Markov Matrix's alphabet, ie, "ACGT"
   @param backgd A vector giving the frequency of each state (values should be in (0,1) and sum to 1).
   @result A character representing the sampled state.
*/
char mm_sample_backgd(char *labels, Vector *backgd);

 /** \} \name Markov Matrix read/save file functions.
   \{  */

/** Create a new Markov Matrix from a file
    @param[in] F File containing Markov Matrix data
    @param[in] type Whether the matrix is continuous or Discrete
    @result Newly allocated Markov Matrix with data from the file
*/
MarkovMatrix* mm_new_from_file(FILE *F, mm_type type); 

/** Save the Markov Matrix to a file
    @param[in] F File to save to
    @param[in] M Markov Matrix to save
*/
void mm_pretty_print(FILE *F, MarkovMatrix *M); 

 /** \} \name Markov Matrix misc. functions
   \{  */


/** Define matrix as having real or complex eigenvectors/eigenvalues.
    @param[in,out] M Markov Matrix to set number type on. 
    @param[in] eigentype Type of numbers eigen vectors / eigen values should use
    @note If the matrix is to be diagonalized, exponentiated, etc.,
     significant saving in computation can be had by using real numbers
     rather than complex numbers where possible (e.g., with reversible
     Markov models). 
*/
void mm_set_eigentype(MarkovMatrix *M, number_type eigentype);

/** Validate that a Markov Matrix adheres to basic Markov Matrix requirements.
    Checks the following:
    - Matrix within Markov matrix is not NULL
    - Matrix is square
    - Rows of matrix sum to one (or zero)

@param[in] M Markov Matrix to check
@result 1 on failure; 0 on success
*/
int mm_validate(MarkovMatrix *M); 

/** Get likelihood for going from one given state to another given state.
    @param[in] M Markov Matrix to interrogate
    @param[in] from Beginning state
    @param[in] to Ending state
    @result Likelihood of going from -> to 
 */
double mm_get_by_state(MarkovMatrix *M, char from, char to);

/** Computes discrete matrix P by the formula P = exp(Qt), given Q and t. 
    @param[out] P Result Markov Matrix
    @param[in] Q Input Markov matrix
    @param[in] t Amount to scale Q by
*/
void mm_exp(MarkovMatrix *P, MarkovMatrix *Q, double t);

/** Copy a Markov Matrix into another existing Markov Matrix
    @param dest Where to copy the Markov Matrix to
    @param src Where to copy the Markov Matrix from
 */
void mm_cpy(MarkovMatrix *dest, MarkovMatrix *src);


/** Diagonalize a Markov Matrix.
    @param M Matrix to diagonalize
*/
void mm_diagonalize(MarkovMatrix *M);

/** Scale a Markov Matrix.
    @param M Matrix to scale
    @param scale Amount to scale matrix M by
*/ 
void mm_scale(MarkovMatrix *M, double scale);

/** Renormalize a discrete Markov Matrix so that all rows sum to 1. 
    @param M Matrix to renormalize
    @note validate would detect an un-normalized matrix
*/
void mm_renormalize(MarkovMatrix *M);

/* allows shorthand for element access */
/** Get an individual element of the Markov Matrix
    @param M Matrix to retrieve element from
    @param row Row the element is on within the matrix
    @param col Column the element is on within the matrix
    @result Element of the Matrix
 */
static PHAST_INLINE
double mm_get(MarkovMatrix *M, int row, int col) {
  return(mat_get(M->matrix, row, col));
}

/** Set an individual element of the Markov Matrix
    @param M Matrix to set element in
    @param row Row to set the element in
    @param col Column to set the element in
    @param val Value to set the element to
*/
static PHAST_INLINE
void mm_set(MarkovMatrix *M, int row, int col, double val) {
  mat_set(M->matrix, row, col, val);
}
/** \} */

#endif
