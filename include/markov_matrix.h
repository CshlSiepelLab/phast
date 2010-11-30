/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: markov_matrix.h,v 1.4 2008-11-16 02:32:54 acs Exp $ */

#ifndef MARKOVMAT_H
#define MARKOVMAT_H

#include <matrix.h>
#include <vector.h>
#include <complex_matrix.h>
#include <complex_vector.h>
#include <external_libs.h>

#define NCHARS 256

typedef enum {DISCRETE, CONTINUOUS} mm_type;
typedef enum {REAL_NUM, COMPLEX_NUM} number_type;

typedef struct {
  Matrix *matrix;

  number_type eigentype;

  /* these are used if eigen_type == COMPLEX_NUM (general case) */
  Zmatrix *evec_matrix_z;
  Zmatrix *evec_matrix_inv_z;
  Zvector *evals_z;

  /* these are used if eigen_type == REAL_NUM (always true if reversible,
     results in significant savings in computation) */
  Matrix *evec_matrix_r;
  Matrix *evec_matrix_inv_r;
  Vector *evals_r;

  int diagonalize_error;  /** This is -1 if diagonalization has not been
                              attempted, 0 if diagonalization has succeeded,
                              and 1 if diagonalization has failed */
                   

  int size;
  char *states;
  int inv_states[NCHARS];
  mm_type type;
} MarkovMatrix;

MarkovMatrix* mm_new(int size, const char *states, mm_type type);
MarkovMatrix* mm_new_from_matrix(Matrix *A, const char *states, mm_type type); 
MarkovMatrix* mm_new_from_file(FILE *F, mm_type type); 
void mm_free(MarkovMatrix *M); 
void mm_free_eigen(MarkovMatrix *M);
void mm_set_eigentype(MarkovMatrix *M, number_type eigentype);
int mm_validate(MarkovMatrix *M); 
double mm_get_by_state(MarkovMatrix *M, char from, char to);
void mm_pretty_print(FILE *F, MarkovMatrix *M); 
void mm_exp(MarkovMatrix *P, MarkovMatrix *Q, double t);
int mm_sample_state(MarkovMatrix *M, int state);
char mm_sample_char(MarkovMatrix *M, char c);
int mm_sample_vector(Vector *v);
char mm_sample_backgd(char *labels, Vector *backgd);
void mm_cpy(MarkovMatrix *dest, MarkovMatrix *src);
MarkovMatrix *mm_create_copy(MarkovMatrix *src);
MarkovMatrix* mm_new_from_counts(Matrix *counts, const char *states);
void mm_diagonalize(MarkovMatrix *M); 
void mm_scale(MarkovMatrix *M, double scale);
void mm_renormalize(MarkovMatrix *M);

/* allows shorthand for element access */
static PHAST_INLINE
double mm_get(MarkovMatrix *M, int row, int col) {
  return(mat_get(M->matrix, row, col));
}

static PHAST_INLINE
void mm_set(MarkovMatrix *M, int row, int col, double val) {
  mat_set(M->matrix, row, col, val);
}

#endif
