/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: markov_matrix.h,v 1.3 2008-11-12 02:07:59 acs Exp $ */

#ifndef MARKOVMAT_H
#define MARKOVMAT_H

#include <matrix.h>
#include <vector.h>
#include <complex_matrix.h>
#include <complex_vector.h>

#define NCHARS 256

typedef enum {DISCRETE, CONTINUOUS} mm_type;

typedef struct {
  Matrix *matrix;
  Zmatrix *evec_matrix;
  Zmatrix *evec_matrix_inv;
  Zvector *evals;
  int size;
  char *states;
  int inv_states[NCHARS];
  mm_type type;
} MarkovMatrix;

MarkovMatrix* mm_new(int size, char *states, mm_type type);
MarkovMatrix* mm_new_from_matrix(Matrix *A, char *states, mm_type type); 
MarkovMatrix* mm_new_from_file(FILE *F, mm_type type); 
void mm_free(MarkovMatrix *M); 
int mm_validate(MarkovMatrix *M); 
double mm_get(MarkovMatrix *M, int row, int col); 
void mm_set(MarkovMatrix *M, int row, int col, double val);
double mm_get_by_state(MarkovMatrix *M, char from, char to);
void mm_pretty_print(FILE *F, MarkovMatrix *M); 
void mm_exp(MarkovMatrix *P, MarkovMatrix *Q, double t);
void mm_diag_exp(Matrix *dest, Matrix *src, 
		 double t);
int mm_sample_state(MarkovMatrix *M, int state);
char mm_sample_char(MarkovMatrix *M, char c);
int mm_sample_vector(Vector *v);
char mm_sample_backgd(char *labels, Vector *backgd);
void mm_cpy(MarkovMatrix *dest, MarkovMatrix *src);
MarkovMatrix *mm_create_copy(MarkovMatrix *src);
MarkovMatrix* mm_new_from_counts(Matrix *counts, char *states);
void mm_diagonalize(MarkovMatrix *M); 
void mm_scale(MarkovMatrix *M, double scale);
void mm_renormalize(MarkovMatrix *M);

/* allows shorthand for element access */
extern inline
double mm_get(MarkovMatrix *M, int row, int col) {
  return(mat_get(M->matrix, row, col));
}

extern inline
void mm_set(MarkovMatrix *M, int row, int col, double val) {
  mat_set(M->matrix, row, col, val);
}

#endif
