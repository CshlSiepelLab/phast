/* $Id: markov_matrix.h,v 1.1.1.1 2004-06-03 22:43:11 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

#ifndef MARKOVMAT_H
#define MARKOVMAT_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>

#define NCHARS 256

typedef enum {DISCRETE, CONTINUOUS} mm_type;

typedef struct {
  gsl_matrix *matrix;
  gsl_matrix_complex *evec_matrix;
  gsl_matrix_complex *evec_matrix_inv;
  gsl_vector_complex *evals;
  int size;
  char *states;
  int inv_states[NCHARS];
  mm_type type;
} MarkovMatrix;

MarkovMatrix* mm_new(int size, char *states, mm_type type);
MarkovMatrix* mm_new_from_matrix(gsl_matrix *A, char *states, mm_type type); 
MarkovMatrix* mm_new_from_file(FILE *F, mm_type type); 
void mm_free(MarkovMatrix *M); 
int mm_validate(MarkovMatrix *M); 
double mm_get(MarkovMatrix *M, int row, int col); 
void mm_set(MarkovMatrix *M, int row, int col, double val);
double mm_get_by_state(MarkovMatrix *M, char from, char to);
void mm_pretty_print(FILE *F, MarkovMatrix *M); 
void mm_exp(MarkovMatrix *P, MarkovMatrix *Q, double t);
void mm_diag_exp(gsl_matrix *dest, gsl_matrix *src, 
		 double t);
int mm_sample_state(MarkovMatrix *M, int state);
char mm_sample_char(MarkovMatrix *M, char c);
int mm_sample_vector(gsl_vector *v);
char mm_sample_backgd(char *labels, gsl_vector *backgd);
void mm_cpy(MarkovMatrix *dest, MarkovMatrix *src);
MarkovMatrix *mm_create_copy(MarkovMatrix *src);
MarkovMatrix* mm_new_from_counts(gsl_matrix *counts, char *states);
void mm_diagonalize(MarkovMatrix *M); 
void mm_scale(MarkovMatrix *M, double scale);
void mm_renormalize(MarkovMatrix *M);

/* allows shorthand for element access */
extern inline
double mm_get(MarkovMatrix *M, int row, int col) {
  return(gsl_matrix_get(M->matrix, row, col));
}

extern inline
void mm_set(MarkovMatrix *M, int row, int col, double val) {
  gsl_matrix_set(M->matrix, row, col, val);
}

#endif
