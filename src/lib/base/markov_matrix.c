/* $Id: markov_matrix.c,v 1.2 2004-06-23 21:37:14 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

/* functions for manipulating continuous and discrete Markov matrices */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"
#include "markov_matrix.h"
#include <gsl/gsl_complex_math.h>
#include "misc.h"

#define SUM_EPSILON 0.0001
#define ELEMENT_EPSILON 0.00001
#define MAXALPHA 1000

MarkovMatrix* mm_new(int size, char *states, mm_type type) {
  int i, alph_size;
  MarkovMatrix *M = (MarkovMatrix*)smalloc(sizeof(MarkovMatrix));
  M->evec_matrix = M->evec_matrix_inv = NULL;
  M->evals = NULL;
  M->matrix = gsl_matrix_calloc(size, size);
  M->size = size;
  alph_size = states == NULL ? size : strlen(states);
  M->states = (char*)smalloc((alph_size+1) * sizeof(char));

  if (states == NULL) {
    for (i = 0; i < alph_size; i++) 
      M->states[i] = i < 26 ? ('A' + (char)i) : 'Z'; 
                                /* avoid non-alphanumeric */
    /* FIXME: not actually using these for anything; probably should
       just rip them out */

    M->states[alph_size] = '\0';
  }
  else 
    strcpy(M->states, states);

  M->type = type;

  /* build inverse table, for lookup of state number from state
   * character */
  for (i = 0; i < NCHARS; i++) M->inv_states[i] = -1;
  for (i = 0; M->states[i] != '\0'; i++) M->inv_states[(int)M->states[i]] = i;

  return M;
}

MarkovMatrix* mm_new_from_matrix(gsl_matrix *A, char *states, mm_type type) {
  MarkovMatrix *M = mm_new(A->size1, states, type);
  M->matrix = A;
  mm_validate(M);
  return M;
}

/* discrete MM only */
MarkovMatrix* mm_new_from_counts(gsl_matrix *counts, char *states) {
  int i, j;
  double rowsum;
  MarkovMatrix *M = mm_new(counts->size1, states, DISCRETE);
  M->matrix = gsl_matrix_calloc(counts->size1, counts->size2);
  for (i = 0; i < counts->size1; i++) {
    for (j = 0, rowsum = 0; j < counts->size2; j++) 
      rowsum += gsl_matrix_get(counts, i, j);
    if (rowsum == 0) gsl_matrix_set(M->matrix, i, i, 1); /* keeps valid */
    else 
      for (j = 0; j < counts->size2; j++) 
        gsl_matrix_set(M->matrix, i, j, 
                       safediv(gsl_matrix_get(counts, i, j), rowsum));
                                /* safediv no longer necessary? */
  }
  mm_validate(M);
  return M;  
}

void mm_free(MarkovMatrix *M) {
  if (M->matrix != NULL)
    gsl_matrix_free(M->matrix);
  if (M->evec_matrix != NULL)
    gsl_matrix_complex_free(M->evec_matrix);
  if (M->evec_matrix_inv != NULL)
    gsl_matrix_complex_free(M->evec_matrix_inv);
  if (M->evals != NULL)
    gsl_vector_complex_free(M->evals);
  if (M->states != NULL)
    free(M->states);
  free(M);
}

/* returns 1 on failure, 0 on success */
int mm_validate(MarkovMatrix *M) {
  int i, j;
  double targetval;

  if (M->matrix == NULL) {
    fprintf(stderr, "Error validating Markov matrix: matrix undefined.\n");
    return 1;
  }

  /* ensure square */
  if (M->matrix->size1 != M->matrix->size2) {
    fprintf(stderr, "Warning validating Markov matrix: matrix is not square (%d x %d).  Will proceed using smaller dimension.\n", M->matrix->size1, M->matrix->size2);    
  }
  M->size = M->matrix->size1 <= M->matrix->size2 ? M->matrix->size1 : 
    M->matrix->size2;

  /* ensure rows sum to one (or zero) */
  targetval = M->type == DISCRETE ? 1 : 0;
  for (i = 0; i < M->size; i++) {
    double sum = 0;
    for (j = 0; j < M->size; j++) 
      sum += mm_get(M, i, j);
    if (abs(sum - targetval) > SUM_EPSILON) {
      fprintf(stderr, "Error validating Markov matrix: rows do not sum to %.1f (+-%f).\n", targetval, SUM_EPSILON);
      return 1;
    }
  }
  
  return 0;
}

/* allows shorthand for element access */
double mm_get(MarkovMatrix *M, int row, int col) {
  return(gsl_matrix_get(M->matrix, row, col));
}

void mm_set(MarkovMatrix *M, int row, int col, double val) {
  gsl_matrix_set(M->matrix, row, col, val);
}

/* element access by state character 
   FIXME: this won't work with higher-order models */
double mm_get_by_state(MarkovMatrix *M, char from, char to) {
  return(gsl_matrix_get(M->matrix, M->inv_states[(int)from], 
                        M->inv_states[(int)to]));
}

MarkovMatrix *mm_new_from_file(FILE *F, mm_type type) {
  /* expects list of 1-character labels for states to appear as first
     line; defines alphabet and order */
  int i, sz = 0;
  char line[2*MAXALPHA], *states;
  MarkovMatrix *M;

  fgets(line, 2*MAXALPHA, F);
  states = smalloc(sizeof(char) * strlen(line));
  for (i = 0; i < strlen(line); i++) {
    if (isalnum(line[i])) {
      states[sz++] = line[i];
    }
  }
  states[sz] = '\0';

  M = mm_new(sz, states, type);
  free(states);

  /* read matrix, according to dimension of labels */
  M->matrix = gsl_matrix_alloc(sz, sz);
  gsl_matrix_fscanf(F, M->matrix);

  mm_validate(M);
  return M;
}

void mm_pretty_print(FILE *F, MarkovMatrix *M) {
  int j;

  for (j = 0; j < M->size; j++) fprintf(F, "%c\t", M->states[j]);
  fprintf(F, "\n");

  gsl_matrix_pretty_print(F, M->matrix);
}

/* computes discrete matrix P by the formula P = exp(Qt),
   given Q and t */
void mm_exp(MarkovMatrix *P, MarkovMatrix *Q, double t) {

  static gsl_matrix_complex *Eexp = NULL; /* reuse these if possible */
  static gsl_matrix_complex *tmp = NULL;
  static int last_size = 0;
  int n = Q->size;
  int i, j;

  assert(P->size == Q->size && t >= 0);

  if (t == 0) {
    gsl_matrix_set_identity(P->matrix);
    return;
  }

  if (last_size != Q->size && Eexp != NULL) {
    gsl_matrix_complex_free(Eexp);
    gsl_matrix_complex_free(tmp);
    Eexp = NULL;
  }

  if (Eexp == NULL) {
    Eexp = gsl_matrix_complex_alloc(Q->size, Q->size);
    tmp = gsl_matrix_complex_alloc(Q->size, Q->size);
    last_size = Q->size;
  }

  /* Diagonalize (if necessary) */
  if (Q->evec_matrix == NULL || Q->evals == NULL || 
      Q->evec_matrix == NULL) 
    mm_diagonalize(Q);

  /* Now compute P(t) = S exp(Dt) S^-1.  Start by computing exp(Dt) S^-1 */
  for (i = 0; i < n; i++) {
    gsl_complex exp_dt_i =
      gsl_complex_exp(gsl_complex_mul_real(gsl_vector_complex_get(Q->evals, i), t));
    for (j = 0; j < n; j++) 
      gsl_matrix_complex_set(tmp, i, j, gsl_complex_mul(exp_dt_i, gsl_matrix_complex_get(Q->evec_matrix_inv, i, j)));
  }

  /* Now multiply by S (on the left) */
  mat_mult_complex_to_real(P->matrix, Q->evec_matrix, tmp);
}

/* given a diagonal matrix with elements l_i,i create a new matrix
   with elements exp(l_i,i * t) */
/* void mm_diag_exp(gsl_matrix *dest, gsl_matrix *src,  */
/*         double t) { */
/*   int i; */
/*   assert((src->size1 == src->size2) && (dest->size1 == dest->size2) */
/*      && (src->size1 == dest->size1)); */
/*   gsl_matrix_set_zero(dest); */
/*   for (i = 0; i < src->size1; i++) { */
/*     double expon = gsl_matrix_get(src, i, i) * t; */
/*     gsl_matrix_set(dest, i, i, exp(expon)); */
/*   } */
/* } */

/* given a state, draw the next state from the multinomial
 * distribution defined by the corresponding row in the matrix */
int mm_sample_state(MarkovMatrix *M, int state) {

  /* this is a bit inefficient, but we'll keep it for now, because it
   * simplifies the code */
  gsl_vector_view v = gsl_matrix_row(M->matrix, state);
  return (mm_sample_vector(&v.vector));
}

/* given a multinomial vector, make a draw, and return the
 * corresponding index. Breaking this out from mm_sample_state allows
 * easy sampling from backgd distribution */ 
int mm_sample_vector(gsl_vector *v) {
  double r;
  int i;
  double sum;

  r = rand()*1.0/RAND_MAX;

  sum = 0;
  for (i = 0; i < v->size; i++) {
    sum += gsl_vector_get(v, i);
    if (sum > r) break;
  }

  if (i == v->size) {
    /* this seems to happen occasionally: check for bad prob vector or
       bad draw */
    if (abs(1 - sum) > SUM_EPSILON) 
      die("ERROR: unnormalized probability vector in mm_sample_vector (sum = %f)\n", sum);
    else if (r > 1)
      die("ERROR: bad random draw in mm_sample_vector (r = %f)\n", r);
    else i = v->size-1;         /* prob. just slight rounding error */
  }

  return i;
}

/* similar to above, but by character; for convenience in sampling
 * from background distrib */
char mm_sample_backgd(char *labels, gsl_vector *backgd) {
  return labels[mm_sample_vector(backgd)];
}


/* given a state, draw the next state from the multinomial
 * distribution defined by the corresponding row in the matrix, but
 * identify both the "from" and "to" states by character label */
char mm_sample_char(MarkovMatrix *M, char c) {
  return(M->states[mm_sample_state(M, M->inv_states[(int)c])]);
}

/* WARNING: assumes matrices in dest already allocated and of correct
 * size.  Also assumes type, states, and size are the same */
void mm_cpy(MarkovMatrix *dest, MarkovMatrix *src) {
  gsl_matrix_memcpy(dest->matrix, src->matrix);
  if (src->evec_matrix != NULL)
    gsl_matrix_complex_memcpy(dest->evec_matrix, src->evec_matrix);
  if (src->evals != NULL)
    gsl_vector_complex_memcpy(dest->evals, src->evals);
  if (src->evec_matrix_inv != NULL)
    gsl_matrix_complex_memcpy(dest->evec_matrix_inv, src->evec_matrix_inv);
}

MarkovMatrix *mm_create_copy(MarkovMatrix *src) {
  MarkovMatrix *retval = mm_new(src->size, src->states, src->type);
  if (src->evec_matrix != NULL)
    retval->evec_matrix = gsl_matrix_complex_alloc(src->size, src->size);
  if (src->evals != NULL)
    retval->evals = gsl_vector_complex_alloc(src->size);
  if (src->evec_matrix_inv != NULL)
    retval->evec_matrix_inv = gsl_matrix_complex_alloc(src->size, src->size);
  mm_cpy(retval, src);
  return(retval);
}

void mm_diagonalize(MarkovMatrix *M) { 
  if (M->evec_matrix == NULL)
    M->evec_matrix = gsl_matrix_complex_alloc(M->size, M->size);
  if (M->evals == NULL)
    M->evals = gsl_vector_complex_alloc(M->size);
  if (M->evec_matrix_inv == NULL)
    M->evec_matrix_inv = gsl_matrix_complex_alloc(M->size, M->size);
  mat_diagonalize(M->matrix, M->evals, M->evec_matrix, M->evec_matrix_inv);
} 

void mm_scale(MarkovMatrix *M, double scale) {
  gsl_matrix_scale(M->matrix, scale);
}

/* Renormalize a discrete Markov matrix so that all rows sum to 1. */
void mm_renormalize(MarkovMatrix *M) {
  int i, j;
  assert(M->type == DISCRETE);
  for (i = 0; i < M->size; i++) {
    double rowsum = 0;
    for (j = 0; j < M->size; j++) 
      rowsum += gsl_matrix_get(M->matrix, i, j);

    if (rowsum == 0)            /* in this case, assume an "absorbing"
                                   state, with prob 1 of a self-transition */
      gsl_matrix_set(M->matrix, i, i, 1);

    else if (rowsum != 1)       /* renormalize */
      for (j = 0; j < M->size; j++) 
        gsl_matrix_set(M->matrix, i, j,
                       gsl_matrix_get(M->matrix, i, j) / rowsum);
                     
  }
}
