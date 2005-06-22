/* $Id: markov_matrix.c,v 1.4 2005-06-22 07:11:19 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

/* functions for manipulating continuous and discrete Markov matrices */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <matrix.h>
#include <markov_matrix.h>
#include <complex.h>
#include <misc.h>
#include <eigen.h>

#define SUM_EPSILON 0.0001
#define ELEMENT_EPSILON 0.00001
#define MAXALPHA 1000

MarkovMatrix* mm_new(int size, char *states, mm_type type) {
  int i, alph_size;
  MarkovMatrix *M = (MarkovMatrix*)smalloc(sizeof(MarkovMatrix));
  M->evec_matrix = M->evec_matrix_inv = NULL;
  M->evals = NULL;
  M->matrix = mat_new(size, size);
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

MarkovMatrix* mm_new_from_matrix(Matrix *A, char *states, mm_type type) {
  MarkovMatrix *M = mm_new(A->nrows, states, type);
  M->matrix = A;
  mm_validate(M);
  return M;
}

/* discrete MM only */
MarkovMatrix* mm_new_from_counts(Matrix *counts, char *states) {
  int i, j;
  double rowsum;
  MarkovMatrix *M = mm_new(counts->nrows, states, DISCRETE);
  M->matrix = mat_new(counts->nrows, counts->ncols);
  for (i = 0; i < counts->nrows; i++) {
    for (j = 0, rowsum = 0; j < counts->ncols; j++) 
      rowsum += mat_get(counts, i, j);
    if (rowsum == 0) mat_set(M->matrix, i, i, 1); /* keeps valid */
    else 
      for (j = 0; j < counts->ncols; j++) 
        mat_set(M->matrix, i, j, 
                       safediv(mat_get(counts, i, j), rowsum));
                                /* safediv no longer necessary? */
  }
  mm_validate(M);
  return M;  
}

void mm_free(MarkovMatrix *M) {
  if (M->matrix != NULL)
    mat_free(M->matrix);
  if (M->evec_matrix != NULL)
    zmat_free(M->evec_matrix);
  if (M->evec_matrix_inv != NULL)
    zmat_free(M->evec_matrix_inv);
  if (M->evals != NULL)
    zvec_free(M->evals);
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
  if (M->matrix->nrows != M->matrix->ncols)
    die("ERROR: Matrix is not square (%d x %d)", 
	(int)M->matrix->nrows, (int)M->matrix->ncols);    
  M->size = M->matrix->nrows;

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
  return(mat_get(M->matrix, row, col));
}

void mm_set(MarkovMatrix *M, int row, int col, double val) {
  mat_set(M->matrix, row, col, val);
}

/* element access by state character 
   FIXME: this won't work with higher-order models */
double mm_get_by_state(MarkovMatrix *M, char from, char to) {
  return(mat_get(M->matrix, M->inv_states[(int)from], 
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
  M->matrix = mat_new_from_file(F, sz, sz);
  mm_validate(M);
  return M;
}

void mm_pretty_print(FILE *F, MarkovMatrix *M) {
  int j;

  for (j = 0; j < M->size; j++) fprintf(F, "%c\t", M->states[j]);
  fprintf(F, "\n");

  mat_print(M->matrix, F);
}

/* computes discrete matrix P by the formula P = exp(Qt),
   given Q and t */
void mm_exp(MarkovMatrix *P, MarkovMatrix *Q, double t) {

  static Zmatrix *Eexp = NULL; /* reuse these if possible */
  static Zmatrix *tmp = NULL;
  static int last_size = 0;
  int n = Q->size;
  int i, j;

  assert(P->size == Q->size && t >= 0);

  if (t == 0) {
    mat_set_identity(P->matrix);
    return;
  }

  if (last_size != Q->size && Eexp != NULL) {
    zmat_free(Eexp);
    zmat_free(tmp);
    Eexp = NULL;
  }

  if (Eexp == NULL) {
    Eexp = zmat_new(Q->size, Q->size);
    tmp = zmat_new(Q->size, Q->size);
    last_size = Q->size;
  }

  /* Diagonalize (if necessary) */
  if (Q->evec_matrix == NULL || Q->evals == NULL || 
      Q->evec_matrix == NULL) 
    mm_diagonalize(Q);

  /* Now compute P(t) = S exp(Dt) S^-1.  Start by computing exp(Dt) S^-1 */
  for (i = 0; i < n; i++) {
    Complex exp_dt_i =
      z_exp(z_mul_real(zvec_get(Q->evals, i), t));
    for (j = 0; j < n; j++) 
      zmat_set(tmp, i, j, z_mul(exp_dt_i, zmat_get(Q->evec_matrix_inv, i, j)));
  }

  /* Now multiply by S (on the left) */
  zmat_mult_real(P->matrix, Q->evec_matrix, tmp);
}

/* given a diagonal matrix with elements l_i,i create a new matrix
   with elements exp(l_i,i * t) */
/* void mm_diag_exp(Matrix *dest, Matrix *src,  */
/*         double t) { */
/*   int i; */
/*   assert((src->nrows == src->ncols) && (dest->nrows == dest->ncols) */
/*      && (src->nrows == dest->nrows)); */
/*   mat_zero(dest); */
/*   for (i = 0; i < src->nrows; i++) { */
/*     double expon = mat_get(src, i, i) * t; */
/*     mat_set(dest, i, i, exp(expon)); */
/*   } */
/* } */

/* given a state, draw the next state from the multinomial
 * distribution defined by the corresponding row in the matrix */
int mm_sample_state(MarkovMatrix *M, int state) {

  /* this is a bit inefficient, but we'll keep it for now, because it
   * simplifies the code */
  Vector *v = mat_get_row(M->matrix, state);
  return (mm_sample_vector(v));
}

/* given a probability vector, make a draw, and return the
 * corresponding index. */ 
int mm_sample_vector(Vector *v) {
  double r;
  int i;
  double sum;

  r = rand()*1.0/RAND_MAX;

  sum = 0;
  for (i = 0; i < v->size; i++) {
    sum += vec_get(v, i);
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
char mm_sample_backgd(char *labels, Vector *backgd) {
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
  mat_copy(dest->matrix, src->matrix);
  if (src->evec_matrix != NULL)
    zmat_copy(dest->evec_matrix, src->evec_matrix);
  if (src->evals != NULL)
    zvec_copy(dest->evals, src->evals);
  if (src->evec_matrix_inv != NULL)
    zmat_copy(dest->evec_matrix_inv, src->evec_matrix_inv);
}

MarkovMatrix *mm_create_copy(MarkovMatrix *src) {
  MarkovMatrix *retval = mm_new(src->size, src->states, src->type);
  if (src->evec_matrix != NULL)
    retval->evec_matrix = zmat_new(src->size, src->size);
  if (src->evals != NULL)
    retval->evals = zvec_new(src->size);
  if (src->evec_matrix_inv != NULL)
    retval->evec_matrix_inv = zmat_new(src->size, src->size);
  mm_cpy(retval, src);
  return(retval);
}

void mm_diagonalize(MarkovMatrix *M) { 
  if (M->evec_matrix == NULL)
    M->evec_matrix = zmat_new(M->size, M->size);
  if (M->evals == NULL)
    M->evals = zvec_new(M->size);
  if (M->evec_matrix_inv == NULL)
    M->evec_matrix_inv = zmat_new(M->size, M->size);
  mat_diagonalize(M->matrix, M->evals, M->evec_matrix, M->evec_matrix_inv);
} 

void mm_scale(MarkovMatrix *M, double scale) {
  mat_scale(M->matrix, scale);
}

/* Renormalize a discrete Markov matrix so that all rows sum to 1. */
void mm_renormalize(MarkovMatrix *M) {
  int i, j;
  assert(M->type == DISCRETE);
  for (i = 0; i < M->size; i++) {
    double rowsum = 0;
    for (j = 0; j < M->size; j++) 
      rowsum += mat_get(M->matrix, i, j);

    if (rowsum == 0)            /* in this case, assume an "absorbing"
                                   state, with prob 1 of a self-transition */
      mat_set(M->matrix, i, i, 1);

    else if (rowsum != 1)       /* renormalize */
      for (j = 0; j < M->size; j++) 
        mat_set(M->matrix, i, j,
                       mat_get(M->matrix, i, j) / rowsum);
                     
  }
}
