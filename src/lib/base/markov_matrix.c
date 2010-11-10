/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: markov_matrix.c,v 1.11 2009-02-19 21:29:24 agd27 Exp $ */

/* functions for manipulating continuous and discrete Markov matrices */

#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <matrix.h>
#include <markov_matrix.h>
#include <complex.h>
#include <misc.h>
#include <eigen.h>
#include <prob_vector.h>

#define SUM_EPSILON 0.0001
#define ELEMENT_EPSILON 0.00001
#define MAXALPHA 1000

MarkovMatrix* mm_new(int size, const char *states, mm_type type) {
  int i, alph_size;
  MarkovMatrix *M = (MarkovMatrix*)smalloc(sizeof(MarkovMatrix));
  M->evec_matrix_z = M->evec_matrix_inv_z = NULL;
  M->evec_matrix_r = M->evec_matrix_inv_r = NULL;
  M->evals_z = NULL;
  M->evals_r = NULL;
  M->matrix = mat_new(size, size);
  mat_zero(M->matrix);
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
  M->eigentype = COMPLEX_NUM;       /* default */

  /* build inverse table, for lookup of state number from state
   * character */
  for (i = 0; i < NCHARS; i++) M->inv_states[i] = -1;
  for (i = 0; M->states[i] != '\0'; i++) M->inv_states[(int)M->states[i]] = i;

  return M;
}

MarkovMatrix* mm_new_from_matrix(Matrix *A, const char *states, mm_type type) {
  MarkovMatrix *M = mm_new(A->nrows, states, type);
  mat_free(M->matrix);
  M->matrix = A;
  mm_validate(M);
  return M;
}

/* discrete MM only */
MarkovMatrix* mm_new_from_counts(Matrix *counts, const char *states) {
  int i, j;
  double rowsum;
  MarkovMatrix *M = mm_new(counts->nrows, states, DISCRETE);
  M->matrix = mat_new(counts->nrows, counts->ncols);
  mat_zero(M->matrix);
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
  if (M->states != NULL)
    sfree(M->states);
  mm_free_eigen(M);
  sfree(M);
}

/* free eigenvector and eigenvalue matrices/vectors, if allocated */
void mm_free_eigen(MarkovMatrix *M) {
  if (M->evec_matrix_r != NULL)
    mat_free(M->evec_matrix_r);
  if (M->evec_matrix_inv_r != NULL)
    mat_free(M->evec_matrix_inv_r);
  if (M->evals_r != NULL)
    vec_free(M->evals_r);
  if (M->evec_matrix_z != NULL)
    zmat_free(M->evec_matrix_z);
  if (M->evec_matrix_inv_z != NULL)
    zmat_free(M->evec_matrix_inv_z);
  if (M->evals_z != NULL)
    zvec_free(M->evals_z);
  M->evec_matrix_r = M->evec_matrix_inv_r = NULL;
  M->evals_r = NULL;
  M->evec_matrix_z = M->evec_matrix_inv_z = NULL;
  M->evals_z = NULL;
}

/* define matrix as having real or complex eigenvectors/eigenvalues.
   If the matrix is to be diagonalized, exponentiated, etc.,
   significant saving in computation can be had by using real numbers
   rather than complex numbers where possible (e.g., with reversible
   Markov models).  */
void mm_set_eigentype(MarkovMatrix *M, number_type eigentype) {
  M->eigentype = eigentype;
  mm_free_eigen(M);             /* instantiate new ones on demand */
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
    double sum = 0.0, diag;
    for (j = 0; j < M->size; j++) 
      if (i != j)
	sum += mm_get(M, i, j);
    diag = mm_get(M, i, i);
    if (abs(sum + diag - targetval) > SUM_EPSILON) {
      fprintf(stderr, "Error validating Markov  matrix: rows do not sum to %.1f (+-%f). %f %f\n", targetval, SUM_EPSILON, sum, diag);
      return 1;
    } /*else if (sum != -diag)
	mm_set(M, i, i, -sum);*/
  }
  return 0;
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

  if (fgets(line, 2*MAXALPHA, F)==NULL)
    die("ERROR reading markov matrix\n");
  states = smalloc(sizeof(char) * strlen(line));
  for (i = 0; i < strlen(line); i++) {
    if (isalnum(line[i])) {
      states[sz++] = line[i];
    }
  }
  states[sz] = '\0';

  M = mm_new(sz, states, type);
  sfree(states);

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


void mm_exp_taylor(MarkovMatrix *P, MarkovMatrix *Q, double t) {
  Matrix *Qt, *lastQ, *newQ, *lastP;
  Matrix *sumP;
  double diff, fac=1.0, d;
  int rep, i, j, n = P->size;

  sumP = mat_new(n, n);
  Qt = mat_create_copy(Q->matrix);
  mat_scale(Qt, t);
  
  mat_set_identity(sumP);
  
  mat_plus_eq(sumP, Qt);
  
  lastQ = mat_create_copy(Qt);
  newQ = mat_new(n, n);
  lastP = mat_new(n, n);
  
  for (rep=2; 1; rep++) {
    mat_mult(newQ, Qt, lastQ);  //lastQ = (Qt)^(rep-1), newQ=(Qt)^rep
    mat_copy(lastQ, newQ);    // now lastQ is (Qt)^rep too
    fac /= (double)rep;
    if (fac == 0.0) break;
    mat_scale(newQ, fac);      //scale newQ by 1/rep!
    mat_copy(lastP, sumP);     //put previous result in lastP
    mat_plus_eq(sumP, newQ);     //add new term
    //check convergence
    diff = 0;
    for (i=0; i < n; i++) {
      for (j=0; j < n; j++) {
	d = mat_get(sumP, i, j) - mat_get(lastP, i, j);
	diff += (d*d);
      }
    }
    if (diff == 0.0 || (diff < 1.0e-8 && rep >= 10)) break;
  }
  //  printf("done manual exponentiation rep=%i diff=%e\n", rep, diff);
  mat_copy(P->matrix, sumP);
  mat_free(sumP);
  mat_free(lastP);
  mat_free(lastQ);
  mat_free(newQ);
  mat_free(Qt);
  return;
}


/* general version allowing for complex eigenvalues/eigenvectors */
void mm_exp_complex(MarkovMatrix *P, MarkovMatrix *Q, double t) {

  static Zmatrix *Eexp = NULL; /* reuse these if possible */
  static Zmatrix *tmp = NULL;
  static int last_size = 0;
  int n = Q->size;
  int i, j;

  if (!(P->size == Q->size && t >= 0))
    die("ERROR mm_exp_complex: got P->size=%i, Q->sizse=%i, t=%f\n",
	P->size, Q->size, t);

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
    set_static_var((void**)&Eexp);
    tmp = zmat_new(Q->size, Q->size);
    last_size = Q->size;
  }

  /* Diagonalize (if necessary) */
  if (Q->evec_matrix_z == NULL || Q->evals_z == NULL || 
      Q->evec_matrix_inv_z == NULL) 
    mm_diagonalize(Q);

  /* Diagonalization failed: use taylor expansion instead */
  if (Q->evec_matrix_z == NULL || Q->evals_z == NULL ||
      Q->evec_matrix_inv_z == NULL) {
    mm_exp_taylor(P, Q, t);
    return;
  }
    
  /* Compute P(t) = S exp(Dt) S^-1.  Start by computing exp(Dt) S^-1 */
  for (i = 0; i < n; i++) {
    Complex exp_dt_i =
      z_exp(z_mul_real(zvec_get(Q->evals_z, i), t));
    for (j = 0; j < n; j++) 
      zmat_set(tmp, i, j, z_mul(exp_dt_i, zmat_get(Q->evec_matrix_inv_z, i, j)));
  }

  /* Now multiply by S (on the left) */
  zmat_mult_real(P->matrix, Q->evec_matrix_z, tmp);

}

/* version that assumes real eigenvalues/eigenvectors */
void mm_exp_real(MarkovMatrix *P, MarkovMatrix *Q, double t) {
  static Vector *exp_evals = NULL; /* reuse if possible */
  static int last_size = -1;
  int n = Q->size;
  int i;

 if (!(P->size == Q->size && t >= 0))
   die("ERROR mm_exp_real: got P->size=%i, Q->sizse=%i, t=%f\n",
       P->size, Q->size, t);

  if (t == 0) {
    mat_set_identity(P->matrix);
    return;
  }

  if (exp_evals == NULL || last_size != Q->size) {
    if (exp_evals != NULL) 
      vec_free(exp_evals);

    exp_evals = vec_new(Q->size);
    set_static_var((void**)&exp_evals);
    last_size = Q->size;
  }

  /* Diagonalize (if necessary) */
  if (Q->evec_matrix_r == NULL || Q->evals_r == NULL || 
      Q->evec_matrix_inv_r == NULL) 
    mm_diagonalize(Q);

  if (Q->evec_matrix_r == NULL || Q->evals_r == NULL ||
      Q->evec_matrix_inv_r == NULL) {
    mm_exp_taylor(P, Q, t);
    return;
  }

  /* Compute P(t) = S exp(Dt) S^-1 */
  for (i = 0; i < n; i++) 
    exp_evals->data[i] = exp(Q->evals_r->data[i] * t);

  mat_mult_diag(P->matrix, Q->evec_matrix_r, exp_evals, Q->evec_matrix_inv_r);
}

/* computes discrete matrix P by the formula P = exp(Qt),
   given Q and t */
void mm_exp(MarkovMatrix *dest, MarkovMatrix *src, double t) {
  if (src->eigentype == REAL_NUM)
    mm_exp_real(dest, src, t);
  else
    mm_exp_complex(dest, src, t);
}

/* given a state, draw the next state from the multinomial
 * distribution defined by the corresponding row in the matrix */
int mm_sample_state(MarkovMatrix *M, int state) {
  Vector *v = mat_get_row(M->matrix, state);
  int retval = pv_draw_idx(v);
  vec_free(v);
  return retval;
}

/* as above but by character */
char mm_sample_backgd(char *labels, Vector *backgd) {
  return labels[pv_draw_idx(backgd)];
}


/* given a state, draw the next state from the multinomial
 * distribution defined by the corresponding row in the matrix, but
 * identify both the "from" and "to" states by character label */
char mm_sample_char(MarkovMatrix *M, char c) {
  return(M->states[mm_sample_state(M, M->inv_states[(int)c])]);
}

/* WARNING: assumes matrices in dest already allocated and of correct
 * size.  Also assumes type, states, size, and eigentype are the same */
void mm_cpy(MarkovMatrix *dest, MarkovMatrix *src) {
  mat_copy(dest->matrix, src->matrix);
  if (src->eigentype == COMPLEX_NUM) {
    if (src->evec_matrix_z != NULL)
      zmat_copy(dest->evec_matrix_z, src->evec_matrix_z);
    if (src->evals_z != NULL)
      zvec_copy(dest->evals_z, src->evals_z);
    if (src->evec_matrix_inv_z != NULL)
      zmat_copy(dest->evec_matrix_inv_z, src->evec_matrix_inv_z);
  }
  else {
    if (src->evec_matrix_r != NULL)
      mat_copy(dest->evec_matrix_r, src->evec_matrix_r);
    if (src->evals_r != NULL)
      vec_copy(dest->evals_r, src->evals_r);
    if (src->evec_matrix_inv_r != NULL)
      mat_copy(dest->evec_matrix_inv_r, src->evec_matrix_inv_r);
  }
}

MarkovMatrix *mm_create_copy(MarkovMatrix *src) {
  MarkovMatrix *retval = mm_new(src->size, src->states, src->type);
  retval->eigentype = src->eigentype;
  if (src->eigentype == COMPLEX_NUM) {
    if (src->evec_matrix_z != NULL)
      retval->evec_matrix_z = zmat_new(src->size, src->size);
    if (src->evals_z != NULL)
      retval->evals_z = zvec_new(src->size);
    if (src->evec_matrix_inv_z != NULL)
      retval->evec_matrix_inv_z = zmat_new(src->size, src->size);
  }
  else {
    if (src->evec_matrix_r != NULL)
      retval->evec_matrix_r = mat_new(src->size, src->size);
    if (src->evals_r != NULL)
      retval->evals_r = vec_new(src->size);
    if (src->evec_matrix_inv_r != NULL)
      retval->evec_matrix_inv_r = mat_new(src->size, src->size);
  }
  mm_cpy(retval, src);
  return(retval);
}

void mm_diagonalize_complex(MarkovMatrix *M) { 
  if (M->evec_matrix_z == NULL)
    M->evec_matrix_z = zmat_new(M->size, M->size);
  if (M->evals_z == NULL)
    M->evals_z = zvec_new(M->size);
  if (M->evec_matrix_inv_z == NULL)
    M->evec_matrix_inv_z = zmat_new(M->size, M->size);
  if (mat_diagonalize(M->matrix, M->evals_z, M->evec_matrix_z, M->evec_matrix_inv_z)) {
    zmat_free(M->evec_matrix_z);
    zvec_free(M->evals_z);
    zmat_free(M->evec_matrix_inv_z);
    M->evec_matrix_z = NULL;
    M->evals_z = NULL;
    M->evec_matrix_inv_z = NULL;
  }
} 

void mm_diagonalize_real(MarkovMatrix *M) { 
  /* use existing routines then "cast" complex matrices/vectors as real */

  /* keep temp storage around -- this function will be called many
     times repeatedly */
  static Zmatrix *evecs_z = NULL;
  static Zmatrix *evecs_inv_z = NULL;
  static Zvector *evals_z = NULL;
  static int size = -1;

  if (evecs_z == NULL || size != M->size) {
    if (evecs_z != NULL) {
      zmat_free(evecs_z);
      zmat_free(evecs_inv_z);
      zvec_free(evals_z);
      evecs_z = NULL;
    }

    evecs_z = zmat_new(M->size, M->size);
    set_static_var((void**)&evecs_z);
    evecs_inv_z = zmat_new(M->size, M->size);
    evals_z = zvec_new(M->size);
    size = M->size;
  }

  if (1 == mat_diagonalize(M->matrix, evals_z, evecs_z, evecs_inv_z)) 
    goto mm_diagonalize_real_fail;

  if (M->evec_matrix_r == NULL) {
    M->evec_matrix_r = mat_new(M->size, M->size);
    M->evals_r = vec_new(M->size);
    M->evec_matrix_inv_r = mat_new(M->size, M->size);
  }

  if (zvec_as_real(M->evals_r, evals_z, FALSE) ||
      zmat_as_real(M->evec_matrix_r, evecs_z, FALSE) ||
      zmat_as_real(M->evec_matrix_inv_r, evecs_inv_z, FALSE))
    goto mm_diagonalize_real_fail;
  return;
  
 mm_diagonalize_real_fail:
  //by setting eigenvalues to NULL, mm_exp will call mm_exp_taylor
  //instead of using eigenvalues.
  if (M->evec_matrix_r != NULL)
    mat_free(M->evec_matrix_r);
  if (M->evals_r != NULL)
    vec_free(M->evals_r);
  if (M->evec_matrix_inv_r != NULL)
    mat_free(M->evec_matrix_inv_r);
  M->evec_matrix_r = M->evec_matrix_inv_r = NULL;
  M->evals_r = NULL;

} 

void mm_diagonalize(MarkovMatrix *M) { 
  if (M->eigentype == COMPLEX_NUM) 
    mm_diagonalize_complex(M);
  else 
    mm_diagonalize_real(M);
}

void mm_scale(MarkovMatrix *M, double scale) {
  mat_scale(M->matrix, scale);
}

/* Renormalize a discrete Markov matrix so that all rows sum to 1. */
void mm_renormalize(MarkovMatrix *M) {
  int i, j;
  if (M->type != DISCRETE)
    die("ERROR mm_renormalize:  M->type should be discrete\n");
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
