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
#include <external_libs.h>

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
  M->diagonalize_error = -1;
  M->matrix = mat_new(size, size);
  mat_zero(M->matrix);
  M->size = size;
  alph_size = states == NULL ? size : (int)strlen(states);
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
  M->diagonalize_error = -1;
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
    if (fabs(sum + diag - targetval) > SUM_EPSILON) {
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


//internal function used by mm_exp_higham; computes some even powers
//of Qt and returns them in an array.  Qt should be square.
Matrix **mm_get_QtPow(int max_m, Matrix *Qt) {
  Matrix **QtPow = smalloc((max_m+1) * sizeof(Matrix*));
  int i, n = Qt->nrows;
  for (i=0; i<=max_m; i++) {
    if (i==1 || i%2==0) QtPow[i] = mat_new(n, n);
    else QtPow[i] = NULL;
    if (i==0)
      mat_set_identity(QtPow[0]);
    else if (i==1)
      mat_copy(QtPow[1], Qt);
    else if (i==2)
      mat_mult(QtPow[2], Qt, Qt);
    else if (i==4)
      mat_mult(QtPow[4], QtPow[2], QtPow[2]);
    else if (i==6)
      mat_mult(QtPow[6], QtPow[2], QtPow[4]);
    else if (i==8)
      mat_mult(QtPow[8], QtPow[4], QtPow[4]);
    else if (i==10)
      mat_mult(QtPow[10], QtPow[4], QtPow[6]);
    else if (i==12)
      mat_mult(QtPow[12], QtPow[6], QtPow[6]);
  }
  return QtPow;
}

//set P->matrix to exp(t * Q->matrix) using algorithm 2.3 from
//"The Scaling and Squaring Method for the Matrix Exponential Revisited"
//by Nicholas J. Higham.
void mm_exp_higham(MarkovMatrix *P, MarkovMatrix *Q, double t, int do_mu) {
#ifdef SKIP_LAPACK
  die("ERROR: LAPACK required for mm_exp_higham routine.\n");
#else
  double b[14] = {64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
		  1187353796428800.0, 129060195264000.0, 10559470521600.0,
		  670442572800.0, 33522128640.0, 1323241920.0,
		  40840800.0, 960960.0, 16380.0, 182.0, 1.0};
  double theta[14], norm, sum, *coef, mu=0.0, val;
  int i, j, n = P->size, mvals[5] = {3, 5, 7, 9, 13}, mlen=5, m, mi, max_i, s;
  LAPACK_DOUBLE mat[n*n], scale[n], matU[n*n], matV[n*n];
  LAPACK_INT ln, ilo, ihi, info, ipiv[n];

  char job = 'B', side='R';
  Matrix *Qt = mat_create_copy(Q->matrix), **QtPow, **matArr, *U, *V;
  int do_balancing=0;  //balancing is not working yet, keep this at 0!
  mat_scale(Qt, t);

  theta[3]=1.495585217958292e-2;
  theta[5]=2.539398330063230e-1;
  theta[7]=9.504178996162932e-1;
  theta[9]=2.097847961257068e0;
  theta[13]=5.371920351148152e0;

  // % Preprocessing to reduce the norm.
  //A ← A − μI, where μ = trace(A)/n.
  if (do_mu) {
    mu = 0.0;
    for (i=0; i < n; i++)
      mu += mat_get(Qt, i, i);
    mu /= (double)n;
    for (i=0; i < n; i++)
      mat_set(Qt, i, i, mat_get(Qt, i, i) - mu);
  }

  //A ← D −1 AD, where D is a balancing transformation (or set D = I if
  //balancing does not reduce the 1-norm of A).
  ln = (LAPACK_INT)n;
  if (do_balancing) {
    mat_to_lapack(Qt, mat);
    ln = (LAPACK_INT)n;
#ifdef R_LAPACK
    F77_CALL(dgebal)(&job, &ln, mat, &ln, &ilo, &ihi, scale, &info);
#else
    dgebal_(&job, &ln, mat, &ln, &ilo, &ihi, scale, &info);
#endif
    if (info != 0)
      die("Error in balancing matrix in lapack routine dgebal info=%i\n", info);
    mat_from_lapack(Qt, mat);
  }

  //compute 1-norm of Qt
  norm = 0.0;
  for (i=0; i < n; i++) {
    sum = 0.0;
    for (j=0; j < n; j++)
      sum += fabs(mat_get(Qt, j, i));
    if (sum > norm) norm = sum;
  }
  U = mat_new(n, n);
  V = mat_new(n, n);
  s = 0;

  for (mi=0; mi < mlen; mi++) {
    m = mvals[mi];
    if (norm <= theta[m]) {
      QtPow = mm_get_QtPow(m-1, Qt);
      max_i = (m-1)/2;
      mat_zero(P->matrix);
      matArr = smalloc((max_i+1)*sizeof(Matrix*));
      coef = smalloc((max_i+1)*sizeof(double));
      for (i=0; i <= max_i; i++) {
	matArr[i] = QtPow[2*i];
	coef[i] = b[2*i+1];
      }
      mat_linear_comb_many(Qt, max_i+1, matArr, coef);
      mat_mult(U, QtPow[1], Qt);

      for (i=0; i<=max_i; i++)
	coef[i] = b[2*i];
      mat_linear_comb_many(V, max_i+1, matArr, coef);
      goto mm_exp_higham_solve_linear_eq;
    }
  }

  s = (int)ceil(log2(norm/theta[13]));
  if (s > 0) mat_scale(Qt, 1.0/pow(2.0, s));

  m=7;
  QtPow = mm_get_QtPow(m, Qt);
  matArr = smalloc(5 * sizeof(Matrix*));
  coef = smalloc(5 * sizeof(double));

  //calculate U
  //b13*A[6] + b11*A[4] + b9*A[2]
  U = mat_new(n, n);
  matArr[0] = QtPow[6];
  coef[0] = b[13];
  matArr[1] = QtPow[4];
  coef[1] = b[11];
  matArr[2] = QtPow[2];
  coef[2] = b[9];
  mat_linear_comb_many(Qt, 3, matArr, coef);

  //multiply by A[6]
  mat_mult(U, QtPow[6], Qt);

  //above + b7*A[6] + b5*A[4] + b3*A[2] + b1*I
  matArr[3] = QtPow[0];
  matArr[4] = U;
  coef[0] = b[7];
  coef[1] = b[5];
  coef[2] = b[3];
  coef[3] = b[1];
  coef[4] = 1.0;
  mat_linear_comb_many(Qt, 5, matArr, coef);

  //multiply whole thing by QtPow[1]
  mat_mult(U, QtPow[1], Qt);

  //calculate V
  //b12*A[6] + b[10]*A[4] + b8*A[2]
  V = mat_new(n, n);
  coef[0] = b[12];
  coef[1] = b[10];
  coef[2] = b[8];
  mat_linear_comb_many(Qt, 3, matArr, coef);

  //multiply by A[6]
  mat_mult(V, QtPow[6], Qt);

  //above + b[6]*A[6] + b[4]*A[4] + b[2]*A[2] + b0*I
  mat_copy(Qt, V);
  matArr[4] = Qt;
  coef[0] = b[6];
  coef[1] = b[4];
  coef[2] = b[2];
  coef[3] = b[0];
  coef[4] = 1.0;
  mat_linear_comb_many(V, 5, matArr, coef);

 mm_exp_higham_solve_linear_eq:
  //set U' = -U + V
  //set V' = U + V
  //solve U'X = V' for X
  mat_copy(Qt, U);
  mat_plus_eq(V, U);
  mat_copy(U, V);
  mat_scale(Qt, -2.0);
  mat_plus_eq(U, Qt);

  sfree(matArr);
  sfree(coef);

  mat_to_lapack(U, matU);
  mat_to_lapack(V, matV);
#ifdef R_LAPACK
  F77_CALL(dgesv)(&ln, &ln, matU, &ln, ipiv, matV, &ln, &info);
#else
  dgesv_(&ln, &ln, matU, &ln, ipiv, matV, &ln, &info);
#endif
  if (info !=0)
    die("Error solving U'X=V' in mm_exp_higham");
  mat_from_lapack(P->matrix, matV);

  //check (for testing only)
  mat_mult(Qt, U, P->matrix);

  mat_free(U);
  mat_free(V);

  // X = r13^(2^s) by repeated squaring
  for (i=0; i < s; i++) {
    mat_copy(Qt, P->matrix);
    mat_mult(P->matrix, Qt, Qt);
  }

  if (do_balancing) {
    mat_to_lapack(P->matrix, mat);
    // undo balancing- check- this function is meant for eigenvalues,
    // not sure if job should be 'R' or 'L', or if it will work at all.
#ifdef R_LAPACK
    F77_CALL(dgebak)(&job, &side, &ln, &ilo, &ihi, scale, &ln, mat, &ln, &info);
#else
    dgebak_(&job, &side, &ln, &ilo, &ihi, scale, &ln, mat, &ln, &info);
#endif
    mat_from_lapack(P->matrix, mat);
  }

  if (do_mu)
    mat_scale(P->matrix, exp(mu));

  //check
  for (i=0; i < n; i++) {
    double sum = 0.0;
    for (j=0; j < n; j++) {
      val = mat_get(P->matrix, i, j);
      if (isinf(val) || isnan(val))
	die("got nan/inf val in matrix in mm_exp_higham t=%f i=%i j=%i\n", t, i, j);
      else if (val < 0 && val > -1.0e-10)
	mat_set(P->matrix, i, j, 0.0);
      else if (val > 1 && val-1 < 1.0e-10)
	mat_set(P->matrix, i, j, 1.0);
      else if (val < 0 || val > 1)
	die("got val %e at pos %i %i in mm_exp_higham", val, i, j);
      sum += val;
    }
    if (fabs(sum-1.0) > 1.0e-4)
      die("got sum %e at pos %i in mm_exp_higham", sum, i);
  }

  mat_free(Qt);
  for (i=0; i < m; i++)
    if (QtPow[i] != NULL) mat_free(QtPow[i]);
  sfree(QtPow);
#endif
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
  if (Q->diagonalize_error != 1 &&
      (Q->evec_matrix_z == NULL || Q->evals_z == NULL ||
       Q->evec_matrix_inv_z == NULL))
    mm_diagonalize(Q);

  /* Diagonalization failed: use higham expansion instead */
  if (Q->evec_matrix_z == NULL || Q->evals_z == NULL ||
      Q->evec_matrix_inv_z == NULL) {
    mm_exp_higham(P, Q, t, 1);
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
  if (Q->diagonalize_error != 1 &&
      (Q->evec_matrix_r == NULL || Q->evals_r == NULL ||
       Q->evec_matrix_inv_r == NULL))
    mm_diagonalize(Q);

  if (Q->evec_matrix_r == NULL || Q->evals_r == NULL ||
      Q->evec_matrix_inv_r == NULL) {
    mm_exp_higham(P, Q, t, 1);
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
  if (strlen(labels) != backgd->size) die("mm_sample_backgd: got num_labels=%i but backgd->size=%i\n", strlen(labels), backgd->size);
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
    M->diagonalize_error = 1;
  }
  else M->diagonalize_error = 0;
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

  M->diagonalize_error = 0;
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
  //by setting eigenvalues to NULL, mm_exp will call mm_exp_higham
  //instead of using eigenvalues.
  if (M->evec_matrix_r != NULL)
    mat_free(M->evec_matrix_r);
  if (M->evals_r != NULL)
    vec_free(M->evals_r);
  if (M->evec_matrix_inv_r != NULL)
    mat_free(M->evec_matrix_inv_r);
  M->evec_matrix_r = M->evec_matrix_inv_r = NULL;
  M->evals_r = NULL;
  M->diagonalize_error = 1;
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
