/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: misc.h,v 1.27 2008-12-10 18:09:17 agd27 Exp $ */

#ifndef MISC_H
#define MISC_H

#include <math.h>
#include <stdio.h>
#include <matrix.h>
#include <ctype.h>
#include <lists.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <external_libs.h>
struct hash_table;

#define TRUE 1
#define FALSE 0

/* for various reasons, it's often useful to represent infinity and
   negative infinity as very large numbers */
#define INFTY 999999999
#define NEGINFTY -INFTY

#define SUM_LOG_THRESHOLD -10   /* see log_sum */

#define exp2(x) (pow(2,x))
#define log2(x) ((x) <= 0 ? NEGINFTY : log(x) / M_LN2)
#define log10(x) ((x) <= 0 ? NEGINFTY : log(x) / M_LN10)
#define logit(x) ( 1 / (1 + exp(-(x))) )

/* CAREFUL: multiple eval! */
#ifndef F2C_INCLUDE
#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif

#define int_to_ptr(i) ((void*) (long) (i))
#define ptr_to_int(p) ((int) (long) (p))


#define safediv(x, y) ((y) != 0 ? (x) / (y) : ((x) == 0 ? 0 : ((x) > 0 ? INFTY : NEGINFTY)))

#define AA_ALPHABET "ARNDCQEGHILKMFPSTWYV$"

#define IS_PURINE(b) (toupper(b) == 'A' || toupper(b) == 'G')
#define IS_PYRIMIDINE(b) (toupper(b) == 'C' || toupper(b) == 'T')

/* raise integer to small integral power */
static PHAST_INLINE
int int_pow(int x, int y) { 
  int retval = 1, i;
  for (i = 0; i < y; i++) retval *= x;
  return retval;
}

/* fast computation of floor(log2(x)), where x is a positive integer */
static PHAST_INLINE
int log2_int(unsigned x) {
  int i;
  for (i = 0; ; i++) {
    x >>= 1;
    if (x == 0) return i;
  }
}

/* efficiently compute log of sum of values, which themselves are
   stored as logs: that is, return log(sum_i exp(l_i)).  The largest
   of the elements of l (call it maxval) is factored out, so that
   log(sum_i(exp(l_i))) = maxval + log(1 + sum_i(exp(l_i-maxval))),
   where the new sum is taken over 2 <= i < n.  All of the quantities
   in the exp must be negative, and those smaller than some reasonable
   threshold can be ignored. [Thanks to David Haussler for showing me
   this trick].  WARNING: sorts list as side effect. */
static PHAST_INLINE
double log_sum(List *l) {
  double maxval, expsum;
  int k;

  if (lst_size(l) > 1)
    lst_qsort_dbl(l, DESCENDING);

  maxval = lst_get_dbl(l, 0);  
  expsum = 1;
  k = 1;

  while (k < lst_size(l) && lst_get_dbl(l, k) - maxval > SUM_LOG_THRESHOLD)
    expsum += exp2(lst_get_dbl(l, k++) - maxval);
  
  return maxval + log2(expsum);        
}

/* Same as above, but base e */
static PHAST_INLINE
double log_sum_e(List *l) {
  double maxval, expsum;
  int k;

  if (lst_size(l) > 1)
    lst_qsort_dbl(l, DESCENDING);

  maxval = lst_get_dbl(l, 0);  
  expsum = 1;
  k = 1;

  while (k < lst_size(l) && lst_get_dbl(l, k) - maxval > SUM_LOG_THRESHOLD)
    expsum += exp(lst_get_dbl(l, k++) - maxval);
  
  return maxval + log(expsum);        
}

/* return n! */
static PHAST_INLINE
int permutations(int n) {
  int i, retval = 1;
  for (i = 2; i <= n; i++) retval *= i;
  return retval;
}

/* return n-choose-k */
static PHAST_INLINE 
int combinations(int n, int k) {
  int i, retval = 1;
  for (i = 0; i < k; i++) retval *= (n - i);
  return retval / permutations(k);
}

/* compute relative entropy in bits of q with respect to p, both
   probability vectors of dimension d */
static PHAST_INLINE
double rel_entropy(double *p, double *q, int d) {
  int i;
  double H = 0;
  for (i = 0; i < d; i++) {
    if (p[i] == 0) continue;
    if (q[i] == 0) return INFTY;    
    H += p[i] * (log2(p[i]) - log2(q[i]));
  }
  return H;
}

/* symmetric version of relative entropy */
static PHAST_INLINE
double sym_rel_entropy(double *p, double *q, int d) {
  double re1 = rel_entropy(p, q, d), re2 = rel_entropy(q, p, d);
  return min(re1, re2);
}
#if defined(__MINGW32__)
int random();
void srandom(int seed);
#endif
void set_seed(int seed);
void choose(int *selections, int N, int k);
void permute(int *permutation, int N);
char* get_codon_mapping(char *alphabet);
int tuple_index(char *tuple, int *inv_alph, int alph_size);
void get_tuple_str(char *tuple_str, int tuple_idx, int tuple_size, 
                   char *alphabet);
Matrix* read_subst_mat(FILE *F, char *alph);
FILE* fopen_fname(const char *fname, char *mode);
#ifdef RPHAST
#undef Rf_error
#undef die
#include <R.h>
#define die Rf_error
#define phast_warning Rf_warning
#undef printf
#define printf Rprintf
int rphast_fprintf(FILE *f, const char *format, ...);
#undef fprintf
#define fprintf rphast_fprintf
#define checkInterrupt() R_CheckUserInterrupt()
#define checkInterruptN(i, n) if ((i)%(n) == 0) R_CheckUserInterrupt()
void sfree(void *ptr);
#else
#define sfree free
void phast_warning(const char *warnfmt, ...);
void die(const char *warnfmt, ...);
#define checkInterrupt()
#define checkInterruptN(i,n)
#define unif_rand(void) (1.0*random()/RAND_MAX)
#endif
List *get_arg_list(char *arg);
List *remaining_arg_list(char *argv[], int argc, int optind);
List *get_arg_list_int(char *arg);
List *get_arg_list_dbl(char *arg);
int get_arg_int(char *arg);
double get_arg_dbl(char *arg);
int get_arg_int_bounds(char *arg, int min, int max);
double get_arg_dbl_bounds(char *arg, double min, double max);
void *smalloc(size_t size);
void *srealloc(void *ptr, size_t size);
void set_static_var(void **ptr);
char *copy_charstr(const char *word);
double normalize_probs(double *p, int size);
int is_transition(char b1, char b2);
int is_indel(char b1, char b2);
void unif_draw(int n, double min, double max, double *draws, int antithetics);
int bn_draw(int N, double p);
int bn_draw_fast(int n, double pp);
void mn_draw(int n, double *p, int d, int *counts);
int draw_index(double *p, int size);
struct hash_table *make_name_hash(char *mapstr);
double exp_draw(double b);
double gamma_pdf(double x, double a, double b);
double gamma_cdf(double x, double a, double b, int lower_tail);
double gamma_draw(double a, double b);
double chisq_pdf(double x, double dof);
double chisq_cdf(double x, double dof, int lower_tail);
double half_chisq_cdf(double x, double dof, int lower_tail);
double d_beta(double x, double a, double b);
double beta_draw(double a, double b);
void dirichlet_draw(int k, double *alpha, double *theta);
int next_comb(int n, int k, int *index);
double incomplete_gamma(double a, double x, char type);
double d_poisson(double lambda, int k);
double cum_poisson(double lambda, int k);
double cum_poisson_c(double lambda, int k);
double cum_norm(double mu, double sigma, double a);
double cum_norm_c(double mu, double sigma, double a);
double inv_cum_norm(double p);
double bvn_p(double x, double y, double mu_x, double mu_y, double sigma_x,
             double sigma_y, double rho);
void norm_confidence_interval(double mu, double sigma, double interval_size, 
                              double *min_x, double *max_x);
void print_seq_fasta(FILE *F, char *seq, char *name, int len);
double get_elapsed_time(struct timeval *start_time);
int file_exists(char *filename);
char **get_iupac_map();
int **build_iupac_inv_map(int *inv_states, int alph_size);
void free_iupac_inv_map(int **iim);
void *alloc_n_dimensional_array(int ndim, int *dimsize, size_t size);
void free_n_dimensional_array(void *data, int ndim, int *dimsize);
#endif
