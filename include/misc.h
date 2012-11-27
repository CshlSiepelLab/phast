/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/**
   @file misc.h
   Miscellaneous definitions and functions: Argument handling, PDF/CDF, Distribution draws, Codon mapping, etc.
   @ingroup base
*/

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

/** Infinity, for various reasons, it's often useful to represent infinity and
   negative infinity as very large numbers */
#define INFTY 999999999
/** Negative Infinity, for various reasons, it's often useful to represent infinity and negative infinity as very large numbers */
#define NEGINFTY -INFTY

/** Threshold for log_sum function
    @see log_sum
 */
#define SUM_LOG_THRESHOLD -50  

/** Shortcut for 2^x */
#define exp2(x) (pow(2,x))
/** Log base 2.  Negative input results in -inf */
#define log2(x) ((x) <= 0 ? NEGINFTY : log(x) / M_LN2)
/** Log base 10.  Negative input results in -inf */
#define log10(x) ((x) <= 0 ? NEGINFTY : log(x) / M_LN10)
/** Logit */
#define logit(x) ( 1 / (1 + exp(-(x))) )

/* CAREFUL: multiple eval! */
#ifndef F2C_INCLUDE
/** Max (if we don't have f2c) */
#define max(x, y) ((x) > (y) ? (x) : (y))
/** Min (if we don't have f2c) */
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif

/** Convert integer to pointer */
#define int_to_ptr(i) ((void*) (long) (i))
/** Convert pointer to integer */
#define ptr_to_int(p) ((int) (long) (p))

/** Safe divide, checks for div by 0 so no arithmetic errors are thrown */
#define safediv(x, y) ((y) != 0 ? (x) / (y) : ((x) == 0 ? 0 : ((x) > 0 ? INFTY : NEGINFTY)))

/** Amino Acid alphabet */
#define AA_ALPHABET "ARNDCQEGHILKMFPSTWYV$"

/** Test if a base is Purine */
#define IS_PURINE(b) (toupper(b) == 'A' || toupper(b) == 'G')
/** Test if a base is Pyrimidine */
#define IS_PYRIMIDINE(b) (toupper(b) == 'C' || toupper(b) == 'T')

/** Raise integer to small integral power
    @param x Integer to raise
    @param y Power to raise to
    @result x^y
 */
static PHAST_INLINE
int int_pow(int x, int y) { 
  int retval = 1, i;
  for (i = 0; i < y; i++) retval *= x;
  return retval;
}

/** \name Log based calculation functions 
\{ */

/** Fast computation of floor(log2(x)), where x is a positive integer 
    @param x Integer to take log2 of then floor
*/
static PHAST_INLINE
int log2_int(unsigned x) {
  int i;
  for (i = 0; ; i++) {
    x >>= 1;
    if (x == 0) return i;
  }
}

/** Efficiently compute log of sum of values.  
   @param l List of doubles containing values
   @result log of sum of values passed in list
   @note Thanks to David Haussler for showing me this trick.  
   @warning Sorts list as side effect. 
*/
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

/** Efficiently compute log (base e) of sum of values.  
   @param l List of doubles containing values
   @result log (base e) of sum of values passed in list
   @note Thanks to David Haussler for showing me this trick.  
   @warning Sorts list as side effect. 
*/
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

/** \} */

/** Compute relative entropy in bits of q with respect to p, both
   probability vectors of dimension d.
   @param q Probability distribution vector to compute relative entry of
   @param p Probability distribution vector with respect of
   @param d Dimension of probability vectors
 */
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

/** Compute symmetric relative entropy in bits of q with respect to p, both probability vectors of dimension d.
   @param q Probability distribution vector to compute relative entry of
   @param p Probability distribution vector with respect of
   @param d Dimension of probability vectors
 */
static PHAST_INLINE
double sym_rel_entropy(double *p, double *q, int d) {
  double re1 = rel_entropy(p, q, d), re2 = rel_entropy(q, p, d);
  return min(re1, re2);
}
#if defined(__MINGW32__)
int random();
void srandom(int seed);
#endif
/** Specify the Random Number Generator seed number.
    @param seed Starting number for RNG
 */
void set_seed(int seed);

/** \name Combination & Permutation functions
\{ */


/** Return n! */
static PHAST_INLINE
int permutations(int n) {
  int i, retval = 1;
  for (i = 2; i <= n; i++) retval *= i;
  return retval;
}

/** Return n-choose-k */
static PHAST_INLINE 
int combinations(int n, int k) {
  int i, retval = 1;
  for (i = 0; i < k; i++) retval *= (n - i);
  return retval / permutations(k);
}


/** Randomly choose k elements from a list of N.
    @param[in,out] selections Result array of size N to be populated with 0 (not chosen) or 1 (chosen). Elements initialized to -1 will be consitered "ineligable" and skipped.
    @param[in] N Number of elements to choose from
    @param[in] k Number of elements to choose
 */ 
void choose(int *selections, int N, int k);

/** Next combination (used for enumerating combinations).
    Call repeatedly to enumerate combinations. 
    @param n Number of elements
    @param k Size of index array
    @param index Array preallocated to size k used to hold output. At end, contains indices in [0, n-1] representing the next element in the series of all possible combinations of n elements. On the first call, set indices[0] = -1 and the array will be initialized appropriately.  
   @result TRUE on success, FALSE when no more choices possible
*/
int next_comb(int n, int k, int *index);

/** Produce a random permutation of the designated size.
    @param permutation Result array of size N to be populated with 0 to N-1 in a random order
    @param N Number of elements in the permutation
*/
void permute(int *permutation, int N);

/** \} */

/** Create map from Codons to Amino Acids.
    @param alphabet Possible characters in a sequence (assumed to contain the characters 'A', 'C', 'G',
   and 'T')
    @result Result char array of size (alphabet size)^3 containing codons mapped to corresponding amino acids according to the universal genetic code.
    @note Stop codons will be mapped to '$' characters
    @note Tuples not in {A,C,G,T}^3 will be mapped to null characters ('\0')
    @code
     char *alph = "ACGT";
     get_codon_mapping(alph);
    
 
//    Result  K N K N T TTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV$Y$YSSSS$CWCLFLF
//            | | | | | |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//    Codons  a a a a a aaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt
//            a a a a c cccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt
//            a c g t a cgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt
//      
//    AAA=K, AAC=N, AAG=K, AAT=N, ACA=T.... etc.

    @endcode
*/
char* get_codon_mapping(char *alphabet);

/** Translate a string tuple (for example "AAC") to an array index (43) for the 'codons' array used in get_codon_mapping.
   @param[in] tuple Three characters to map to array index
   @param[in] inv_alph Inverse alphabet that maps chars (as int) to indices of alphabet array
   @param[in] alph_size Size of the alphabet
   @result Array index used for 'codons' array OR -1 if char not in alphabet encountered
   @note A "digital" indexing scheme is assumed in which the right-most character in 
         the tuple is considered the least significant digit.
   @note Inverse of get_tuple_str
   @see get_tuple_str
   @see get_codon_mapping
  */
int tuple_index(char *tuple, int *inv_alph, int alph_size);

/** Translate an array index (for example 43) from the 'codons' array used in get_codon_mapping to a string tuple ("AAC")
   @param[out] tuple_str String corresponding to tuple index
   @param[in] tuple_idx Index of tuple in Codons array
   @param[in] tuple_size Amount of chars in a tuple
   @param[in] alphabet Valid characters in a tuple
   @note Inverse of tuple_index
   @see tuple_index
*/
void get_tuple_str(char *tuple_str, int tuple_idx, int tuple_size, 
                   char *alphabet);
/** Read a substitution matrix from the specified file.
   @param F File descriptor of file containing substitution matrix
   @param alph (Optional)  If 'alph' equals the empty string (must be preallocated to adequate
               size), then the sequence of characters that defines the matrix rows
               and columns will be stored in it; otherwise, 'alph' will be taken
               to define the desired alphabet, and the order of rows and columns
               will be rearranged as necessary to be consistent with it
   @note format expected to be that used by BLAST, as output by the NCBI "pam" program
   @note Characters in the file but not in 'alph' will be ignored
   @warning Only minor testing performed.  Check carefully if use
   without predefined alphabet or with alphabets that do not match in
   order 

*/
Matrix* read_subst_mat(FILE *F, char *alph);

/** Open a file by filename and get file descriptor.
    @param fname Full path to file
    @param mode Open mode i.e. w, r, r+, w+, etc.
    @result File descriptor 
    @note Exits with error message if open unsuccessful
 */
FILE* phast_fopen(const char *fname, char *mode);

/** Open a file by filename and get file descriptor.
    @param fname Full path to file
    @param mode Open mode i.e. w, r, r+, w+, etc.
    @result File descriptor 
 */
FILE* phast_fopen_no_exit(const char *fname, char *mode);


void phast_fclose(FILE *f);
#ifdef RPHAST
#undef Rf_error
#undef die
#include <R.h>
#define die Rf_error
#define phast_warning Rf_warning
#undef printf
#define printf Rprintf
#undef stdout
#undef stderr
extern FILE *rphast_stdout;
extern FILE *rphast_stderr;
#define stdout rphast_stdout
#define stderr rphast_stderr

/** Write text to file or R stdout/stderr.
    @param f File descriptor of file to write to OR stdout to write to R console OR stderr to write to R error console.
    @param format Format of the string to write like sprintf
    @param ... Parameters like sprintf
    @result Success == 1  
  */
int rphast_fprintf(FILE *f, const char *format, ...);
#undef fprintf
#define fprintf rphast_fprintf
#define checkInterrupt() R_CheckUserInterrupt()
#define checkInterruptN(i, n) if ((i)%(n) == 0) R_CheckUserInterrupt()

#else 

/** Display a warning on the console.
    @param warnfmt Format of the string to write to console (like printf)
*/
void phast_warning(const char *warnfmt, ...);

/** Display info about unrecoverable error to console and end program.
    @param warnfmt Format of the string to write to console (like printf)
*/
void die(const char *warnfmt, ...);
#define checkInterrupt()
#define checkInterruptN(i,n)
#define unif_rand(void) (1.0*random()/RAND_MAX)
#endif

/** \name Program argument handling functions
\{ */
/** Parse comma separated value string or file reference into list.
   @param arg Comma separated value string or file reference (using the "*" convention)
   @result List of values from arg
   @note: List and all Strings are newly allocated (should be freed externally) 
*/
List *get_arg_list(char *arg);

/** Returns remaining command-line arguments as a List of Strings. 
  @param argv As passed to main
  @param argc As passed to main
  @param optind Index of first unprocessed argument
  @result List of strings containing unprocessed command line arguments
*/
List *remaining_arg_list(char *argv[], int argc, int optind);

/** Parse comma separated value string or file reference into list of integers.
    @param arg CSV list or file reference
    @result List of integers
*/
List *get_arg_list_int(char *arg);

/** Parse comma separated value string or file reference into list of integers.
    @param arg CSV list or file reference
    @result List of doubles
*/
List *get_arg_list_dbl(char *arg);

/** Argument conversion with error checking (int)
   @param arg String to convert to int
   @result Int parsed from string
 */
int get_arg_int(char *arg);

/** Argument conversion with error checking (double)
   @param arg String to convert to double
   @result Double parsed from string
 */
double get_arg_dbl(char *arg);

/** Argument conversion with error and bounds checking (int) 
   @param arg String to convert to int
   @param min Minimum acceptable int
   @param max Maximum acceptable int
   @result Int parsed from string
*/
int get_arg_int_bounds(char *arg, int min, int max);

/** Argument conversion with error and bounds checking (double) 
   @param arg String to convert to double
   @param min Minimum acceptable double
   @param max Maximum acceptable dobule
   @result Double parsed from string
*/
double get_arg_dbl_bounds(char *arg, double min, double max);

/** \} */

/** Safe malloc 
    @param size Size of memory to allocate
*/
void *smalloc(size_t size);

/** Safe re-malloc 
    @param ptr Pointer to memory to reallocate
    @param size New size
*/
void *srealloc(void *ptr, size_t size);

void register_open_file(FILE *F);
void unregister_open_file(FILE *F);

#ifdef USE_PHAST_MEMORY_HANDLER
/** Safe memory free
    @param ptr Pointer to object to free
*/
void sfree(void *ptr);
#else
#define sfree free
#endif

void set_static_var(void **ptr);

/** Copy a string
   @param word String to copy
   @result Word coppied by value
*/
char *copy_charstr(const char *word);

/** Normalize a probability vector.  
   @pre All values nonnegative.  
   @param p Probability vector to normalize
   @param size Number of elements in vector p
   @result Normalization constant 
*/
double normalize_probs(double *p, int size);

/** Test if b1->b2 is a transition.
   @param b1 DNA base to start with
   @param b2 DNA base to end with
   @result 1 if a change from b1 to b2 is a transition, otherwise 0
 */
int is_transition(char b1, char b2);

/** Test if b1->b2 is an indel.
   @param b1 DNA base to start with
   @param b2 DNA base to end with
   @result 1 if a change from b1 to b2 involves an insertion/deletion, otherwise 0
 */
int is_indel(char b1, char b2);


/** Create HashTable from String. 
    Parse string defining mapping from old names to new and store as hash  
    @param Input String to translate to HashMap e.g. "hg17->human; mm5->mouse; rn3->rat" 
    @result HashTable containing data parsed from string
    @note Can use "=" or "->" characters to indicate mapping. 
*/
struct hash_table *make_name_hash(char *mapstr);

/** \name PDF/CDF functions 
\{ */

/** Evaluate PDF of gamma distribution with parameters a and b
    @param x Gamma distributed variable
    @param a Shape parameter of gamma distribution
    @param b Rate parameter
    @result Evaluation of PDF of Gamma dist.
*/
double gamma_pdf(double x, double a, double b);

/** Evaluate CDF of gamma distribution with parameters a and b
    @param x Gamma distributed variable
    @param a Shape parameter of gamma distribution
    @param b Rate parameter of gamma distribution
    @param lower_tail If == TRUE returns P(X<=x), else returns P(X>=x)
    @result evaluation of CDF of Gamma dist.
*/
double gamma_cdf(double x, double a, double b, int lower_tail);

/** Evaluate pdf of chi-square distribution with dof degrees of freedom.
    @param x Chi-Square distributed variable
    @param dof Degrees of freedom
    @result Evaluation of pdf
*/
double chisq_pdf(double x, double dof);

/** Evaluate cdf of chi-square distribution with dof degrees of freedom.
    @param x Chi-Square distributed variable
    @param dof Degrees of freedom
    @param lower_tail If lower_tail == TRUE returns P(X<=x), else returns P(X>=x)
    @result Evaluation of cdf
*/
double chisq_cdf(double x, double dof, int lower_tail);

/** Evaluate cdf of a 50:50 mixture of a chisq distribution with dof
   degrees of freedom and a point mass at 0.
   @param x Chi-Square distributed variable
   @param dof Degrees of freedom
   @param lower_tail If lower_tail == TRUE returns P(X<=x), else returns P(X>=x)  
   @note This function is useful in likelihood ratio tests of bounded parameters. 
*/
double half_chisq_cdf(double x, double dof, int lower_tail);

/** \} \name Density evaluation functions
\{ */

/** Evaluate density of Beta distribution. 
    @param x Beta distributed variable
    @param a Shape parameter
    @param b Shape parameter
    @result Density of Beta distribution
*/
double d_beta(double x, double a, double b);

/** Evaluate density of poisson distribution.
    @param lambda Expected number of occurrences during a given interval.
    @param k Number of occurrences of an event during a given interval we want to find the probability of
    @result Probability that k occurrences happen in a given interval, when we expect lambda occurrences during that interval
*/
double d_poisson(double lambda, int k);


/** Evaluate density of bivariate normal, p(x, y | mu_x, mu_y, sigma_x, sigma_y, rho).
   @param mu_x Marginal mean
   @param mu_y Marginal mean
   @param sigma_x Marginal standard deviation
   @param sigma_y Marginal standard deviations
   @param rho Correlation coefficient 
   @result Estimate of bivariate normal density
*/
double bvn_p(double x, double y, double mu_x, double mu_y, double sigma_x,
             double sigma_y, double rho);

/** \} \name Distribution Draw functions 
\{ */

/** Make a draw from a beta distribution with parameters 'a' and 'b'.
    @param a Shape parameter
    @param b Shape parameter
    @result Draw from distribution
*/
double beta_draw(double a, double b);

/** Make a draw from a k-dimensional Dirichlet distribution.
   @param k Dimensionality of distribution
   @param theta Scale of distribution
   @param alpha Distribution parameters given by alpha[0], ..., alpha[k-1].  
   @result Draw from distribution
   @note This is accomplished by sampling from gamma 
   distributions with parameters alpha[0], ..., alpha[k-1]
   and renormalizing.
 */
void dirichlet_draw(int k, double *alpha, double *theta);


/** Make 'n' draws from a uniform distribution on the interval [min,
   max], optionally with antithetics.    
   @pre Call srandom
   @param[in] n Number of draws to make from uniform distribution
   @param[in] min Low end of interval to make draws on
   @param[in] max High end of interval to make draws on
   @param[out] draws Resulting draws
   @param[in] antithetics (Optional) If == 1 use antithetics
   @note Designed for use with real (floating-point) numbers.  
 */
void unif_draw(int n, double min, double max, double *draws, int antithetics);

/** Make a draw from a binomial distribution with parameters 'N' and
   'p'.  
   @pre Call srandom
   @param N Amount of samples in binomial distribution
   @param p Probability of each sample yielding success
   @result Number of successful samples
   @note Computational complexity is O(N) */
int bn_draw(int N, double p);

/** Make a draw from a binomial distribution with parameters 'n' and
   'pp'.  
   @pre Call srandom
   @param N Amount of samples in binomial distribution
   @param pp Probability of occurrence of the event???
   @result Number of successful samples???
   @note Computational complexity is O(N) 
   @note This version uses rejection sampling unless n < 25 or n * p <
   1/25.   */
int bn_draw_fast(int n, double pp);

/** Make 'n' draws from a multinomial distribution defined by
   probability vector 'p' with dimension 'd'. 
   @pre Call srandom
   @param n Number of draws to make from multinomial distribution
   @param p Probability vector defining multinomial distribution
   @param d Dimensionality of multinomial distribution
   @param counts Number of counts for each category
   @note Sum of elements in 'counts' will equal 'n'.  
   @note Probability vector is assumed to be normalized.
 */
void mn_draw(int n, double *p, int d, int *counts);

/** Make a draw from an exponential distribution with parameter
   (expected value) 'b' 
   @param b Defines exponential distribution
   @result Draw from exponential distribution
*/
double exp_draw(double b);


/** Make a draw from a gamma distribution with parameters 'a' and 'b'.
   @pre Call srandom
   @param a Shape parameter of gamma distribution; If a == 1, exp_draw is
   called.  If a > 1, Best's (1978) rejection algorithm is used, and
   if a < 1, rejection sampling from the Weibull distribution is
   performed
   @param b Rate parameter of gamma distribution
   @result Draw from distribution
   @note Both algorithms used as described in "Non-Uniform 
   Random Variate Generation" by Luc Devroye, available 
   online at http://cg.scs.carleton.ca/~luc/rnbookindex.html
*/
double gamma_draw(double a, double b);

/** Given a probability vector, draw an index. 
   @pre Call srandom externally
   @param p Probability vector
   @param size Number of elements in probability vector
   @result Index drawn from probability vector
*/
int draw_index(double *p, int size);


/** \} */



/** Incomplete Gamma function.
    @param a Shape of integrand
    @param x End of integral if type='p' or Beginning of integral if type='q'
    @param type 'p' means first half of integral, 'q' means second half (see p. 171)
    @result Integration result
 */
double incomplete_gamma(double a, double x, char type);

/** \name Cumulative distribution probability functions
 \{ */

/** Evaluate P(x <= k | lambda) (first 1/2 of integral)
    Return P(x <= k | lambda), for a variable x that obeys a Poisson
   distribution with parameter lambda
    @param lambda Expected number of occurrences during a given interval.
    @param k Number of occurrences of an event during a given interval we want to find the probability of
    @result Probability that k occurrences happen in or prior to a given interval, when we expect lambda occurrences during that interval
*/
double cum_poisson(double lambda, int k);

/** Evaluate P(x <= k | lambda) (last 1/2 of integral)
    Return P(x <= k | lambda), for a variable x that obeys a Poisson
   distribution with parameter lambda
    @param lambda Expected number of occurrences during a given interval.
    @param k Number of occurrences of an event during a given interval we want to find the probability of
    @result Probability that k occurrences happen in or prior to a given interval, when we expect lambda occurrences during that interval
*/
double cum_poisson_c(double lambda, int k);

/** Return P(x <= a | mu, sigma) for a variable a that obeys a normal ???
   distribution with mean mu and s.d. sigma. 
   @param mu Mean of distribution
   @param sigma Standard deviation of distribution
   @param a Number of occurrences in a given interval
   @result Probability that 'a' occurrences will happen within or prior to a given interval
*/
double cum_norm(double mu, double sigma, double a);

/** Return P(x >= a | mu, sigma) for a variable a that obeys a normal ???
   distribution with mean mu and s.d. sigma.  
   Use this function instead of 1-cum_norm when a is large (better precision) 
   @param mu Mean of distribution
   @param sigma Standard deviation of distribution
   @param a Number of occurrences in a given interval
   @result Probability that 'a' occurrences happen in or after a given interval
*/
double cum_norm_c(double mu, double sigma, double a);

/** Return inverse of standard normal, i.e, inv_cum_norm(p) = a such
   that cum_norm(0, 1, a) = p.  
   @param p Probability that 'a' inquiries happen during an interval ???
   @result Number of occurrences that relate to the specified probability
   @note The function is approximated using an algorithm by Peter Acklam
    given at http://home.online.no/~pjacklam/notes/invnorm/. 
*/
double inv_cum_norm(double p);
/** \} */


/** Compute min and max of (central) confidence interval.
   @param mu Mean of normal distribution
   @param sigma Standard deviation of normal distribution
   @param interval_size Size of confidence interval (between 0 and 1) 
   @param min_x Min of central confidence interval
   @param max_x Max of central confidence interval
  */
void norm_confidence_interval(double mu, double sigma, double interval_size, 
                              double *min_x, double *max_x);

/** Save sequence in file as fasta format.
    @param F File descriptor to save sequence with
    @param seq Sequence chars to write out
    @param name Name of the sequence   
    @param len Length of the sequence chars
*/
void print_seq_fasta(FILE *F, char *seq, char *name, int len);

/** Return time in seconds since start_time.
    @param timeval Struct used to get current time value
    @param start_time Time in the past that an event started at
    @result Number of seconds since the time defined in start_time
*/
double get_elapsed_time(struct timeval *start_time);

/** Check to see if a file is present and readable on the file system
  @param filename path of file to check for
  @result 1 if exists and readable
*/
int file_exists(char *filename);


/** \name IUPAC mapping functions 
\{ */

/** Return a static mapping from IUPAC ambiguity characters to the bases
   that they represent.
   @code
   e.g.
   retval['R'] = "AG";
   retval['Y'] = "CT";
   @endcode
   @result Ambiguity characters -> bases map
*/
char **get_iupac_map();

/** Build an inverse mapping that allows an IUPAC ambiguity character
   to be mapped to an array, alph_size elements, with 1s or 0s
   indicating presence or absence of each character in the set
   associated with the IUPAC character.  This array obeys the indexing
   of the provided inv_states.  For use in computing likelihoods with
   ambiguity characters 
   @param inv_states Inverse states array mapping states to unique numbers
   @param alph_size Number of entries in inv_states (size of alphabet)
   @code
   e.g. 
      inv_states['A'] = 0
      inv_states['T'] = 1
      inv_states['G'] = 2
      inv_states['C'] = 3

      alph_size = 4
                       ATGC
      retval[(int)'R']  =  '1010'
      retval[(int)'Y']  =  '0101'
   @endcode
*/
int **build_iupac_inv_map(int *inv_states, int alph_size);

/** Free inverse iupac mapping object. 
    @param iim Inverse IUPAC mapping object to free.
*/
void free_iupac_inv_map(int **iim);

/** \} \name 'n' Dimensional array functions
\{  */

/** Allocate an 'n' dimensional array.
    @param ndim Number of dimensions in array
    @param dimsize Array of int describing size of each dimension
    @param Size of 1 dimensional array
    @result Pointer to 'n' dimensional array
*/
void *alloc_n_dimensional_array(int ndim, int *dimsize, size_t size);

/** Free an 'n' dimensional array.
    @param Pointer to 'n' dimensional array
    @param ndim Number of dimensions in array
    @param dimsize Array of int describing size of each dimension
 */
void free_n_dimensional_array(void *data, int ndim, int *dimsize);

/** \} */

int get_nlines_in_file(FILE *F);
#endif
