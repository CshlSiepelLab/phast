/* $Id: misc.c,v 1.11 2005-07-17 22:20:12 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

#include <stdlib.h>
#include <misc.h>
#include <lists.h>
#include <time.h>
#include <stacks.h>
#include <stringsplus.h>
#include <assert.h>
#include <stdarg.h>
#include <hashtable.h>

#define NCODONS 64

int int_pow(int x, int y) {
  int retval = 1, i;
  for (i = 0; i < y; i++) retval *= x;
  return retval;
}

/* fill an array with 1s or zeroes, indicating a random choice of k
   elements from a list of N.  The array 'selections' must already be
   allocated to be of length N, and should be initialized.
   Initialization may be all zeroes, if all elements are "eligible"
   selections, or may include -1s designating prohibited elements.  If
   k >= N (or, more precisely, the number of eligible elements in N),
   then all (eligible) items will be marked as selected.  */
void choose(int *selections, int N, int k) {
  int i;
  List *eligible = lst_new_int(N);
  for (i = 0; i < N; i++) {
    if (selections[i] != -1) {  /* we'll consider anything that's not
                                   a -1 to be eligible, and we'll make
                                   sure that it is marked 0 */
      lst_push_int(eligible, i);
      selections[i] = 0;
    }
  }

  srand(time(NULL));
  for (i = 0; i < k && lst_size(eligible) > 0; i++) {
    int randidx = rint(1.0 * (lst_size(eligible)-1) * rand() / RAND_MAX);
    int item = lst_get_int(eligible, randidx);
    selections[item] = 1;
    
    /* replace selected element with last one in list, then shorten list */
    item = stk_pop_int(eligible);
    if (randidx != lst_size(eligible)) /* it's possible that the
                                          selected item was the last
                                          one in the list */
        lst_set_int(eligible, randidx, item);    
  }
  lst_free(eligible);
}

/* produce a random permutation of the designated size; 'permutation'
   must be allocated externally  */
void permute(int *permutation, int N) {
  int i, element, randidx;
  List *eligible = lst_new_int(N);
  for (i = 0; i < N; i++) lst_push_int(eligible, i);

  srand(time(NULL));
  for (i = 0; i < N; i++) {
    randidx = rint(1.0 * (lst_size(eligible)-1) * rand() / RAND_MAX);
    permutation[i] = lst_get_int(eligible, randidx);
    
    /* replace selected element with last one in list, then shorten list */
    element = stk_pop_int(eligible);
    if (randidx != lst_size(eligible)) /* it's possible that the
                                          selected item was the last
                                          one in the list */
      lst_set_int(eligible, randidx, element);    
  }
  lst_free(eligible);  
}

/* Given an alphabet (assumed to contain the characters 'A', 'C', 'G',
   and 'T'), return a character array of size alph_size^3 such that
   each codon (tuple of three characters) maps to the corresponding
   amino acid, according to the universal genetic code.  Stop codons
   will be mapped to '$' characters, and tuples not in {A,C,G,T}^3 to
   null characters ('\0'). */
char* get_codon_mapping(char *alphabet) {
  int alph_size = strlen(alphabet);
  int nstates = int_pow(alph_size, 3);
  char *retval = (char*)smalloc(nstates * sizeof(char)); 
                                /* NOTE: no null terminator */
  int inv_alph[256];
  int i;

  static char *codons[] = 
    { "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
      "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
      "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
      "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
      "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG",
      "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
      "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
      "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG" };
  
  static char *aas[] = 
    { "F", "F", "L", "L", "S", "S", "S", "S",
      "Y", "Y", "$", "$", "C", "C", "$", "W",
      "L", "L", "L", "L", "P", "P", "P", "P",
      "H", "H", "Q", "Q", "R", "R", "R", "R",
      "I", "I", "I", "M", "T", "T", "T", "T",
      "N", "N", "K", "K", "S", "S", "R", "R",
      "V", "V", "V", "V", "A", "A", "A", "A",
      "D", "D", "E", "E", "G", "G", "G", "G" };

  /* create inverse alphabet mapping */
  for (i = 0; i < 256; i++) inv_alph[i] = -1;
  for (i = 0; alphabet[i] != '\0'; i++) inv_alph[(int)alphabet[i]] = i;
  
  for (i = 0; i < nstates; i++) retval[i] = '\0';
  for (i = 0; i < NCODONS; i++) 
    retval[tuple_index(codons[i], inv_alph, alph_size)] = aas[i][0];

  return retval;
}

/* given an inverse alphabet (mapping every character to an index),
   return the state number corresponding to the specified tuple of
   characters.  A "digital" indexing scheme is assumed in which the
   right-most character in the tuple is considered the least
   significant digit.  Returns -1 if a character is encountered that
   is not in the specified alphabet. */
int tuple_index(char *tuple, int *inv_alph, int alph_size) {
  int retval = 0, i;
  int tuple_size = strlen(tuple);
  for (i = 0; i < tuple_size; i++) {
    if (inv_alph[(int)tuple[tuple_size-i-1]] < 0) return -1;
    retval += inv_alph[(int)tuple[tuple_size-i-1]] * int_pow(alph_size, i);
                                /* i == 0 => least sig. dig; 
                                   i == tuple_size-1 => most sig. dig */
  }
  return retval;
}

/* inverse of above: given an alphabet and a tuple idx, return the
   corresponding string.  The number of characters in the tuple is
   also required  */
void get_tuple_str(char *tuple_str, int tuple_idx, int tuple_size, 
                   char *alphabet) {
  int k, remainder, alph_size = strlen(alphabet);
  remainder = tuple_idx;
  for (k = 0; k < tuple_size; k++) {
    tuple_str[tuple_size-k-1] = alphabet[remainder % alph_size];
    remainder /= alph_size;
  }
}

/* Read a substitution matrix from the specified file (format expected
   to be that used by BLAST, as output by the NCBI "pam" program).  If
   'alph' equals the empty string (must be preallocated to adequate
   size), then the sequence of characters that defines the matrix rows
   and columns will be stored in it; otherwise, 'alph' will be taken
   to define the desired alphabet, and the order of rows and columns
   will be rearranged as necessary to be consistent with it
   (characters in the file but not in 'alph' will be ignored).*/
/* WARNING: only minor testing perfomed.  Check carefully if use
   without predefined alphabet or with alphabets that do not match in
   order */
Matrix* read_subst_mat(FILE *F, char *alph) {
  Matrix *retval = NULL;
  int i = 0, j, size = 0, file_size = 0, row_idx;
  List *fields = lst_new_ptr(100);
  String *line = str_new(STR_MED_LEN), *rowstr;
  int predefined_alph = strlen(alph) == 0 ? 0 : 1;
  int inv_alph[256];
  char file_alph[STR_SHORT_LEN];
  
  /* if alphabet is predefined, we need to be able to map characters
     to states */
  if (predefined_alph) {
    for (i = 0; i < 256; i++) inv_alph[i] = -1;
    for (i = 0; alph[i] != '\0'; i++) inv_alph[(int)alph[i]] = i;
  }

  for (i = 0; (retval == NULL || i < file_size) && 
         str_readline(line, F) != EOF; ) {
    str_double_trim(line);
    if (str_starts_with_charstr(line, "#") || line->length == 0) continue;

    /* assume first non-comment, non-empty line defines alphabet
       and matrix size */
    if (retval == NULL) {
      str_remove_all_whitespace(line);
      strcpy(file_alph, line->chars);
      if (!predefined_alph) strcpy(alph, line->chars);
      file_size = strlen(file_alph); 
      size = strlen(alph);
      retval = mat_new(size, size);
    }
    else {
      str_split(line, NULL, fields);
      if (lst_size(fields) != file_size + 1) {
        fprintf(stderr, "ERROR: unexpected number of columns for row %d.\n",
                i+1);
        exit(1);
      }
      rowstr = lst_get_ptr(fields, 0);
      if (rowstr->chars[0] != file_alph[i]) {
        fprintf(stderr, "ERROR: unexpected row label in row %d\n", i+1);
        exit(1);
      } 
      str_free(rowstr);

      row_idx = (predefined_alph ? inv_alph[(int)file_alph[i]] : i);
      i++;
      if (row_idx == -1) continue;


      for (j = 0; j < file_size; j++) {
        double val;
        int col_idx = (predefined_alph ? inv_alph[(int)file_alph[j]] : j);
        if (col_idx != -1) {
          if (str_as_dbl(lst_get_ptr(fields, j+1), &val) != 0) {
            fprintf(stderr, "ERROR: non-numeric matrix element in subst. matrix ('%s')\n", 
                    ((String*)lst_get_ptr(fields, j+1))->chars);
            exit(1);
          }
          mat_set(retval, row_idx, col_idx, val);
        }
        str_free(lst_get_ptr(fields, j+1));
      }
    }
  }

  if (i != file_size) {
    fprintf(stderr, "ERROR: too few rows in subst. matrix.\n");
    exit(1);
  }

  lst_free(fields);
  str_free(line);
  return retval;
}

/* simple wrapper for fopen that opens specified filename or aborts
   with appropriate error message.  Saves typing in mains for
   command-line programs */
FILE* fopen_fname(char *fname, char *mode) {
  FILE *F = NULL;
  if (!strcmp(fname, "-")) {
    if (strcmp(mode, "r") == 0)
      return stdin;
    else if (strcmp(mode, "w+") == 0)
      return stdout;
    else die("ERROR: bad args to fopen_fname.\n");
  }
  if ((F = fopen(fname, mode)) == NULL) {
    fprintf(stderr, "ERROR: cannot open %s.\n", fname);
    exit(1);
  }
  return F;
}

/* print error message and die with exit 1; saves typing in mains */
void die(char *warnfmt, ...) {
  va_list args;
 
  va_start(args, warnfmt);
  vfprintf(stderr, warnfmt, args);
  va_end(args);
  exit(1);
}

/* returns List of Strings derived from an argument that may either be
   a literal comma-separated list or a reference to a file (using the
   "*" convention).  Note: List and all Strings are newly allocated
   (should be freed externally) */
List *get_arg_list(char *arg) {
  String *argstr = str_new_charstr(arg);
  List *l = lst_new_ptr(10);
  if (str_starts_with_charstr(argstr, "*")) {
    FILE *F;
    String *fname_str;
    if ((F = fopen(&argstr->chars[1], "r")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open file %s.\n", &argstr->chars[1]);
      exit(1);
    }
    fname_str = str_new(STR_MED_LEN);
    str_slurp(fname_str, F);
    str_split(fname_str, NULL, l);
    fclose(F);
    str_free(fname_str);
  }
  else {
    /* if contains commas, assume comma delimited, otherwise assume
       whitespace delimited */
    char *delim = NULL; int i;
    for (i = 0; i < argstr->length && argstr->chars[i] != ','; i++);
    if (i < argstr->length) delim = ",";      
    str_split(argstr, delim, l);
  }

  str_free(argstr);
  return l;
}

/** Returns remaining command-line arguments as a List of Strings.  */
List *remaining_arg_list(char *argv[], /**< As passed to main  */
                         int argc, /**< As passed to main  */
                         int optind /**< Index of first unprocessed argument */
                ) {
  List *l = lst_new_ptr(argc - optind);
  int i;
  for (i = optind; i < argc; i++) 
    lst_push_ptr(l, str_new_charstr(argv[i]));
  return l;
}

List *get_arg_list_int(char *arg) {
  List *l = get_arg_list(arg);
  List *retval = str_list_as_int(l);
  lst_free_strings(l);
  lst_free(l);
  return retval;
}

List *get_arg_list_dbl(char *arg) {
  List *l = get_arg_list(arg);
  List *retval = str_list_as_dbl(l);
  lst_free_strings(l);
  lst_free(l);
  return retval;
}

/* argument conversion with error checking (int) */
int get_arg_int(char *arg) {
  char *endptr;
  int retval = (int)strtol(arg, &endptr, 0);
  if (*endptr != '\0') die("ERROR: cannot parse integer '%s'\n", arg);
  return retval;
}

/* argument conversion with error checking (double) */
double get_arg_dbl(char *arg) {
  char *endptr;
  double retval = (double)strtod(arg, &endptr);
  if (*endptr != '\0') die("ERROR: cannot parse floating point '%s' at command line\n", arg);
  return retval;
}

/* argument conversion with error and bounds checking (int) */
int get_arg_int_bounds(char *arg, int min, int max) {
  char *endptr;
  int retval = (int)strtol(arg, &endptr, 0);
  if (*endptr != '\0') die("ERROR: cannot parse integer '%s' at command line\n", arg);
  if (retval < min || retval > max) 
    die("ERROR: integer %d at command line outside allowable range %d-%d.\n", 
        retval, min, max);
  return retval;
}

/* argument conversion with error checking (double) */
double get_arg_dbl_bounds(char *arg, double min, double max) {
  char *endptr;
  double retval = (double)strtod(arg, &endptr);
  if (*endptr != '\0') die("ERROR: cannot parse floating point '%s' at command line\n", arg);
  if (retval < min || retval > max) 
    die("ERROR: floating point %f at command line outside allowable range %f-%f.\n", 
        retval, min, max);
  return retval;
}

/* safe malloc and realloc */
void *smalloc(size_t size) {
  void *retval = malloc(size);
  if (retval == NULL) {
    fprintf(stderr, "FATAL ERROR: out of memory.\n");
    assert(0);
  }
  return retval;
}

void *srealloc(void *ptr, size_t size) {
  void *retval = realloc(ptr, size);
  if (retval == NULL && ptr != NULL && size != 0) {
    fprintf(stderr, "FATAL ERROR: out of memory.\n");
    assert(0);
  }
  return retval;
}

/* efficiently compute log of sum of values, which themselves are
   stored as logs: that is, return log(sum_i exp(l_i)).  The largest
   of the elements of l (call it maxval) is factored out, so that
   log(sum_i(exp(l_i))) = maxval + log(1 + sum_i(exp(l_i-maxval))),
   where the new sum is taken over 2 <= i < n.  All of the quantities
   in the exp must be negative, and those smaller than some reasonable
   threshold can be ignored. [Thanks to David Haussler for showing me
   this trick].  WARNING: sorts list as side effect. */
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

/* same as above, but base e */
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

/* return 1 if a change from b1 to b2 is a transition, and 0 otherwise */
int is_transition(char b1, char b2) {
  b1 = toupper(b1); b2 = toupper(b2);
  return ((b1 == 'A' && b2 == 'G') || 
          (b1 == 'G' && b2 == 'A') || 
          (b1 == 'T' && b2 == 'C') || 
          (b1 == 'C' && b2 == 'T'));
}

/** Normalize a probability vector.  Assumes all values are
   nonnegative.  Returns normalization constant */
double normalize_probs(double *p, int size) {
  int i;
  double sum = 0;
  for (i = 0; i < size; i++) sum += p[i];
  for (i = 0; i < size; i++) p[i] /= sum;  
  return sum;
}


/** Make 'n' draws from a uniform distribution on the interval [min,
   max], optionally with antithetics.  Store in 'draws'.  Designed for use with
   real (floating-point) numbers.  Be sure to call srand externally. */ 
void unif_draw(int n, double min, double max, double *draws, int antithetics) {
  int i;
  double range = max - min;
  for (i = 0; i < n; i++) {
    draws[i] = min + range * rand()/RAND_MAX;
    if (antithetics) {
      draws[i+1] = min + (max - draws[i]);
      i++;
    }
  }
}

/** Make 'n' draws from a binomial distribution with parameters 'N'
   and 'p'.  Store numbers of successes in 'draws'.  Be sure to call
   srand externally.  WARNING: computational complexity is O(n*N) --
   see Numerical Recipes for a better way (rejection sampling) */
void bn_draw(int n, int N, double p, int *draws) {
  int i, j;
  double *unif_draws;
  assert(n >= 1 && N >= 1);
  unif_draws = smalloc(N * sizeof(double));
  for (i = 0; i < n; i++) {
    unif_draw(N, 0, 1, unif_draws, FALSE);
                                /* antithetics can have undesirable
                                   effect; e.g., if p = 0.5, it will
                                   always be true that draws[i] =
                                   N/2 */
    draws[i] = 0;
    for (j = 0; j < N; j++) if (unif_draws[j] < p) draws[i]++;
                                /* number of uniform draws less than p
                                   is binomial */
  }
  free(unif_draws);
}

/** Make 'n' draws from a multinomial distribution defined by
   probability vector 'p' with dimension 'd'.  Record the counts for
   each category in 'counts'.  Sum of elements in 'counts' will equal
   'n'.  Probability vector is assumed to be normalized.  WARNING:
   computational complexity is O(n * d) -- need a better
   implementation of bn_draw to handle large n efficiently.  Be sure
   to call srand externally. */
void mn_draw(int n, double *p, int d, int *counts) {
  int i, nremaining = n;
  double cum_p = 0;

  /* recursively make draws from binomial */
  for (i = 0; i < d-1; i++) {
    if (p[i] == 0 || nremaining == 0) {
      counts[i] = 0; /* save time and avoid possible div by 0 */
      continue;
    }
    bn_draw(1, nremaining, p[i] / (1-cum_p), &counts[i]);
    nremaining -= counts[i];
    cum_p += p[i];
  }
  /* last category is constrained by prev */
  counts[d-1] = nremaining;
}

/** Parse string defining mapping from old names to new and store as
   hash (e.g., "hg17->human; mm5->mouse; rn3->rat").  Can use "="
   character as well as "->" to indicate mapping.  */
struct hash_table *make_name_hash(char *mapstr) {
  Hashtable *retval = hsh_new(20);
  Regex *map_re = str_re_new("^[[:space:]]*([A-Za-z0-9_]+)[[:space:]]*(->|=)[[:space:]]*([A-Za-z0-9_]+)[[:space:]]*");
  List *mappings = lst_new_ptr(20), *names = lst_new_ptr(3);
  String *s = str_new_charstr(mapstr);
  int i;

  str_split(s, ";", mappings);
  for (i = 0; i < lst_size(mappings); i++) {
    String *oldname, *newname;
    if (str_re_match(lst_get_ptr(mappings, i), map_re, names, 3) < 0)
      die("ERROR: cannot parse mapping ('%s')\n", lst_get_ptr(mappings, i));
    oldname = lst_get_ptr(names, 1);
    newname = lst_get_ptr(names, 3);
    hsh_put(retval, oldname->chars, strdup(newname->chars));
    lst_free_strings(names);
  }
  lst_free_strings(mappings);
  lst_free(mappings);
  lst_free(names);
  str_free(s);
  str_re_free(map_re);

  return retval;
}

/** Evaluate pdf of gamma distribution with parameters a and b */
double gamma_pdf(double x, double a, double b) {
  return 1/(gamma(a) * pow(b, a)) * pow(x, a-1) * exp(-x/b);
}

/* make a draw from an exponential distribution with parameter
   (expected value) 'b' */
double exp_draw(double b) {
  return -log(1.0 * rand() / RAND_MAX) * b;
}

/* make a draw from a gamma distribution with parameters 'a' and
   'b'. Be sure to call srand externally.  If a == 1, exp_draw is
   called.  If a > 1, Best's (1978) rejection algorithm is used, and
   if a < 1, rejection sampling from the Weibull distribution is
   performed, both as described in "Non-Uniform Random Variate
   Generation" by Luc Devroye, available online at
   http://cgm.cs.mcgill.ca/~luc/rnbookindex.html */
double gamma_draw(double a, double b) {
  double retval = -1;

  assert(a > 0);

  if (a == 1) return exp_draw(b);

  else if (a > 1) {
    while (retval == -1) {
      double U, V, W, X, Y, Z;
      double d = a - 1, c = 3 * a - 0.75;

      U = 1.0*rand()/RAND_MAX;	/* uniform on [0, 1]; used for draw */
      V = 1.0*rand()/RAND_MAX;	/* also uniform on [0, 1]; used for
				   rejection/acceptance */
      W = U * (1 - U);
      Y = sqrt(c / W) * (U - 0.5); 
      X = d + Y;
      /* Y is a scaled version of a random variate from a t distribution
	 with 2 dof; X is a random deviate from a shifted distribution
	 whose ratio with a gamma(a, 1) can be bounded by a constant */

      if (X < 0) continue;	/* truncate because of nonnegativity
				   of gamma */
    
      Z = 64 * W * W * W * V * V;
      if (log(Z) <= 2 * (d * log(X / d) - Y))
	retval = X;
      /* there's some algebra behind this, but underneath it's just
	 the standard rejection sampling idea of accepting with
	 probability p(x)/(M*q(x)), where p(x) is the density of the
	 desired distrib, q(x) is the density of the distrib from which
	 you have sampled, and M is a constant such that p(x) / q(x) <=
	 M for all x */
    }
  }

  else {			/* use Weibull */
    double c = 1/a, d = pow(a, a/(1-a)) * (1-a);
    while (retval == -1) {
      double E, Z, X;
      E = exp_draw(1);
      Z = exp_draw(1);
      X = pow(Z, c);		/* X is Weibull(a) */
      if (Z + E >= d + X)	/* note: wrong in book, correct
				   formula in errata */
	retval = X;
    }
  }

  /* so far we only have a draw from gamma(a, 1); multiply by b to
     obtain a draw from gamma(a, b) */
  return retval * b;
}

/* make a draw from a k-dimensional Dirichlet distribution, with
   parameters given by alpha[0], ..., alpha[k-1].  This is
   accomplished by sampling from gamma distributions with parameters
   alpha[0], ..., alpha[k-1] and renormalizing */
void dirichlet_draw(int k, double *alpha, double *theta) {
  int i;
  for (i = 0; i < k; i++) 
    theta[i] = gamma_draw(alpha[i], 1);
  normalize_probs(theta, k);
}

/* return n! */
int permutations(int n) {
  int i, retval = 1;
  for (i = 2; i <= n; i++) retval *= i;
  return retval;
}

/* return n-choose-k */
int combinations(int n, int k) {
  int i, retval = 1;
  for (i = 0; i < k; i++) retval *= (n - i);
  return retval / permutations(k);
}

/* compute relative entropy in bits of q with respect to p, both
   probability vectors of dimension d */
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
double sym_rel_entropy(double *p, double *q, int d) {
  double re1 = rel_entropy(p, q, d), re2 = rel_entropy(q, p, d);
  return min(re1, re2);
}

/***************************************************************************/
/* for debugging: these functions can be called dynamically in gdb to
   print the contents of 1d and 2d arrays */
/***************************************************************************/

void sq_matrix_pretty_print(double **mat, int D) {
  int i, j;
  for (i = 0; i < D; i++) {
    for (j = 0; j < D; j++) 
      printf("%8.6f ", mat[i][j]);
    printf("\n");
  }
}

void int_vector_print(int *vect, int D) {
  int i;
  for (i = 0; i < D; i++) printf("%d: %d\n", i, vect[i]);
}

void dbl_vector_print(double *vect, int D) {
  int i;
  for (i = 0; i < D; i++) printf("%d: %f\n", i, vect[i]);
}
