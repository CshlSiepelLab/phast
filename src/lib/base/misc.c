/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: misc.c,v 1.37 2008-12-10 18:09:17 agd27 Exp $ */

#include <stdlib.h>
#include <misc.h>
#include <lists.h>
#include <stacks.h>
#include <stringsplus.h>
#include <stdarg.h>
#include <hashtable.h>
#include <unistd.h>
#include <assert.h>

#define NCODONS 64

//avoid conflict with R
#undef choose

#ifdef RPHAST
FILE *rphast_stdout=(FILE*)0;
FILE *rphast_stderr=(FILE*)1;
#endif

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

  /* if RPHAST RNG is seeded externally */
#ifndef RPHAST
  srandom((unsigned int)time(NULL));
#endif
  for (i = 0; i < k && lst_size(eligible) > 0; i++) {
    int randidx = (int)rint(1.0 * (lst_size(eligible)-1) * unif_rand());
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
#ifndef RPHAST
#if defined(__MINGW32__)
/*these functions inplemented by lib iberty in mingw do not work correctly
the C standard rand() and srand() have been subsituted and give an answer close
but not exactly the same as the linux version.*/
int random()
{
  return rand();
}

void srandom(int seed)
{
  srand(seed);
}
#endif
#endif
/* produce a random permutation of the designated size; 'permutation'
   must be allocated externally  */
void permute(int *permutation, int N) {
  int i, element, randidx;
  List *eligible = lst_new_int(N);
  for (i = 0; i < N; i++) lst_push_int(eligible, i);

  /* if RPHAST RNG is seeded externally */
#ifndef RPHAST
  srandom((unsigned int)time(NULL));
#endif

  for (i = 0; i < N; i++) {
    randidx = (int)rint(1.0 * (lst_size(eligible)-1) * unif_rand());
    if (!(randidx >= 0 && randidx < lst_size(eligible)))
      die("ERROR permute: randidx=%i, should be in [0, %i)\n",
	  randidx, lst_size(eligible));
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
  int alph_size = (int)strlen(alphabet);
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
  int tuple_size = (int)strlen(tuple);
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
  int k, remainder, alph_size = (int)strlen(alphabet);
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
      file_size = (int)strlen(file_alph); 
      size = (int)strlen(alph);
      retval = mat_new(size, size);
      mat_zero(retval);
    }
    else {
      str_split(line, NULL, fields);
      if (lst_size(fields) != file_size + 1) {
        die("ERROR: unexpected number of columns for row %d.\n",
	    i+1);
      }
      rowstr = lst_get_ptr(fields, 0);
      if (rowstr->chars[0] != file_alph[i]) {
        die("ERROR: unexpected row label in row %d\n", i+1);
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
            die("ERROR: non-numeric matrix element in subst. matrix ('%s')\n", 
		((String*)lst_get_ptr(fields, j+1))->chars);
          }
          mat_set(retval, row_idx, col_idx, val);
        }
        str_free(lst_get_ptr(fields, j+1));
      }
    }
  }

  if (i != file_size) {
    die("ERROR: too few rows in subst. matrix.\n");
  }

  lst_free(fields);
  str_free(line);
  return retval;
}

/* simple wrapper for fopen that opens specified filename or aborts
   with appropriate error message.  Saves typing in mains for
   command-line programs */
FILE* phast_fopen_no_exit(const char *fname, const char *mode) {
  FILE *F = NULL;
  if (!strcmp(fname, "-")) {
    if (mode[0]=='r') 
      return stdin;
    else if (mode[0]=='w')
      return stdout;
    else die("ERROR: bad args to phast_fopen.\n");
  }
  F = fopen(fname, mode);
  if (F != NULL) register_open_file(F);
  return F;
}

FILE* phast_fopen(const char *fname, const char *mode) {
  FILE *F = phast_fopen_no_exit(fname, mode);
  if (F == NULL)
    die("ERROR: cannot open %s.\n", fname);
  return F;
}

void phast_fclose(FILE *f) {
  if (f != stdout && f!=stderr) {
    fclose(f);
    unregister_open_file(f);
  }
}

/* print error message and die with exit 1; saves typing in mains */
#ifndef RPHAST
void die(const char *warnfmt, ...) {
  va_list args;

  va_start(args, warnfmt);
  vfprintf(stderr, warnfmt, args);
  va_end(args);
#ifdef PHAST_DEBUG
  assert(0);
#endif
  exit(1);
}

void phast_warning(const char *warnfmt, ...) {
  va_list args;
  va_start(args, warnfmt);
  vfprintf(stderr, warnfmt, args);
  va_end(args);
}

void set_seed(int seed) {
  if (seed < 0) {
    struct timeval now;
    gettimeofday(&now, NULL);
    srandom((unsigned int)now.tv_usec);
  } else srandom(seed);
}

#endif

#ifdef RPHAST

//seed is ignored in RPHAST mode!
void set_seed(int seed) {
  GetRNGstate();
}


int rphast_fprintf(FILE *f, const char *format, ...) {
  va_list args;
  va_start(args, format);
  if (f == stdout)
    Rvprintf(format, args);
  else if (f == stderr)
    REvprintf(format, args);
  else return vfprintf(f, format, args);
  return 1;
}

#endif

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
    F = phast_fopen(&argstr->chars[1], "r");
    fname_str = str_new(STR_MED_LEN);
    str_slurp(fname_str, F);
    str_split(fname_str, NULL, l);
    phast_fclose(F);
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
List *remaining_arg_list(char *argv[], int argc, int optind ) {
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

//rphast versions defined in rph_util.c
#ifndef USE_PHAST_MEMORY_HANDLER
/* safe malloc and realloc */
void *smalloc(size_t size) {
  void *retval = malloc(size);
  if (retval == NULL)
    die("FATAL ERROR: out of memory.\n");
  return retval;
}

void *srealloc(void *ptr, size_t size) {
  void *retval = realloc(ptr, size);
  if (retval == NULL && ptr != NULL && size != 0)
    die("FATAL ERROR: out of memory.\n");
  return retval;
}
#endif

/* make a copy of word, allocating just enough space.*/
char *copy_charstr(const char *word) {
  int len = (int)strlen(word);
  char *retval = smalloc((len+1)*sizeof(char));
  strcpy(retval, word);
  return retval;
}


/* return 1 if a change from b1 to b2 is a transition, and 0 otherwise */
int is_transition(char b1, char b2) {
  b1 = (char)toupper(b1); b2 = (char)toupper(b2);
  return ((b1 == 'A' && b2 == 'G') || 
          (b1 == 'G' && b2 == 'A') || 
          (b1 == 'T' && b2 == 'C') || 
          (b1 == 'C' && b2 == 'T'));
}

/* return 1 if a change from b1 to b2 involves an insertion/deletion,
   0 otherwise */
int is_indel(char b1, char b2) {
  return (b1 == '-' || b2 == '-');
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
   real (floating-point) numbers.  Be sure to call srandom externally. */ 
void unif_draw(int n, double min, double max, double *draws, int antithetics) {
  int i;
  double range = max - min;
  for (i = 0; i < n; i++) {
    draws[i] = min + range * unif_rand();
    if (antithetics) {
      draws[i+1] = min + (max - draws[i]);
      i++;
    }
  }
}

/** Make a draw from a binomial distribution with parameters 'N' and
   'p'.  Be sure to call srandom externally.  WARNING: computational
   complexity is O(N) -- see Numerical Recipes for a better way
   (rejection sampling) */
int bn_draw(int N, double p) {
  int j, retval = 0;
  double *unif_draws;
  if (N < 1)
    die("ERROR bn_draw: got N=%i\n", N);
  unif_draws = smalloc(N * sizeof(double));
  unif_draw(N, 0, 1, unif_draws, FALSE);
                                /* antithetics can have undesirable
                                   effect; e.g., if p = 0.5, it will
                                   always be true that retval =
                                   N/2 */
  for (j = 0; j < N; j++) if (unif_draws[j] < p) retval++;
                                /* number of uniform draws less than p
                                   is binomial */
  sfree(unif_draws);
  return retval;
}

/* Make a draw from a binomial distribution with parameters 'n' and
   'pp'.  This version uses rejection sampling unless n < 25 or n * p <
   1/25.  It takes constant expected time.  Be sure to call srandom
   externally. */
int bn_draw_fast(int n, double pp) {
  int j;
  static int nold = -1;
  double am, em, g, angle, p, bn1, sq, t, y;
  static double pold = -1, pc, plog, pclog, en, oldg;

  if (n < 25) return bn_draw(n, pp);

  p = (pp <= 0.5 ? pp : 1.0 - pp); /* can assume p less than 0.5, and
                                      adjust return value as
                                      necessary */

  am = n * p;

  if (am < 1.0) {               /* if fewer than one event out of 25
                                   is expected, distribution is quite
                                   accurately Poisson; use direct
                                   Poisson method */
    g = exp(-am);
    t = 1.0;
    for (j = 0; j <= n; j++) {
      t *= unif_rand();
      if (t < g) break;
    }
    bn1 = min(j, n);
  }
  else {                        /* use rejection method */
    if (n != nold) {            /* compute first time only */
      en = n;
      oldg = lgamma(en + 1);
      nold = n;
    }
    if (p != pold) {            /* compute first time only */
      pc = 1 - p;
      plog = log(p);
      pclog = log(pc);
      pold = p;
    }
    sq = sqrt(2 * am * pc);     /* rejection method with Lorentzian
                                   comparison function */
    do {
      do {
        angle = M_PI * unif_rand();
        y = tan(angle);
        em = sq * y + am;
      } 
      while(em < 0 || em >= en + 1); /* reject */
      em = floor(em);
      t = 1.2*sq*(1+y*y)* 
        exp(oldg - lgamma(em+1) - lgamma(en-em+1) + 
            em*plog + (en-em)*pclog);
    } 
    while (unif_rand() > t);
    bn1 = em;
  }

  if (p != pp) bn1 = n - bn1;   /* undo symmetry transformation */
  return (int)bn1;
}

/* structure and comparison function used in mn_draw */
struct mndata {
  int idx;
  double p;
  double count;
};

int mn_compare(const void* ptr1, const void* ptr2) {
  const struct mndata *d1 = ptr1, *d2 = ptr2;
  if (d1->p == d2->p) return 0;
  else if (d1->p < d2->p) return 1;
  else return -1;
}

/** Make 'n' draws from a multinomial distribution defined by
   probability vector 'p' with dimension 'd'.  Record the counts for
   each category in 'counts'.  Sum of elements in 'counts' will equal
   'n'.  Probability vector is assumed to be normalized.  Be sure
   to call srandom externally. */
void mn_draw(int n, double *p, int d, int *counts) {
  int i, nremaining = n;
  double rem_p = 1;

  /* sort in descending order by probability */
  struct mndata *data = smalloc(d * sizeof(struct mndata));
  for (i = 0; i < d; i++) {
    data[i].idx = i;
    data[i].p = p[i];
    data[i].count = 0;
  }
  qsort(data, d, sizeof(struct mndata), mn_compare);

  /* recursively make draws from binomial */
  for (i = 0; i < d-1; i++) {
    if (data[i].p == 0 || nremaining == 0) 
      break;
    data[i].count = bn_draw_fast(nremaining, data[i].p / rem_p);
    nremaining -= data[i].count;
    rem_p -= data[i].p;
  }
  /* last category is constrained by prev */
  data[d-1].count = nremaining;

  /* now populate counts */
  for (i = 0; i < d; i++)
    counts[data[i].idx] = (int)data[i].count;

  sfree(data);
}

/** Given a probability vector, draw an index.  Call srandom externally */
int draw_index(double *p, int size) {
  int i;
  double sum = 0, r = unif_rand();
  
  for (i = 0; i < size; i++) {
    sum += p[i];
    if (r < sum) break;
  }
  if (i == size) i = size-1;	/* to be safe */
  return i;
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
    hsh_put(retval, oldname->chars, copy_charstr(newname->chars));
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
  double tgamma(double d);
  if (x < 0)
    die("ERROR gamma_pdf got x=%f\n", x);
  return 1/(tgamma(a) * pow(b, a)) * pow(x, a-1) * exp(-x/b);
}

/** Evaluate cdf of gamma distribution with parameters a and b.  If
    lower_tail == TRUE returns P(X<=x), else returns P(X>=x) */
double gamma_cdf(double x, double a, double b, int lower_tail) {
  if (x < 0)
    die("ERROR gamma_cdf got x=%f\n", x);
  return incomplete_gamma(a, x/b, lower_tail ? 'p' : 'q');
}

/* Evaluate pdf of chi-square distribution with dof degrees of freedom */
double chisq_pdf(double x, double dof) {
  if (x < 0)
    die("ERROR chisq_pdf got x=%f\n", x);
  return gamma_pdf(x, dof/2, 2);
}

/* Evaluate cdf of chi-square distribution with dof degrees of
    freedom.  If lower_tail == TRUE returns P(X<=x), else returns
    P(X>=x) */
double chisq_cdf(double x, double dof, int lower_tail) {
  if (x < 0)
    die("ERROR chisq_cdf got x=%f\n", x);
  return gamma_cdf(x, dof/2, 2, lower_tail);
}

/* Evaluate cdf of a 50:50 mixture of a chisq distribution with dof
   degrees of freedom and a point mass at 0.  If lower_tail == TRUE
   returns P(X<=x), else returns P(X>=x).  This function is useful in
   likelihood ratio tests of bounded parameters. */
double half_chisq_cdf(double x, double dof, int lower_tail) {
  double retval = 0.5 * chisq_cdf(x, dof, lower_tail);
  if (lower_tail || x == 0) retval += 0.5;
  return(retval);
}

/* make a draw from an exponential distribution with parameter
   (expected value) 'b' */
double exp_draw(double b) {
  return -log(unif_rand()) * b; /* inversion method */
}

/* make a draw from a gamma distribution with parameters 'a' and
   'b'. Be sure to call srandom externally.  If a == 1, exp_draw is
   called.  If a > 1, Best's (1978) rejection algorithm is used, and
   if a < 1, rejection sampling from the Weibull distribution is
   performed, both as described in "Non-Uniform Random Variate
   Generation" by Luc Devroye, available online at
   http://cgm.cs.mcgill.ca/~luc/rnbookindex.html */
double gamma_draw(double a, double b) {
  double retval = -1;

  if (a <= 0)
    die("ERROR gamma_draw got a=%f\n", a);

  if (a == 1) return exp_draw(b);

  else if (a > 1) {
    while (retval == -1) {
      double U, V, W, X, Y, Z;
      double d = a - 1, c = 3 * a - 0.75;

      U = unif_rand();  	/* uniform on [0, 1]; used for draw */
      V = unif_rand();	        /* also uniform on [0, 1]; used for
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

/* evaluate density of Beta distribution */
double d_beta(double x, double a, double b) {
  double lb;
  if (!(x >= 0 && x <= 1 && a >= 0 && b >= 0))
    die("ERROR d_beta got x=%f, a=%f, b=%f\n", x, a, b);
  lb = lgamma(a+b) - lgamma(a) - lgamma(b) + (a-1) * log(x) + (b-1) * log(1-x);
  return (exp(lb));
}

/* make a draw from a beta distribution with parameters 'a' and 'b'.  */
double beta_draw(double a, double b) {
  double x = gamma_draw(a, 1);
  double y = gamma_draw(b, 1);
  return x / (x + y);
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

/* incomplete gamma function */
double incomplete_gamma(double a, 
			double x, 
			char type /* 'p' means first half of integral,
				     'q' means second half; see p. 171 */
                        ) {
  double gln;
  int n;
  double retval = -1;

  if (!(x >= 0 && a > 0 && (type == 'p' || type == 'q')))
    die("ERROR incomplete_gamma got x=%f, a=%f, type=%c\n", x, a, type);

  gln = lgamma(a);

  if (x < a + 1) {		/* use series representation */
    double ap, del, sum;

    ap = a;
    del = sum = 1/a;
    for (n = 1; n <= 200; n++) {
      ap += 1.0;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum) * 3.0e-7) {
        retval = sum * exp(-x + a * log(x) - gln);
        break;
      }
    }
  }
  else {			/* use continued fraction representation */
    double gold = 0.0, g=-1.0, fac = 1.0, b1 = 1.0;
    double b0 = 0, anf, ana, an, a1, a0 = 1.0;

    a1 = x;
    for (n = 1; n <= 200; n++) {
      an = n;
      ana = an - a;
      a0 = (a1 + a0 * ana) * fac;
      b0 = (b1 + b0 * ana) * fac;
      anf = an * fac;
      a1 = x * a0 + anf * a1;
      b1 = x * b0 + anf * b1;
      if (a1) {
        fac = 1/a1;
        g = b1 * fac;
        if (fabs((g - gold) / g) < 3.0e-7) {
          retval = 1 - exp(-x + a * log(x) - gln) * g;
          break;
        }
      }
      gold = g;
    }
  }

  if (retval == -1)		/* failed to converge */
    fprintf(stderr, "WARNING: 'a' too large in incomplete_gamma.\n");

  else if (type == 'q') 
    retval = 1 - retval;

  return retval;
}

/* return P(x = k | lambda), for a variable x that obeys a Poisson
   distribution with parameter lambda */
double d_poisson(double lambda, int k) {
  if (!(lambda >= 0 && k >= 0))
    die("ERROR d_poisson got lambda=%f, k=%i\n", lambda, k);
  return exp(-lambda + k * log(lambda) - lgamma(k+1));
}

/* return P(x <= k | lambda), for a variable x that obeys a Poisson
   distribution with parameter lambda */
double cum_poisson(double lambda, int k) {
  if (!(lambda >= 0 && k >= 0))
    die("ERROR cum_poisson got lambda=%f, k=%i\n", lambda, k);
  return incomplete_gamma(k+1, lambda, 'q');
}

/* return P(x > k | lambda), for a variable x that obeys a Poisson
   distribution with parameter lambda */
double cum_poisson_c(double lambda, int k) {
  if (!(lambda >= 0 && k >= 0))
    die("ERROR cum_poisson_c got lambda=%f, k=%i\n", lambda, k);
  return incomplete_gamma(k+1, lambda, 'p');
}

/* return P(x <= a | mu, sigma) for a variable a that obeys a normal
   distribution with mean mu and s.d. sigma */
double cum_norm(double mu, double sigma, double a) {
  if (mu != 0 || sigma != 1)
    a = (a - mu) / sigma;
  if (a >= 0) 
    return 0.5 * (1 + erf(a / sqrt(2)));
  else 
    return cum_norm_c(0, 1, -a);
}

/* return P(x >= a | mu, sigma) for a variable a that obeys a normal
   distribution with mean mu and s.d. sigma.  Use this function
   instead of 1-cum_norm when a is large (better precision) */
double cum_norm_c(double mu, double sigma, double a) {
  if (mu != 0 || sigma != 1)
    a = (a - mu) / sigma;
  if (a >= 0) 
    return 0.5 * erfc(a / sqrt(2));
  else 
    return cum_norm(0, 1, -a);
}

/* return inverse of standard normal, i.e., inv_cum_norm(p) = a such
   that cum_norm(0, 1, a) = p.  The function is approximated using an
   algorithm by Peter Acklam given at
   http://home.online.no/~pjacklam/notes/invnorm/.  */
double inv_cum_norm(double p) {
  double p_low, p_high, q, r, x, e, u;
  
  static double a[] = {
    0,
    -3.969683028665376e+01,
    2.209460984245205e+02,
    -2.759285104469687e+02,
    1.383577518672690e+02,
    -3.066479806614716e+01,
    2.506628277459239e+00};
  
  static double b[] = {
    0,
    -5.447609879822406e+01,
    1.615858368580409e+02,
    -1.556989798598866e+02,
    6.680131188771972e+01,
    -1.328068155288572e+01};
  
  static double c[] = {
    0,
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
    4.374664141464968e+00,
    2.938163982698783e+00};
  
  static double d[] = {
    0,
    7.784695709041462e-03,
    3.224671290700398e-01,
    2.445134137142996e+00,
    3.754408661907416e+00};
  
  if (!(p > 0 && p < 1))
    die("ERROR inv_cum_norm got p=%f\n", p);
  
  p_low = 0.02425;
  p_high = 1 - p_low;

  /* rational approximation for lower region */
  if (p < p_low) {
    q = sqrt(-2*log(p));
    x = (((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6]) /
      ((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1);
  }

  /* rational approximation for central region */
  else if (p >= p_low && p <= p_high) {
    q = p - 0.5;
    r = q*q;
    x = (((((a[1]*r+a[2])*r+a[3])*r+a[4])*r+a[5])*r+a[6])*q /
      (((((b[1]*r+b[2])*r+b[3])*r+b[4])*r+b[5])*r+1);
  }

  /* rational approximation for upper region */
  else {			/* p_high < p < 1 */
    q = sqrt(-2*log(1-p));
    x = -(((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6]) /
      ((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1);
  }

  /* refinement to full machine precision using Halley's rational
     method */
  e = cum_norm(0, 1, x) - p;
  u = e * sqrt(2*M_PI) * exp(x*x/2);
  x = x - u/(1 + x*u/2);

  return x;
}

/* compute min and max of (central) confidence interval of specified
   size (between 0 and 1) assuming a normal distribution with mean mu
   and s.d. sigma */
void norm_confidence_interval(double mu, double sigma, double interval_size, 
                              double *min_x, double *max_x) {
  double a;
  if (!(interval_size > 0 && interval_size < 1))
    die("ERROR norm_confidence_interval got interval_size=%f\n", interval_size);
  a = inv_cum_norm((1 - interval_size) / 2) * sigma; /* a will be negative */
  *min_x = mu + a;
  *max_x = mu - a;
}

/* density function for bivariate normal, p(x, y | mu_x, mu_y,
   sigma_x, sigma_y, rho), where mu_x and mu_y are the marginal means,
   sigma_x and sigma_y are the marginal standard deviations, and rho
   is the correlation coefficient */
double bvn_p(double x, double y, double mu_x, double mu_y, double sigma_x,
             double sigma_y, double rho) {
  double rho2 = rho*rho;
  x = (x - mu_x) / sigma_x;
  y = (y - mu_y) / sigma_y;
  return (1/(2*M_PI * sigma_x * sigma_y * sqrt(1-rho2)) *
          exp(-0.5 * 1/(1-rho2) * (x*x - 2*rho*x*y + y*y)));
}


/* Call repeatedly to enumerate combinations.  On successful exit, the
   array index (preallocate to size k) will contain indices in [0,
   n-1] representing the next element in the series of all possible
   combinations of n elements.  On the first call, set indices[0] = -1
   and the array will be initialized appropriately.  Returns TRUE on success,
   FALSE when no more choices possible */
int next_comb(int n, int k, int *index) {
  int i;

  if (!(n > 0 && k > 0 && k <= n))
    die("ERROR next_comb got n=%i k=%i\n", n, k);

  if (index[0] == -1) {
    for (i = 0; i < k; i++) index[i] = i;
    return TRUE;
  }

  /* basic idea is to advance least significant "digit" that can
     safely be advanced, ensuring that index[i] < index[i+1] for all i
     in [0, k-1] */

  /* scan backwards for first "digit" that can safely be advanced */
  for (i = k - 1; i >= 0; i--)
    if ((i == k - 1 && index[i] < n - 1) || 
	(i < k - 1 && index[i] < index[i+1] - 1))
      break;
  
  if (i < 0) return FALSE;	/* have enumerated all possibilities */

  /* advance digit, then "sweep" forward, resetting subsequent digits
     to lowest allowable value */
  index[i]++;
  for (i++; i < k; i++)
    index[i] = index[i-1] + 1;

  return TRUE;
}

/* print a single sequence in FASTA format */
void print_seq_fasta(FILE *F, char *seq, char *name, int len) {
  int j, k;
  fprintf(F, "> %s\n", name);
  for (j = 0; j < len; j += 70) {
    for (k = 0; k < 70 && j + k < len; k++) 
      fprintf(F, "%c", seq[j+k]);
    fprintf(F, "\n");
  }
}

/* return elapsed time in seconds since start_time */
double get_elapsed_time(struct timeval *start_time) {
  struct timeval now;
  gettimeofday(&now, NULL);
  return now.tv_sec - start_time->tv_sec + 
    (now.tv_usec - start_time->tv_usec)/1.0e6;
}

/* check to see if a file is present and readable on the filesystem */
int file_exists(char *filename) {
 return (access(filename, F_OK) == 0);
}

/* build static mapping from IUPAC ambiguity characters to the bases
   that they represent */
static char **build_iupac_map() {
  char **retval = smalloc(256 * sizeof(char*));
  int i;

  for (i = 0; i < 256; i++) 
    retval[i] = NULL;

  retval['R'] = "AG";
  retval['Y'] = "CT";
  retval['S'] = "CG";
  retval['W'] = "AT";
  retval['K'] = "GT";
  retval['M'] = "AC";
  retval['D'] = "AGT";
  retval['H'] = "ACT";
  retval['B'] = "CGT";
  retval['V'] = "ACG";  

  return(retval);
}

/* accessor for static mapping */
char **get_iupac_map() {
  static char **iupac_map = NULL;
  if (iupac_map == NULL) {
    iupac_map = build_iupac_map();
    set_static_var((void**)(&iupac_map));
  }    
  return iupac_map;
}

/* build an inverse mapping that allows an IUPAC ambiguity character
   to be mapped to an array, alph_size elements, with 1s or 0s
   indicating presence or absence of each character in the set
   associated with the IUPAC character.  This array obeys the indexing
   of the provided inv_states.  For use in computing likelihoods with
   ambiguity characters */
int **build_iupac_inv_map(int *inv_states, int alph_size) {
  int i, j, k;
  int **retval = smalloc(256 * sizeof(int*));
  char **iupac_map = get_iupac_map();

  memset(retval, 0, 256 * sizeof(int*));

  for (i = 0; i < 256; i++) {
    if (iupac_map[i] != NULL) {
      retval[i] = smalloc(alph_size * sizeof(int));
      for (j = 0; j < alph_size; j++) retval[i][j] = 0;
      for (j = 0; j < strlen(iupac_map[i]); j++) {
        if (inv_states[(int)iupac_map[i][j]] < 0) {
          /* in this case, we have an alphabet inconsistent with
             IUPAC; erase the previous settings, return an array of
             NULL pointers, which will cause IUPAC characters to be
             treated like 'N's in phylogenetic analysis */
          for (k = 0; k <= i; k++) {
            if (retval[k] != NULL) {
              sfree(retval[k]);
              retval[k] = NULL;
            }
          }          
          return retval;
        }
        retval[i][inv_states[(int)iupac_map[i][j]]] = 1;
      }
    }
  }
  return retval;
}

void free_iupac_inv_map(int **iim) {
  int i;
  for (i = 0; i < 256; i++) {
    if (iim[i] != NULL)
      sfree(iim[i]);
  }
  sfree(iim);
}


void *alloc_n_dimensional_array(int ndim, int *dimsize, size_t size) {
  int i;
  void **rv;
  if (ndim == 1)
    return smalloc(dimsize[0]*size);
  rv = smalloc(dimsize[0]*sizeof(void*));
  for (i=0; i < dimsize[0]; i++) 
    rv[i] = alloc_n_dimensional_array(ndim-1, &(dimsize[1]), size);
  return (void*)rv;
}

void free_n_dimensional_array(void *data, int ndim, int *dimsize) {
  int i;
  if (ndim > 1) {
    for (i=0; i < dimsize[0]; i++) {
      free_n_dimensional_array(((void**)data)[i], ndim-1, &(dimsize[1]));
    }
  }
  sfree(data);
}


int get_nlines_in_file(FILE *F) {
  char buffer[BUFFERSIZE];
  String *line = str_new(STR_MED_LEN);
  int abort = 0,  lines=0;
 
  do {
    buffer[BUFFERSIZE - 2] = '\n'; 
    if (fgets(buffer, BUFFERSIZE, F) == NULL)
      abort = 1;
    else {
      if (buffer[BUFFERSIZE - 2] == '\n' || buffer[BUFFERSIZE - 2] == '\0') {
       str_append_charstr(line, buffer);
       str_double_trim(line);
       lines = lines + (line->length > 0) ;
      } 
      str_append_charstr(line, buffer);
    }
    //Determine how many characters of the buffer were used
   
  } while (!abort);
  fseek ( F , 1 , SEEK_SET );
  return lines;

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
