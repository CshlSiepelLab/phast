/* $Id: misc.c,v 1.3 2004-06-09 17:10:29 acs Exp $
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

#define NCODONS 64

int int_pow(int x, int y) {
  int retval = 1, i;
  for (i = 0; i < y; i++) retval *= x;
  return retval;
}

/* fill an array with 1s or zeroes, indicating a random choice of K
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
gsl_matrix* read_subst_mat(FILE *F, char *alph) {
  gsl_matrix *retval = NULL;
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
      retval = gsl_matrix_calloc(size, size);
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
          gsl_matrix_set(retval, row_idx, col_idx, val);
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
    assert(strcmp(mode, "r") == 0);
    return stdin;
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

/* normalize a probability vector.  Assumes all values are
   nonnegative.  Returns normalization constant */
double normalize_probs(double *p, int size) {
  int i;
  double sum = 0;
  for (i = 0; i < size; i++) sum += p[i];
  for (i = 0; i < size; i++) p[i] /= sum;  
  return sum;
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
