/* String-handling functions, with automatic memory management and
   basic regex support.
   
   $Id: stringsplus.c,v 1.1.1.1 2004-06-03 22:43:11 acs Exp $
   Written by Adam Siepel, Summer 2002
   Copyright 2002, Adam Siepel, University of California 
*/

#include "stringsplus.h"
#include "misc.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <regex.h>

extern reg_syntax_t re_syntax_options;
extern unsigned regs_allocated;


String *str_new(int starting_nchars) {
  String *s = (String*)smalloc(sizeof(String));
  s->chars = (char*)smalloc(starting_nchars * sizeof(char) + 1);
  str_clear(s);
  s->nchars = starting_nchars;
  return s;
}

String *str_new_charstr(char *str) {
  String *s = str_new(strlen(str));
  str_cpy_charstr(s, str);
  return s;
}

void str_free(String *s) {
  free(s->chars);
  free(s);
}

void str_clear(String *s) {
  s->length = 0;
  s->chars[0] = '\0';
}

void str_concat(String *dest, String *src1, String *src2) {
  str_clear(dest);
  str_append(dest, src1);
  str_append(dest, src2);
}

void str_append(String *s, String *suffix) {
  str_nappend_charstr(s, suffix->chars, suffix->length);
}

/* stops at null terminator if occurs before len */
void str_nappend_charstr(String *s, char *charstr, int len) {
  int i;
  if (s->length + len >= s->nchars)
    str_realloc(s, max(s->length + len, s->nchars * 2) + 1);
                                /* try to avoid heavy srealloc when
                                   frequently appending short
                                   strings  */

  for (i = 0; i < len && charstr[i] != '\0'; i++)
    s->chars[s->length + i] = charstr[i];

  s->length += i;
  s->chars[s->length] = '\0';  
}

/* uses NULL terminator of charstr */
void str_append_charstr(String *s, char *charstr) {
  str_nappend_charstr(s, charstr, strlen(charstr));
}

void str_append_char(String *s, char c) {
  str_nappend_charstr(s, &c, 1);
}

void str_append_int(String *s, int i) {
  char tmp[STR_SHORT_LEN];
  sprintf(tmp, "%d", i);
  str_append_charstr(s, tmp);
}

/* defaults to 9 digits beyond decimal pt */
void str_append_dbl(String *s, double d) {
  char tmp[(int)ceil(log10(d)) + 1 + 10];
  sprintf(tmp, "%.9f", d);
  str_append_charstr(s, tmp);
}


void str_realloc(String *s, int new_nchars) {
  s->chars = (char*)srealloc(s->chars, new_nchars * sizeof (char));
  s->nchars = new_nchars;
}

/* strcpy for Strings */
void str_cpy(String *dest, String *src) {
  str_clear(dest);
  str_append(dest, src);
}

void str_cpy_charstr(String *dest, char *src) {
  str_clear(dest);
  str_append_charstr(dest, src);
}

void str_ncpy_charstr(String *dest, char *src, int len) {
  str_clear(dest);
  str_nappend_charstr(dest, src, len);
}

/* strdup for Strings */
String *str_dup(String *src) {
  String *s = str_new(src->length);
  str_cpy(s, src);
  return s;
}

/* returns starting index or -1 */
int str_index_of(String *s, String *substr) {
  char *ptr = strstr(s->chars, substr->chars);
  if (ptr == NULL) return -1;
  return ((int)(ptr - s->chars) / sizeof(char)); /* is this right? */
}

void str_substring(String *dest, String *src, int startidx, int len) {
  str_clear(dest);
  if (len < 0 || len > src->length - startidx) 
    len = src->length - startidx;
  str_nappend_charstr(dest, &src->chars[startidx], len);
}

int str_readline(String *s, FILE *F) {
  char buffer[BUFFERSIZE];
  int stop = 0, abort = 0;
  
  str_clear(s);

  do {
    buffer[BUFFERSIZE - 2] = '\n'; 
    if (fgets(buffer, BUFFERSIZE, F) == NULL)
      abort = 1;
    else {
      if (buffer[BUFFERSIZE - 2] == '\n' || buffer[BUFFERSIZE - 2] == '\0') {
    stop = 1;
    /* the penultimate character (immediately preceding the null
       terminator) will NOT be a carriage return or null terminator if
       and only if the length of the line exceeds the size of the
       buffer */
      }
      str_append_charstr(s, buffer);
    }
  } while (!stop && !abort);

  return abort ? EOF : 0;
}

void str_slurp(String *s, FILE *F) {
  String *line = str_new(BUFFERSIZE);
  str_clear(s);
  while (str_readline(line, F) != EOF) 
    str_append(s, line);
  str_free(line);
}

int str_equals(String *s1, String *s2) {
  return (str_compare(s1, s2) == 0);
}

int str_equals_charstr(String *s1, char *s2) {
  return (str_compare_charstr(s1, s2) == 0);
}

int str_equals_nocase(String *s1, String *s2) {
  return (str_compare_nocase(s1, s2) == 0);
}

int str_equals_nocase_charstr(String *s1, char *s2) {
  return (str_compare_nocase_charstr(s1, s2) == 0);
}

int str_compare(String *s1, String *s2) {
  return strcmp(s1->chars, s2->chars);
}

int str_compare_charstr(String *s1, char *s2) {
  return strcmp(s1->chars, s2);
}

int str_compare_nocase(String *s1, String *s2) {
  return strcasecmp(s1->chars, s2->chars);
}

int str_compare_nocase_charstr(String *s1, char *s2) {
  return strcasecmp(s1->chars, s2);
}

void str_trim(String *s) {
  int i;
  for (i = s->length - 1; i >= 0 && isspace(s->chars[i]); i--) 
    s->length--;
  s->chars[s->length] = '\0';
}

void str_double_trim(String *s) {
  int i, j;
  str_trim(s);
  for (j = 0; isspace(s->chars[j]) && j < s->length; j++);
  if (j > 0) {
    s->length -= j;
    for (i = 0; i < s->length; i++) s->chars[i] = s->chars[i+j];
    s->chars[s->length] = '\0';
  }
}

void str_remove_all_whitespace(String *s) {
  int i, j;
  for (i = 0, j = 0; i < s->length; i++) {
    if (!isspace(s->chars[i])) {
      if (i != j)               /* avoid unnecessary assignments */
        s->chars[j] = s->chars[i];
      j++;
    }
  }
  s->length = j;
  s->chars[s->length] = '\0';
}

/** Remove quotes from around string, if they exist */
void str_remove_quotes(String *str) {
  if ((str->chars[0] == '\"' || str->chars[0] == '\'') &&
      (str->chars[str->length-1] == '\"' || str->chars[str->length-1] == '\'')) {
    str->chars[0] = ' ';
    str->chars[str->length-1] = ' ';
    str_double_trim(str);
  }
}

int str_split(String *s, char* delim, List *l) {
  int i, j, n;
  int inv_delim[NCHARS];
  String *tok;
  char *real_delim;

  if (delim == NULL)        /* whitespace */
    real_delim = " \t\n\r\f\v"; 
  else
    real_delim = delim;

  /* prepare inv_delim */
  for (i = 0; i < NCHARS; i++) inv_delim[i] = 0;
  for (i = 0; real_delim[i] != '\0'; i++) inv_delim[(int)real_delim[i]] = 1;

  lst_clear(l);
  n = 0;
  for (i = 0; i < s->length; i += n + 1) {
    for (j = i; j < s->length; j++) if (inv_delim[(int)s->chars[j]]) break;
    n = j - i;
    tok = str_new(n);
    str_substring(tok, s, i, n);
    lst_push_ptr(l, tok);       

    if (delim == NULL)      /* gobble whitespace */
      for (j++; j < s->length; j++, n++) 
    if (!inv_delim[(int)s->chars[j]]) break;
  }
  return lst_size(l);
}

int str_as_int(String *s, int *i) {
  char *endptr;
  int tmp = (int)strtol(s->chars, &endptr, 0);
  if (endptr == s->chars) return 1;
  *i = tmp;
  return (endptr - s->chars == s->length ? 0 : 2);
}

int str_as_dbl(String *s, double *d) {
  char *endptr;
  double tmp = (double)strtod(s->chars, &endptr);
  if (endptr == s->chars) return 1;
  *d = tmp;
  return (endptr - s->chars == s->length ? 0 : 2);
}

int str_starts_with(String *s, String *substr) {
  if (substr->length > s->length) return 0;
  return (strncmp(s->chars, substr->chars, substr->length) == 0);
}

int str_starts_with_charstr(String *s, char *substr) {
  int len = strlen(substr);
  if (len > s->length) return 0;
  return (strncmp(s->chars, substr, strlen(substr)) == 0);
}


Regex *str_re_new(char *re_str) {
  const char *retval;
  Regex *re = (Regex*)smalloc(sizeof(Regex));
  re->translate = 0;
  re->fastmap = 0;
  re->buffer = 0;
  re->allocated = 0;

  re_syntax_options =  RE_SYNTAX_EGREP; 
                                /* use egrep-style regex conventions
                                   (see GNU regex manual). */ 

  if ((retval = re_compile_pattern(re_str, strlen(re_str), re)) != NULL) {
    fprintf(stderr, "ERROR: cannot compile regular expression '%s'\n\n\
Error message is as follows: %s\n.", re_str, retval);
    exit(-1);
  }

  return re;
}

void str_re_free(Regex *re) {
  regfree(re);
  free(re);
}

int str_re_match(String *s, Regex *re, List *l, int nsubexp) {
  static struct re_registers regs;
  int retval, i;
  String *substr;

  if (l != NULL)
    lst_clear(l);

  retval = re_match(re, s->chars, s->length, 0, &regs);

  if (retval > 0 && l != NULL) {
    for (i = 0; i <= nsubexp; i++) {
      if (regs.start[i] == -1)
        lst_push_ptr(l, NULL);
      else {
        substr = str_new(regs.end[i] - regs.start[i]);
        str_substring(substr, s, regs.start[i],
                      regs.end[i] - regs.start[i]);
        lst_push_ptr(l, substr);
      }
    }
  }

  return retval;
}

int str_re_search(String *s, Regex *re, int start_offset, List *l, 
                  int nsubexp) {
  static struct re_registers regs;
  int retval, i;
  String *substr;

  if (l != NULL)
    lst_clear(l);

  retval = re_search(re, s->chars, s->length, start_offset, s->length, &regs);

  if (retval >= 0 && l != NULL) {
    for (i = 0; i <= nsubexp; i++) {
      if (regs.start[i] == -1)
        lst_push_ptr(l, NULL);
      else {
        substr = str_new(regs.end[i] - regs.start[i]);
        str_substring(substr, s, regs.start[i],
                      regs.end[i] - regs.start[i]);
        lst_push_ptr(l, substr);
      }
    }
  }

  return retval;
}

void str_re_split(String *s, Regex *re, List *l) {
  List *delim_l;
  String *tok;
  int tok_idx, delim_idx, j, tok_len, delim_len;

  delim_l = lst_new_ptr(1);     /* only care about match to entire
                                   expression (0th register) */
  lst_clear(l);

  for (tok_idx = 0; (delim_idx = str_re_search(s, re, tok_idx, delim_l, 0)) > 0;  
       tok_idx += tok_len + delim_len) {

    tok_len = delim_idx - tok_idx;
    delim_len = ((String*)lst_get_ptr(delim_l, 0))->length;

    tok = str_new(tok_len);
    str_substring(tok, s, tok_idx, tok_len);
    lst_push_ptr(l, tok);

    /* the lines below are a precaution against the possibility that the
       user puts parentheses in the delimiter regex (generally won't
       be necessary or even sensible); we have to be sure we free ALL
       allocated substrings */
    for (j = 0; j < lst_size(delim_l); j++)
      str_free((String*)lst_get_ptr(delim_l, j));
  }

  /* don't forget last token */
  if (tok_idx < s->length) {
    tok = str_new(s->length - tok_idx);
    str_substring(tok, s, tok_idx, s->length);
    lst_push_ptr(l, tok);
  }    

  lst_free(delim_l);
}

/* Extracts portion of String before the final '.' character.  Useful
   for filenames.  If no '.' character is present, root is equal to
   the original string. */
void str_get_name_root(String *prefix, String *src) {
  int i;
  for (i = src->length - 1; i >= 0 && src->chars[i] != '.'; i--);
  if (i == 0) str_cpy(prefix, src);
  else str_substring(prefix, src, 0, i);
}

/* Extracts portion of String after the final '.' character.  Useful
   for filenames.  If no '.' character is present, suffix is the empty
   string */
void str_get_name_suffix(String *prefix, String *src) {
  int i;
  for (i = src->length - 1; i >= 0 && src->chars[i] != '.'; i--);
  if (i == 0) str_clear(prefix);
  else str_substring(prefix, src, i+1, src->length - (i+1));
}

/* Removes portion of String prior to the final '/' or '\' character.
   Useful for filenames.  If no '/' or '\' character is present,
   string remains unchanged */
void str_remove_path(String *str) {
  int i, offset;
  for (i = str->length - 1; i >= 0 && str->chars[i] != '/' && 
         str->chars[i] != '\\'; i--);
  offset = i+1;
  if (offset == 0) return;
  for (i = offset; i < str->length; i++)
    str->chars[i-offset] = str->chars[i];
  str->length -= offset;
  str->chars[str->length] = '\0';
}

/* Returns 1 if list includes specified string, 0 otherwise */
int str_in_list(String *s, List *l) {
  int i;
  for (i = 0; i < lst_size(l); i++)
    if (str_equals(s, lst_get_ptr(l, i))) return 1;
  return 0;
}

/* Computes the overlap of two lists of strings.  Returns 1 if they
   overlap and 0 if they don't.  Designed for short lists (no
   hashing).  The destination list can be one of the source lists. */
int str_list_overlap(List *dest, List *src1, List *src2) {
  int i, j;
  List *l = lst_new_ptr(min(lst_size(src1), lst_size(src2)));
  for (i = 0; i < lst_size(src1); i++) {
    for (j = 0; j < lst_size(src2); j++) {
      if (str_equals(lst_get_ptr(src1, i), lst_get_ptr(src2, j))) {
        lst_push_ptr(l, lst_get_ptr(src1, i));
        break;
      }
    }
  }
  lst_cpy(dest, l);
  lst_free(l);
  return (lst_size(dest) > 0);
}

/* converts a list of strings to a list of integers; aborts on error */
List *str_list_as_int(List *str_list) {
  List *retval = lst_new_int(lst_size(str_list));
  int i, tmp;
  for (i = 0; i < lst_size(str_list); i++) {
    if (str_as_int(lst_get_ptr(str_list, i), &tmp) != 0)
      die("ERROR: bad integer ('%s').\n", 
          ((String*)lst_get_ptr(str_list, i))->chars);
    lst_push_int(retval, tmp);
  }
  return retval;
}

/* converts a list of strings to a list of integers; aborts on error */
List *str_list_as_dbl(List *str_list) {
  List *retval = lst_new_dbl(lst_size(str_list));
  int i;
  double tmp;
  for (i = 0; i < lst_size(str_list); i++) {
    if (str_as_dbl(lst_get_ptr(str_list, i), &tmp) != 0)
      die("ERROR: bad floating-point number ('%s').\n", 
          ((String*)lst_get_ptr(str_list, i))->chars);
    lst_push_dbl(retval, tmp);
  }
  return retval;
}
