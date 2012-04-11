/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* String-handling functions, with automatic memory management and
   basic regex support.
   
   $Id: stringsplus.c,v 1.12 2009-02-19 23:33:48 agd27 Exp $ */

#include <pcre.h>
#include "stringsplus.h"
#include "misc.h"
#include <stdlib.h>
#include <ctype.h>

String *str_new(int starting_nchars) {
  String *s = (String*)smalloc(sizeof(String));
  s->chars = (char*)smalloc((starting_nchars+1) * sizeof(char));
  str_clear(s);
  s->nchars = starting_nchars;
  return s;
}

String *str_new_charstr(const char *str) {
  String *s = str_new(strlen(str));
  str_cpy_charstr(s, str);
  return s;
}

String *str_new_int(int i) {
  char tmp[STR_SHORT_LEN];
  sprintf(tmp, "%d", i);
  return str_new_charstr(tmp);
}

String *str_new_dbl(double d) {
  char tmp[STR_SHORT_LEN];
  sprintf(tmp, "%f", d);
  return str_new_charstr(tmp);
}

void str_free(String *s) {
  sfree(s->chars);
  sfree(s);
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
void str_nappend_charstr(String *s, const char *charstr, int len) {
  int i;
  if (s->length + len > s->nchars)
    str_realloc(s, max(s->length + len, s->nchars * 2));
                                /* try to avoid heavy srealloc when
                                   frequently appending short
                                   strings  */

  for (i = 0; i < len && charstr[i] != '\0'; i++)
    s->chars[s->length + i] = charstr[i];

  s->length += i;
  s->chars[s->length] = '\0';  
}

/* uses NULL terminator of charstr */
void str_append_charstr(String *s, const char *charstr) {
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
  s->chars = (char*)srealloc(s->chars, (new_nchars+1) * sizeof (char));
  s->nchars = new_nchars;
}

/* strcpy for Strings */
void str_cpy(String *dest, String *src) {
  str_clear(dest);
  str_append(dest, src);
}

void str_cpy_charstr(String *dest, const char *src) {
  str_clear(dest);
  str_append_charstr(dest, src);
}

void str_ncpy_charstr(String *dest, const char *src, int len) {
  str_clear(dest);
  str_nappend_charstr(dest, src, len);
}

/* copy_charstr for Strings */
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

/* Peek at next line in file */
int str_peek_next_line(String *s, FILE *F) {
  char buffer[BUFFERSIZE];
  int stop = 0, abort = 0, i=0, buffer_used=0;
  
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
    //Determine how many characters of the buffer were used
    buffer_used=-1;
    i=0;
    while (buffer[i] != '\0')
    {  
      buffer_used++; 
      i++;
    }
    //write back the characters we read since this is only peek
    ungetc('\n', F);
    for(i=buffer_used;i>0; i--)
    {
      ungetc(buffer[i-1], F);
    }
  } while (!stop && !abort);
  return abort ? EOF : 0;
}



/* Recursively peek at line L in file */
int str_peek_line(String *s, FILE *F, int L) {
  char buffer[BUFFERSIZE];
  int stop = 0, abort = 0, i=0, buffer_used=0;
  
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
       if( L > 1)
        {
          str_clear(s);
          str_peek_line(s, F, L-1);
        }else
          str_append_charstr(s, buffer);
      } else 
         str_append_charstr(s, buffer);
      
    }
    //Determine how many characters of the buffer were used
    buffer_used=-1;
    i=0;
    while (buffer[i] != '\0')
    {  
      buffer_used++; 
      i++;
    }
    //write back the characters we read since this is only peek
    ungetc('\n', F);
    for(i=buffer_used;i>0; i--)
    {
      ungetc(buffer[i-1], F);
    }
  } while (!stop && !abort);
  return abort ? EOF : 0;
}


void str_substring(String *dest, String *src, int startidx, int len) {
  str_clear(dest);
  if (startidx < 0 || startidx >= src->length)
    die("ERROR in str_substr: startidx is outside the coordinates of the source string!\n");
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

int str_equals_charstr(String *s1, const char *s2) {
  return (str_compare_charstr(s1, s2) == 0);
}

int str_equals_nocase(String *s1, String *s2) {
  return (str_compare_nocase(s1, s2) == 0);
}

int str_equals_nocase_charstr(String *s1, const char *s2) {
  return (str_compare_nocase_charstr(s1, s2) == 0);
}

int str_compare(String *s1, String *s2) {
  return strcmp(s1->chars, s2->chars);
}

int str_compare_charstr(String *s1, const char *s2) {
  return strcmp(s1->chars, s2);
}

int str_compare_nocase(String *s1, String *s2) {
  return strcasecmp(s1->chars, s2->chars);
}

int str_compare_nocase_charstr(String *s1, const char *s2) {
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

/** Remove quotes from around string, if they exist.  Also removes
    whitespace immediately inside quotes */
void str_remove_quotes(String *str) {
  if ((str->chars[0] == '\"' || str->chars[0] == '\'') &&
      (str->chars[str->length-1] == '\"' || str->chars[str->length-1] == '\'')) {
    str->chars[0] = ' ';
    str->chars[str->length-1] = ' ';
    str_double_trim(str);
  }
}


/** Split on whitespace, but don't split strings that fall inside quotes */
int str_split_with_quotes(String *s, const char *delim, List *l) {
  const char *real_delim = " \t\n\r\f\v";
  char *quote;
  int quoteIdx=-1, inv_delim[NCHARS], i, j, n;
  String *tok;

  lst_clear(l);
  if (s->length == 0) return 0;

  if (delim == NULL)        /* whitespace */
    real_delim = " \t\n\r\f\v"; 
  else
    real_delim = delim;

  /* prepare inv_delim */
  for (i = 0; i < NCHARS; i++) inv_delim[i] = 0;
  for (i = 0; real_delim[i] != '\0'; i++) {
    if (real_delim[i] == '"' || real_delim[i] == '\'')
      die("str_split_whitespace_with_quotes can't split on quotes");
    inv_delim[(int)real_delim[i]] = 1;
  }
  quote = smalloc(s->length*sizeof(char));

  n = 0;
  for (i = 0 ; i < s->length; i += n+1) {
    for (j=i; j < s->length; j++) {
      if (s->chars[j] == '"' ||
	  s->chars[j] == '\'') {
	if (quoteIdx >= 0 && quote[quoteIdx]==s->chars[j])
	  quoteIdx--;
	else quote[++quoteIdx] = s->chars[j];
      } else 
	if (quoteIdx == -1 && inv_delim[(int)s->chars[j]]) break;
    }
    n = j - i;
    tok = str_new(n);
    str_substring(tok, s, i, n);
    lst_push_ptr(l, tok);       

    if (delim == NULL)      /* gobble whitespace */
      for (j++; j < s->length; j++, n++) 
    if (!inv_delim[(int)s->chars[j]]) break;
  }
  sfree(quote);
  return lst_size(l);
}


int str_split(String *s, const char* delim, List *l) {
  int i, j, n;
  int inv_delim[NCHARS];
  String *tok;
  const char *real_delim;

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

int str_starts_with_charstr(String *s, const char *substr) {
  int len = strlen(substr);
  if (len > s->length) return 0;
  return (strncmp(s->chars, substr, len) == 0);
}

int str_ends_with(String *s, String *substr) {
  if (substr->length > s->length) return 0;
  return (strncmp(&s->chars[s->length - substr->length], 
                  substr->chars, substr->length) == 0);
}

int str_ends_with_charstr(String *s, const char *substr) {
  int len = strlen(substr);
  if (len > s->length) return 0;
  return (strncmp(&s->chars[s->length - len], substr, len) == 0);
}

Regex *str_re_new(const char *re_str) {
  Regex *re;
  const char *errstr;
  int erroffset;

  re = pcre_compile(re_str, 0, &errstr, &erroffset, NULL);
  if (re == NULL) {
    die("ERROR: cannot compile regular expression '%s' (%d): %s\n",
	re_str, erroffset, errstr);
  }
  return re;
}


//NOTE Regex are allocated by pcre; do not use sfree
void str_re_free(Regex *re) {
  if (re != NULL)
    free(re);
}


#define OVECCOUNT 300
int str_re_match_sub(String *s, Regex *re, List *l, int offset, int nsubexp, 
		     int *first_match) {
  int i, len, rc, ovector[OVECCOUNT], rv;
  String *substr;

  /* WARNING: lst_clear DOES NOT free memory associated with the contents,
     so must free substrings from previous calls if these are no longer being
     used or there will be a memory leak! */
  if (l != NULL) lst_clear(l);

  rc = pcre_exec(re, NULL, s->chars, s->length, offset, 0, ovector, OVECCOUNT);
  if (rc == PCRE_ERROR_NOMATCH) return -1;
  if (rc < 0) return -2;  //any other error
  if (first_match != NULL) (*first_match) = ovector[0];
  rv = ovector[1]-ovector[0];
  if (rc >= 0 && l != NULL) {
    if (rc == 0) {
      printf("nsubexp=%i rc=%i\n", nsubexp, rc);
      fprintf(stderr, "Warning: pcre_exec only has room for %d captured substrings.  May need to increase OVECCOUNT and re-compile\n", OVECCOUNT/3);
      rc = OVECCOUNT/3;
    }
    for (i = 0; i < rc && i <= nsubexp; i++) {
      if (ovector[2*i]==-1) {
	if (ovector[2*i+1] != -1)
	  die("ERROR str_re_match_sub expected ovector[%i]==-1, got %i\n",
	      2*i+1, ovector[2*i+1]);
	lst_push_ptr(l, NULL);
      } else {
	len = ovector[2*i+1] - ovector[2*i];
	substr = str_new(len);
	str_substring(substr, s, ovector[2*i], len);
	lst_push_ptr(l, substr);
      }
    }
  }
  return rv;
}


int str_re_match(String *s, Regex *re, List *l, int nsubexp) {
  return str_re_match_sub(s, re, l, 0, nsubexp, NULL);
}

int str_re_search(String *s, Regex *re, int start_offset, List *l,
		  int nsubexp) {
  int first_match_idx, rc;
  rc = str_re_match_sub(s, re, l, start_offset, nsubexp, &first_match_idx);
  if (rc < 0) return rc;
  return first_match_idx;
}


/* Reduces to portion of String before final instance of specified
   delimiter.  Useful for filenames.  If no delimiter is present,
   string is left unchanged. */
void str_root(String *str, char delim) {
  int i;
  for (i = str->length - 1; i >= 0 && str->chars[i] != delim; i--);
  if (i < 0) return;
  str->length = i;
  str->chars[str->length] = '\0';
}

/* As str_root but finds the first occurence of delimiter.
   Useful for spliting species name and chromosome in MAF files. */
void str_shortest_root(String *str, char delim) {
  int i;
  for (i = 0; i < str->length && str->chars[i] != delim; i++);
  if (i == str->length) return;
  str->length = i;
  str->chars[str->length] = '\0';
}


/* Reduces to portion of String after final instance of specified
   delimiter .  Useful for filenames.  If no delimiter is present,
   suffix is the empty string */
void str_suffix(String *str, char delim) {
  int i, j;
  for (i = str->length - 1; i >= 0 && str->chars[i] != delim; i--);
  if (i < 0) str_clear(str);
  else {
    str->length -= (i + 1);
    for (j = 0; j < str->length; j++) str->chars[j] = str->chars[j+i+1];
    str->chars[str->length] = '\0';
  }
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

/* Returns 1 if list includes specified string, 0 otherwise.  Sets
   *idx to index of element in list if found */
int str_in_list_idx(String *s, List *l, int *idx) {
  int i;
  *idx = -1;
  for (i = 0; i < lst_size(l); i++) {
    if (str_equals(s, lst_get_ptr(l, i))) {
      *idx = i;
      return 1;
    }
  }
  return 0;
}

/* Returns 1 if list includes specified string, 0 otherwise */
int str_in_list_charstr(const char *s, List *l) {
  int i;
  for (i = 0; i < lst_size(l); i++)
    if (str_equals_charstr(lst_get_ptr(l, i), s)) return 1;
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
  int i, tmp=0;
  for (i = 0; i < lst_size(str_list); i++) {
    if (str_as_int(lst_get_ptr(str_list, i), &tmp) != 0)
      die("ERROR: bad integer ('%s').\n", 
          ((String*)lst_get_ptr(str_list, i))->chars);
    lst_push_int(retval, tmp);
  }
  return retval;
}

/* converts a list of strings to a list of doubles; aborts on error */
List *str_list_as_dbl(List *str_list) {
  List *retval = lst_new_dbl(lst_size(str_list));
  int i;
  double tmp=0;
  for (i = 0; i < lst_size(str_list); i++) {
    if (str_as_dbl(lst_get_ptr(str_list, i), &tmp) != 0)
      die("ERROR: bad floating-point number ('%s').\n", 
          ((String*)lst_get_ptr(str_list, i))->chars);
    lst_push_dbl(retval, tmp);
  }
  return retval;
}

void str_toupper(String *s) {
  int i;
  for (i = 0; i < s->length; i++)
    s->chars[i] = toupper(s->chars[i]);
}

void str_tolower(String *s) {
  int i;
  for (i = 0; i < s->length; i++)
    s->chars[i] = tolower(s->chars[i]);
}
