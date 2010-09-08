/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* stringsplus - String-handling functions, with automatic memory management and basic regex support. */
   
/* $Id: stringsplus.h,v 1.8 2008-11-12 02:07:59 acs Exp $ */

/*
   TODO
   * support "views" of strings (use same memory as another string).
   How to handle null terminators?
   * add support for GNU's regex "fastmaps"
   * add regex substitition
*/


#ifndef STRINGSPLUS_H
#define STRINGSPLUS_H

#include <pcre.h>
#include "lists.h"
#include "stdio.h"

#define BUFFERSIZE 1000		/* size of line-by-line buffer;
                               overflow okay (handled by multiple
                               reads) */

#define NCHARS 256		/* number of ASCII chars */


/* for convenience when initializing strings */
#define STR_SHORT_LEN 50
#define STR_MED_LEN 256
#define STR_LONG_LEN 1000
#define STR_VERY_LONG_LEN 10000

/* String object */
typedef struct {
  int length;			/* length of string */
  char *chars;			/* underlying char* rep; NULL
                           terminator is maintained, so this
                           attribute can be extracted and used
                           as a normal char* string */
  int nchars;			/* number of bytes currently allocated */
} String;


typedef pcre Regex;
				/* COMMENT */

/* Create new String (by size).
   Returns newly allocated String object, with specified starting capacity. */
String *str_new(int starting_nchars);

/* Create new String (by char*).  
   Returns newly allocated String object, equal to specified char*
   string (string will be copied). */
String *str_new_charstr(const char *str);

String *str_new_int(int i);
String *str_new_dbl(double d);

/* Free memory associated with string.  String itself will be freed also. */
void str_free(String *s);

/* Clear contents of String.
   Contents will be cleared but memory will remain allocated. */
void str_clear(String *s);

/* void str_match(String *s, Regexp *r); */ /* tags? */

/* Concatenate Strings.  String dest will be set equal to the
   concatenation of string src2 to string src1.  All memory will be
   copied.  Destination String must be initialized externally. */
void str_concat(String *dest, String *src1, String *src2);

/* Append one String object to another.  All memory is copied.  Destination
   String must be initialized externally. */
void str_append(String *s, String *suffix);

/* Append a char* string to a String object.  All memory is copied.
   Destination String must be initialized externally. */
void str_append_charstr(String *s, const char *charstr);

/* Append a character to a String object.  Destination String must be
   initialized externally. */
void str_append_char(String *s, char c);

/* Append first n bytes of a char* string to a String object.  All
   memory is copied.  Destination String must be initialized
   externally.  Operation will not proceed beyond NULL terminator.  If
   necessary, a NULL terminator is added to the String object's chars
   attribute.  */
void str_nappend_charstr(String *s, const char *charstr, int len);

/* Append a string representation of an integer to a String object.
   Destination String must be initialized externally.  */
void str_append_int(String *s, int i);

/* Append a string representation of an double to a String object.
   Destination String must be initialized externally.  */
void str_append_dbl(String *s, double d);

/* Reallocate encapsulated char* string.  Intended for internal use
   only. */
void str_realloc(String *s, int new_nchars);

/* Copy a String.  Destination String must be initialized externally. */
void str_cpy(String *dest, String *src);

/* Copy a char* string to a String object.  Requires proper NULL
   terminator.  Destination String must be initialized externally.  */
void str_cpy_charstr(String *dest, const char *src);

/* Copy first n bytes of a char* string to a String object.
   Destination String must be initialized externally.  Operation will
   not proceed beyond NULL terminator.  If necessary, a NULL
   terminator is added to the String object's chars attribute.  */
void str_ncpy_charstr(String *dest, const char *src, int len);

/* Create a duplicate copy of a String.
   Returns newly allocated String object, equal to specified "template". */
String *str_dup(String *src);

/* Obtain index of first instance of specified substring.
   Returns starting index or -1 (if substring does not exist). */
int str_index_of(String *s, String *substr);

/* Obtain substring.  String dest will be set equal to substring of
   length len starting at index startidx.  Destination String must be
   initialized externally.  If indices are out of bounds, dest will be
   set equal to the substring from startidx to the end of src. */
void str_substring(String *dest, String *src, int startidx, int len);

/* Read a line from a file.  Destination String must be initialized
   externally.  
   Returns 0 on success and EOF if end-of-file or error. */
int str_readline(String *s, FILE *F); 

/* Read an entire file.   Destination String must be initialized
   externally.  */
void str_slurp(String *s, FILE *F); 

/* Test if two String objects are equal (uses str_compare).
   Returns 1 if equal, 0 otherwise. */
int str_equals(String *s1, String *s2);

/* Test if a String object is equal to a char* string 
   (uses str_compare_charstr).
   Returns 1 if equal, 0 otherwise. */
int str_equals_charstr(String *s1, const char *s2);

/* Test if two String objects are equal, ignoring case (uses
   str_compare_nocase).
   Returns 1 if equal, 0 otherwise. */
int str_equals_nocase(String *s1, String *s2);

/* Test if a String object is equal to a char* string, ignoring case
   (uses str_compare_nocase_charstr).
   Returns 1 if equal, 0 otherwise. */
int str_equals_nocase_charstr(String *s1, const char *s2);

/* Compare two String objects alphanumerically.  Uses strcmp.
   Returns an integer less than zero if s1 < s2, greater than zero if
   s1 > s2, or equal to zero if s1 == s2. */
int str_compare(String *s1, String *s2);

/* Compare a String object and a char* string alphanumerically.  Uses
   strcmp.
   Returns an integer less than zero if s1 < s2, greater than zero if
   s1 > s2, or equal to zero if s1 == s2. */
int str_compare_charstr(String *s1, const char *s2);

/* Compare two String objects alphanumerically, ignoring case.  Uses
   strcasecmp.
   Returns an integer less than zero if s1 < s2, greater than zero if
   s1 > s2, or equal to zero if s1 == s2. */
int str_compare_nocase(String *s1, String *s2);

/* Compare a String object and a char* string alphanumerically,
   ignoring case.  Uses strcasecmp.
   Returns an integer less than zero if s1 < s2, greater than zero if
   s1 > s2, or equal to zero if s1 == s2. */
int str_compare_nocase_charstr(String *s1, const char *s2);

/* Trim trailing whitespace.  All characters c are removed such that
   isspace(c) != 0. */
void str_trim(String* s);

/* Trim leading and trailing whitespace.  All characters c are removed
   such that isspace(c) != 0. */
void str_double_trim(String *s);

/* Remove all whitespace -- leading, trailing, and interior.  All
   characters c are removed such that isspace(c) != 0. */
void str_remove_all_whitespace(String *s);

void str_remove_quotes(String *str);

/* Split on whitespace, unless enclosed in quotes.  Handles
   nested double- and single-quotes appropriately */
int str_split_with_quotes(String *s, const char *delim, List *l);

/* Split a string according to specified delimiters.  Copies of
   resulting substrings will be added to specified List.
   Returns number of fields.

   Warning: substrings added to List l are newly allocated and must be
   freed externally. */
int str_split(String *s,  
              const char* delim,	/* NULL-terminated char* string of
                               valid delimiters or NULL; if NULL,
                               String will be split on
                               whitespace  */
              List *l);	/* List to which to add substrings.
                           Must be allocated externally (use
                           lst_new_ptr(...)).  */

/* Attempt to convert a String to an integer.  
   Returns 0 on success, 1 on failure, and 2 on partial success (only
   a prefix of the string could be converted). */
int str_as_int(String *s, int *i);

/* Attempt to convert a String to a double.  
   Returns 0 on success, 1 on failure, and 2 on partial success (only
   a prefix of the string could be converted). */
int str_as_dbl(String *s, double *d);

/* Tests whether a string starts with the specified substring (passed
   as String object).
   Returns 0 on exact match and 1 otherwise. */
int str_starts_with(String *s, String *substr);

/* Tests whether a string starts with the specified substring (passed
   as char* string).
   Returns 0 on exact match and 1 otherwise. */
int str_starts_with_charstr(String *s, const char *substr);

int str_ends_with(String *s, String *substr);
int str_ends_with_charstr(String *s, const char *substr);

/* Create new regular expression object based on the specified string.
   Character string re_str must be NULL terminated.  Function aborts
   and reports an error message to stderr if the expression cannot be
   compiled.  The PCRE regex package is used, which uses perl regular
   expression syntax.
   Returns a newly allocated and compiled Regex object. */
Regex *str_re_new(const char *re_str);

/* Free resources associated with regular expression object. The
   object itself is freed also. */
void str_re_free(Regex *re);

/* Test whether the specified string matches the specified regex.  If
   List l is non-NULL, it will be populated with substrings
   corresponding to subexpressions (designated with parentheses).  The
   0th such substring corresponds to the entire regex.  A maximum of
   nsubexp will be examined (not including the 0th one).  NULLs will
   be added for all non-matching groups.  The list must be initialized
   externally.  This function uses the pcre_exec function of the PCRE
   regex package.

   Returns number of matched characters (possibly zero) on match, -1
   on no match, and -2 on error.

   Warning: substrings added to List l are newly allocated and must be
   freed externally. */
int str_re_match(String *s, Regex *re, List *l, int nsubexp);

/* Search the specified string for the first instance of the specified
   regex.  The first start_offset characters will be ignored.  If List
   l is non-NULL, it will be populated with substrings corresponding
   to subexpressions, as described under str_re_match.  The list must be
   initialized externally.  This function uses the pcre_exec function
   of the PCRE regex package.
   Returns index of first match, -1 if no match exists, or -2 if an
   internal error occurs. 

   Warning: substrings added to List l are newly allocated and must be
   freed externally. */
int str_re_search(String *s, Regex *re, int start_offset, List *l, 
                  int nsubexp);

void str_root(String *str, char delim);
void str_shortest_root(String *str, char delim);
void str_suffix(String *str, char delim);
void str_remove_path(String *str);
int str_in_list(String *s, List *l);
int str_in_list_charstr(const char *s, List *l);
int str_list_overlap(List *dest, List *src1, List *src2);
List *str_list_as_int(List *str_list);
List *str_list_as_dbl(List *str_list);
int str_in_list_idx(String *s, List *l, int *idx);
void str_toupper(String *s);
void str_tolower(String *s);

#endif
