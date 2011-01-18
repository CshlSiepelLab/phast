/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file stringsplus.h
    String-handling functions, with automatic memory management and basic regex support. 
    @ingroup base
*/
   
/*
   TODO
   * support "views" of strings (use same memory as another string).
   How to handle null terminators?
   * add support for GNU's regex "fastmaps"
   * add regex substitution
*/


#ifndef STRINGSPLUS_H
#define STRINGSPLUS_H

#include <pcre.h>
#include "lists.h"
#include "stdio.h"

/** Size of line-by-line buffer;
overflow okay (handled by multiple reads) */
#define BUFFERSIZE 100000	

/** Number of ASCII chars */
#define NCHARS 256		


/* for convenience when initializing strings */
/** Length of short string */
#define STR_SHORT_LEN 50
/** Length of medium string */
#define STR_MED_LEN 256
/** Length of long string */
#define STR_LONG_LEN 1000
/** Length of very long string */
#define STR_VERY_LONG_LEN 10000

/** String object */
typedef struct {
  int length;			/**< Length of string in characters */
  char *chars;			/**< Pointer to underlying char* rep; NULL
                           terminator is maintained, so this
                           attribute can be extracted and used
                           as a normal char* string */
  int nchars;			/**< Number of bytes currently allocated */
} String;

/** PCRE is another name for Regex */
typedef pcre Regex;
				
/** \name String Allocate/Cleanup functions 
\{ */

/** Create new String (by size).
   @param starting_nchars Number of characters string can support initially
   @result Newly allocated String object, with specified starting capacity. */
String *str_new(int starting_nchars);

/** Create new String (by char*).  
   @param str char* to create String from
   @result Newly allocated String object, equal to specified char*
   string (string will be copied). */
String *str_new_charstr(const char *str);

/** Create new string (by int).
   @param i Int to create String from
   @result Newly allocated String object, equal to specified int
*/
String *str_new_int(int i);

/** Create new string (by double).
   @param d Double to create String from
   @result Newly allocated String object, equal to specified double
*/
String *str_new_dbl(double d);

/** Free memory associated with string.  
  @param s String to free
  @note String itself will be freed also. 
*/
void str_free(String *s);

/** Clear contents of String.
   @param s String to clear
   @note Contents will be cleared but memory will remain allocated. 
*/
void str_clear(String *s);

/** \} */

/* void str_match(String *s, Regexp *r); */ /* tags? */

/** \name String Append functions */

/** Append one String object to another. 
    @pre Destination String s must be initialized externally.  
    @param s Destination string
    @param suffix String to append
    @note All memory is copied.  
*/
void str_append(String *s, String *suffix);

/** Append a char* string to a String object.  
    @pre Destination String s must be initialized externally.  
    @param s Destination string
    @param charstr *char string to append
    @note All memory is copied.   
*/
void str_append_charstr(String *s, const char *charstr);

/** Append a character to a String object.     
    @pre Destination String s must be initialized externally.  
    @param s Destination string
    @param c Char to append
    @note All memory is copied.  
*/
void str_append_char(String *s, char c);

/** Append first n bytes of a char* string to a String object. 
    @pre Destination String s must be initialized externally.  
    @param s Destination string
    @param charstr Char to append
    @param len Number of characters to append from charstr to String
    @note All memory is copied. 
    @note Operation will not proceed beyond NULL terminator. 
 */
void str_nappend_charstr(String *s, const char *charstr, int len);

/** Append an integer to a String object.
    @pre Destination String s must be initialized externally.  
    @param s Destination string
    @param i Integer to append as a string
    @note All memory is copied. 
  */
void str_append_int(String *s, int i);

/** Append a double to a String object.
    @pre Destination String s must be initialized externally.  
    @param s Destination string
    @param d Double to append as a string
    @note All memory is copied.   */
void str_append_dbl(String *s, double d);

/** Concatenate Strings.  
   @pre Destination String must be initialized externally. 
   @param dest Will be set equal to the concatenation of string src2 to string src1.  
   @param src1 First string to be copied into dest
   @param src2 Second string to be appended to dest
   @note All memory will be copied. */
void str_concat(String *dest, String *src1, String *src2);
/** \} */

/* Reallocate encapsulated char* string.  Intended for internal use
   only. */
void str_realloc(String *s, int new_nchars);

/** \name String copy functions
\{ */

/** Copy a String.     
    @pre Destination String dest must be initialized externally.  
    @param dest Destination string
    @param src String to copy
 */
void str_cpy(String *dest, String *src);

/** Copy a char* string to a String object.  
    @pre Destination String s must be initialized externally.  
    @param dest Destination string
    @param src Source string to copy to dest
    @warning Requires proper NULL terminator
 */
void str_cpy_charstr(String *dest, const char *src);

/** Copy first n bytes of a char* string to a String object.
    @pre Destination String dest must be initialized externally.  
    @param dest Destination string
    @param src Source *char String to copy to dest
    @param len Number of bytes copied from src to dest
    @note Operation will not proceed beyond NULL terminator.  
    @note If necessary, a NULL terminator is added to the 
      String object's chars attribute.  
*/
void str_ncpy_charstr(String *dest, const char *src, int len);

/** Create a duplicate copy of a String.
   @param src String to duplicate
   @result newly allocated String object, equal to specified "template".
 */
String *str_dup(String *src);

/** \}  */

/** Obtain index of first instance of specified substring.
   @param s String to search inside of
   @param substr String to search for
   @result Starting index, or -1 (if substring does not exist). 
*/
int str_index_of(String *s, String *substr);

/** Read a line from a file without modifying file read position (peek). 
   @pre Destination String must be initialized externally.
   @param s String to fill with next line from file
   @param F File descriptor to file to peek into
   @results 0 on success, EOF if end-of-file or error.
*/
int str_peek_next_line(String *s, FILE *F);

/** Read any specified line from a file without modifying file read position (peek).
   @pre Destination String must be initialized externally.
   @param s String to fill with the line from file
   @param F File descriptor to file to peek into
   @param L Line number to read from file
   @results 0 on success, EOF if end-of-file or error.
*/
int str_peek_line(String *s, FILE *F, int L);

/** Obtain substring.  
   @pre Destination String s must be initialized externally.  
   @param dest Destination string
   @param src Source String to copy from
   @param startidx Index of Source String to start copying from
   @param len Length of substring to copy

   @note String dest will be set equal to substring of
   length len starting at index startidx.  
   @note If indices are out of bounds, dest will be
   set equal to the substring from startidx to the end of src. */
void str_substring(String *dest, String *src, int startidx, int len);

/** Read a line from a file.
    @pre Destination String s must be initialized externally.  
    @param s Destination string
    @param F File to read from
    @result 0 on success and EOF if end-of-file or error
*/
int str_readline(String *s, FILE *F); 

/** Read an entire file.   
   @pre Destination String s must be initialized externally.  
   @param s Destination string
   @param F File to read from     
 */
void str_slurp(String *s, FILE *F); 

/** \name String comparison functions 
\{ */

/** Test if two String objects are equal 
   @param s1 String to compare
   @param s2 String to compare to
   @result 1 if equal, 0 otherwise
   @note uses str_compare 
*/
int str_equals(String *s1, String *s2);

/** Test if a String object is equal to a char* string 
   @param s1 String to compare
   @param s2 *char string to compare to
   @result 1 if equal, 0 otherwise. 
   @note uses str_compare_charstr
*/
int str_equals_charstr(String *s1, const char *s2);

/** Test if two String objects are equal, ignoring case 
   @param s1 String to compare
   @param s2 String to compare to
   @result 1 if equal, 0 otherwise.
   @note uses str_compare_nocase 
*/
int str_equals_nocase(String *s1, String *s2);

/** Test if a String object is equal to a char* string, ignoring case
   @param s1 String to compare
   @param s2 *char string to compare to
   @result 1 if equal, 0 otherwise. 
   @note uses str_compare_nocase_charstr */
int str_equals_nocase_charstr(String *s1, const char *s2);

/** Compare two String objects alphanumerically.  
   @param s1 String to compare
   @param s2 String to compare to
   @result Negative integer if s1 < s2, Positive integer if s1 > s2,
   or zero if s1 == s2. 
   @note Uses strcmp
*/
int str_compare(String *s1, String *s2);

/** Compare a String object and a char* string alphanumerically.
   @param s1 String to compare
   @param s2 *char String to compare to
   @result Negative integer if s1 < s2, Positive integer if s1 > s2,
   or zero if s1 == s2. 
   @note Uses strcmp */
int str_compare_charstr(String *s1, const char *s2);

/** Compare two String objects alphanumerically, ignoring case. 
   @param s1 String to compare
   @param s2 String to compare to
   @result Negative integer if s1 < s2, Positive integer if s1 > s2,
   or zero if s1 == s2. 
   @note Uses strcmp */
int str_compare_nocase(String *s1, String *s2);

/** Compare a String object and a char* string alphanumerically,
   ignoring case.  
   @param s1 String to compare
   @param s2 *char String to compare to
   @result Negative integer if s1 < s2, Positive integer if s1 > s2,
   or zero if s1 == s2. 
   @note Uses strcasecmp */
int str_compare_nocase_charstr(String *s1, const char *s2);


/** Tests whether a string starts with the specified substring (passed
   as String object).
   @param s String to search for substring in beginning
   @param substr Substring to search for
   @result 0 on exact match and 1 otherwise
 */
int str_starts_with(String *s, String *substr);

/** Tests whether a string starts with the specified substring (passed
   as char* string).
   @param s String to search for substring in beginning
   @param substr char* string to search for 
   @result 0 on exact match and 1 otherwise. 
*/
int str_starts_with_charstr(String *s, const char *substr);

/** Tests whether a string ends with a specified substring
   @param s String to search for substring at end
   @param substr String to search for
   @result 0 on exact match and 1 otherwise
*/
int str_ends_with(String *s, String *substr);

/** Tests whether a string ends with a specified substring
   @param s String to search for substring at end
   @param substr *char String to search for
   @result 0 on exact match and 1 otherwise
*/
int str_ends_with_charstr(String *s, const char *substr);


/** \} */

/** Split on whitespace, unless enclosed in quotes.  
   @note Handles nested double- and single-quotes appropriately 
   @param s String to split
   @param delim Delimiter to split string when encountered
   @param l List containing split up parts of string s at delimiters
   @warning: substrings added to List l are newly allocated and must be
     freed externally.
*/
int str_split_with_quotes(String *s, const char *delim, List *l);

/** Split string on delimiters  
   @pre List l must be allocated externally use lst_new_ptr(...)
   @param s String to split
   @param delim NULL-terminated char* string of valid delimiters
   @param l List to which to add substrings
   @warning: substrings added to List l are newly allocated and must be
     freed externally. */
int str_split(String *s, const char* delim, List *l);

/** \name String Regular Expression (regex) functions 
\{ */
/** Create new regular expression object based on the specified string.
   @note Character string re_str must be NULL terminated.  
   @note Function aborts and reports an error message to stderr if the 
   expression cannot be compiled.  
   @note The PCRE regex package is used, which uses perl regular
   expression syntax.
   @result Newly allocated and compiled Regex object.
 */
Regex *str_re_new(const char *re_str);

/** Free resources associated with regular expression object. 
    @param re Regex object to free
    @note The object itself is freed also. 
*/
void str_re_free(Regex *re);

/** Test whether the specified string matches the specified regex.
   @pre The list 'l' must be initialized externally if non-NULL.  
   @param s String to test if it matches regex
   @param re Regex 
   @param l (Optional) List that holds results of subexpressions (defined with parenthesis). i.e. '22(.)3(..)5' would extract the 3rd digit from the string as first entry and the 5th, and 6th as second entry.
   @param nsubexp Maximum number of subexpressions to examine (not including the 0th one).  
   @result Number of matched characters (possibly zero) on match, -1
   on no match, and -2 on error.
   @note NULLs will be added for all non-matching groups in list 'l'
   @note In the list 'l', the 0th substring corresponds to the entire regex. 
   @note This function uses the pcre_exec function of the PCRE
   regex package.
   @warning Substrings added to List l are newly allocated and must be
   freed externally. */
int str_re_match(String *s, Regex *re, List *l, int nsubexp);

/** Search the specified string for the first instance of the specified
   regex.  
   @pre If used, list 'l' must be allocated externally
   @param start_offset The first start_offset characters will be ignored.
   @param l (Optional) If non-NULL, it will be populated with substrings corresponding
   to subexpressions, as described under str_re_match.  
   @note This function uses the pcre_exec function of the PCRE regex package.
   @result Index of first match, -1 if no match exists, or -2 if an
   internal error occurs. 
   @warning Substrings added to List l are newly allocated and must be
   freed externally. 
   @see str_re_match
*/
int str_re_search(String *s, Regex *re, int start_offset, List *l, 
                  int nsubexp);

/** \} */

/** \name String Modification functions 
\{ */


/** Trim trailing whitespace.  
   All characters c are removed such that isspace(c) != 0.
   @param s String to remove trailing white space from
 */
void str_trim(String* s);

/** Trim leading and trailing whitespace.  
   All characters c are removed such that isspace(c) != 0. 
   @param s String to remove leading and trailing white space from
*/
void str_double_trim(String *s);

/** Remove all whitespace -- leading, trailing, and interior. 
   All characters c are removed such that isspace(c) != 0. 
   @param s String to remove all white space from
*/
void str_remove_all_whitespace(String *s);

/** Remove quotes from string 
  @param str String to remove quotes from
*/
void str_remove_quotes(String *str);


/** Removes portion of String after final instance of specified
   delimiter.
   @param str String to shorten if delimiter exists
   @param delim Keep only characters before last instance of this delimiter
   @code
   //e.g.
   String *tmpstr = str_new_charstr( "/home/user/file.txt");
   str_root(tmpstr, '/');
   //Result = "/home/user"
   @endcode  
   @note Useful for filenames. 
   @note If no delimiter is present, string is left unchanged.
*/
void str_root(String *str, char delim);
/** Removes portion of String after first instance of specified
   delimiter.
   @param str String to shorten if delimiter exists
   @param delim Keep only characters before first instance of this delimiter
   @code
   //e.g.
   String *tmpstr = str_new_charstr( "/home/user/file.txt");
   str_root(tmpstr, 'e');
   //Result = "/hom"
   @endcode  
   @note Useful for splitting species name and chromosome in MAF files.
   @note If no delimiter is present, string is left unchanged.
*/
void str_shortest_root(String *str, char delim);

/** Retain only portion of string after final instance of specified delimiter (if delimiter not present string is cleared)
   @param String to modify
   @param delim Only characters after the final instance of this delimiter are left
   @note Useful for filenames.  
   @note If no delimiter is present, suffix is the empty string 
   @note Very similar to str_root which does not modify string if delimiter is not present
   @see str_root
*/
void str_suffix(String *str, char delim);

/** Remove path from string leaving filename only.
    @param str String containing path to remove
    @note if no '/' or '\' character is present string remains unchanged
    @todo Make sure this works on Windows paths as well as Mac/Linux 
*/
void str_remove_path(String *str);


/** Changes every character in the string to upper case 
    @param s String to make all capitals
*/
void str_toupper(String *s);
/** Changes every character in the string to lower case
   @param s String to make all lower case
*/
void str_tolower(String *s);

/** \} */

/** Test if list contains string
   @param s String to look for in list
   @param l List to look in
   @result 1 if list includes specified string, 0 otherwise
 */
int str_in_list(String *s, List *l);

/** Test if list contains *char string 
   @param s *Char String to look for in list
   @param l List to look in
   @result 1 if list includes specified string, 0 otherwise
*/
int str_in_list_charstr(const char *s, List *l);

/** Computes the overlap of two lists of strings. 
   @param dest List of overlapping strings between source lists  
   @param src1 List of strings 
   @param src2 List of strings
   @result  Returns 1 if they overlap or 0 if they don't
   @note The destination list can be one of the source lists.
   @warning Designed for short lists (no hashing).  
*/
int str_list_overlap(List *dest, List *src1, List *src2);

/** \name String cast functions 
\{ */

/** Attempt to convert a String to an integer.  
   @param s String to convert to an integer
   @param i Integer converted from string
   @result 0 on success, 1 on failure, and 2 on partial success
    (only a prefix of the string could be converted).
 */
int str_as_int(String *s, int *i);

/** Attempt to convert a String to a double.  
   @param s String to convert to a double
   @param d Double converted from string
   @return 0 on success, 1 on failure, and 2 on partial success (only
   a prefix of the string could be converted). */
int str_as_dbl(String *s, double *d);


/** Converts a list of strings to a list of integers; aborts on error.
   @param str_list List of strings that can be converted to integers
   @result List of integers converted from list of strings
*/
List *str_list_as_int(List *str_list);

/** Converts a list of strings to a list of doubles; aborts on error.
   @param str_list List of strings that can be converted to doubles
   @result List of doubles converted from list of strings
*/
List *str_list_as_dbl(List *str_list);
/** \} */

/** Find index of string in list.
  @param[in] s String to search for
  @param[in] l List of strings possibly containing string s
  @param[out] idx Index of string in list (-1 if not found)
  @result 1 if list includes specified string, 0 otherwise. 
 */
int str_in_list_idx(String *s, List *l, int *idx);



#endif
