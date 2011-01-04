#include <lists.h>
#include <tree_model.h>
/* Recursive list structure emulating R's lists.  List of particular list types:
   char*, integer, or double currently supported */

/** @file list_of_lists.h
   Recursive list structure and supporting functions.

   Recursive list structure emulating R's lists.  List of particular list types:
   char*, integer, or double currently supported.  Makes use of singular list in lists.h.
   @ingroup base
   @see lists.h
*/


#ifndef __LIST_OF_LISTS__
#define __LIST_OF_LISTS__
/** List types supported */
typedef enum { INT_LIST, /**< List of Integers */
	       DBL_LIST, /**< List of Doubles */
	       CHAR_LIST,  /**< List of Chars */
		   MSA_PTR_LIST, /**< List of pointers to MSAs */
		   GFF_PTR_LIST, /**< List of pointers to GFFs */
	       LIST_LIST  /**< List of Lists */
	      } list_element_type;

/** Basic List of Lists object used recursively to create List of Lists  
*/
struct list_of_list_struct {
  List *lst;      /**< list of lists;  Each element can be be a List * 
                  containing ints, doubles, or char*, OR an element can
                  be a ListOfList*. */
  List *lstName;  /**< list of char* giving name of each element of lst */
  List *lstType;  /**< list of list_element_type giving type of each element in lst */
  char *class;    /**< NULL if this should be treated as a list; otherwise tells
                  R to coerce this list into this type of element (ie,
                  "matrix", "data.frame", "tm") */
};

typedef struct list_of_list_struct ListOfLists;

/** \name List of List allocation/cleanup functions 
 \{ */

/** Create a new list of lists
  @param approx_size approximate size, this should be as close to the maximum size your list will be so the list doesn't have to be expanded which is computationally intensive
  @result newly created list of lists of size specified
*/
ListOfLists *lol_new(int approx_size);

/** Frees a list of lists
  @param[in,out] lol List of Lists
  @warning uses the lst_free method in lists.h for all lists other than chars
*/
void lol_free(ListOfLists *lol);


/** \name List of List push functions 
 \{ */

/** Add an object to the end of a list of lists object
  @param[in,out] lol List of Lists
  @param[in] data object to be added to list of lists
  @param[in] name name of object being added
  @param[in] listType what type of object is being added 
*/
void lol_push(ListOfLists *lol, void *data,
	      const char *name, 
	      list_element_type listType);

/** Add a list to the end of a list of lists object
  @param[in,out] lol List of Lists
  @param[in] data List to be added to list of lists
  @param[in] name name of list being added
  @param[in] listType what type of list is being added 
  @see lol_push
*/
void lol_push_list(ListOfLists *lol, List *data, 
		   const char *name, list_element_type listType);

/** Add a list to the end of a list of lists object
  @param[in,out] lol List of Lists
  @param[in] data List of Lists to be added to list of lists
  @param[in] name name of list of lists being added
  @see lol_push
*/
void lol_push_lol(ListOfLists *lol, 
		  ListOfLists *data, const char *name);

/** Add a list of doubles to end of list of lists
  @param[in,out] lol List of Lists
  @param[in] vals Array of doubles to be added to list of lists
  @param[in] len number of doubles to be added to the list
  @param[in] name name of list of doubles being added
  @see lol_push
*/
void lol_push_dbl(ListOfLists *lol, double *vals, int len,
		  const char *name);
/** Add a list of integers to end of list of lists
  @param[in,out] lol List of Lists
  @param[in] vals Array of integers to be added to list of lists
  @param[in] len number of integers to be added to the list
  @param[in] name name of list of integers being added
  @see lol_push
*/
void lol_push_int(ListOfLists *lol, int *vals, int len, 
		  const char *name);

/** Add a character vector to end of list of lists
  @param[in,out] lol List of Lists
  @param[in] vals character vector to be added to list of lists
  @param[in] len length of character vector to be added to the list
  @param[in] name name of character vector being added
  @see lol_push
*/
void lol_push_charvec(ListOfLists *lol, char **vals, int len,
		      const char *name);

/** Add a matrix to end of list of lists
  @param[in,out] lol List of Lists
  @param[in] mat Matrix to be added to list of lists
  @param[in] name name of the Matrix being added
  @see lol_push
*/
void lol_push_matrix(ListOfLists *lol, Matrix *mat, 
		     const char *name);

/** Add an MSA to end of list of lists
  @param[in,out] lol List of Lists
  @param[in] msa msa to be added to list of lists
  @param[in] name name of msa being added
  @see lol_push
*/
void lol_push_msa(ListOfLists *lol, MSA *msa, const char *name);

/** Add a Tree Model to end of list of lists
  @param[in,out] lol List of Lists
  @param[in] tm Tree Model to be added to list of lists
  @param[in] name name of Tree Model being added
  @see lol_push
*/
void lol_push_treeModel(ListOfLists *lol, TreeModel *tm,
			const char *name);

/** Add an GFF to end of list of lists
  @param[in,out] lol List of Lists
  @param[in] gff gff to be added to list of lists
  @param[in] name name of gff being added
  @see lol_push
*/
void lol_push_gff(ListOfLists *lol, GFF_Set *gff,
		  const char *name);

/** Add a WIG to end of list of lists
  @param[in,out] lol List of Lists
  @param[in] scores Array of conservation scores to be added to list of lists along with corresponding msa
  @param[in] msa msa to be added to list of lists along with corresponding scores
  @see lol_push
*/
void lol_push_wig(ListOfLists *lol, double *scores, MSA *msa);

/** Add a GFF (via pointer) to end of list of lists 
  @param[in,out] lol List of Lists
  @param[in] gff Feature Set to add
  @param[in] name Name to refer to Feature Set by
*/
void lol_push_gff_ptr(ListOfLists *lol, GFF_Set *gff, const char *name);

/** Add a MSA (via pointer) to end of list of lists 
  @param lol[in,out] lol List of Lists
  @param[in] msa Multiple Alignment to add
  @param[in] name Name to refer to MSA by
*/
void lol_push_msa_ptr(ListOfLists *lol, MSA *msa, const char *name);

/** Add a Double Array to the end of list of lists 
   @param lol[in,out] lol List of Lists
   @param data[in] Array of doubles to add 
   @param name[in] Name to refer to Double Array by
   @param ndim[in] Number of dimensions
   @param dimsize[in] Dimension size
   @param dimname[in] Dimension name
*/
void lol_push_dbl_array(ListOfLists *lol, void *data, char *name, 
			int ndim, int *dimsize, char ***dimname);


/** \} \name List Of Lists RPHAST functions.
    \{ */

/** Set a class of list of lists
  @param[in,out] lol List of Lists
  @param[in] class what type of data R should coerce this list into
  @note NULL if this is simply a list
*/
void lol_set_class(ListOfLists *lol, char *class);



/** \} \name List of Lists Search by Name
\{ */

/** Find a list by name & type.
    @param lol[in] List of Lists to search for list within
	@param lstName[in] Name of list to retrieve
	@param lstType[in] Type of list to retrieve i.e. MSA_PTR_LIST
	@result List identified by name if found, or NULL if not found
	@warning If multiple lists have the same name program will die
*/
List *lol_find_list(ListOfLists *lol, const char *lstName, 
		    list_element_type lstType);

/** Find a List of Lists within a List of Lists.
    @param lol List of Lists to search within
	@param name Name of List of Lists to retrieve
	@result List of Lists identified by name if found, or NULL if not found
	@warning If multiple lists have the same name program will die
*/
ListOfLists *lol_find_lol(ListOfLists *lol, const char *name);
/** \} */



#endif
