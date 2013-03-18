#include <lists.h>
#include <tree_model.h>
#include <ms.h>
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
		   MSA_PTR_LIST, /**< List of data from MSAs */
		   GFF_PTR_LIST, /**< List of data from GFFs */
	       LIST_LIST,  /**< List of Lists */
	       MS_PTR_LIST, /**< List of MSs */
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

/** Append an object to a list of lists
  @param[in,out] lol List of Lists
  @param[in] data object to be added to list of lists
  @param[in] name name of object being added
  @param[in] listType what type of object is being added 
*/
void lol_push(ListOfLists *lol, void *data,
	      const char *name, 
	      list_element_type listType);

/** Append a list to list of lists
  @param[in,out] lol List of Lists
  @param[in] data List to be added to list of lists
  @param[in] name name of list being added
  @param[in] listType what type of list is being added 
  @see lol_push
*/
void lol_push_list(ListOfLists *lol, List *data, 
		   const char *name, list_element_type listType);

/** Append a list of lists to list of lists
  @param[in,out] lol List of Lists
  @param[in] data List of Lists to be added to list of lists
  @param[in] name name of list of lists being added
  @see lol_push
*/
void lol_push_lol(ListOfLists *lol, 
		  ListOfLists *data, const char *name);

/** Append a list of doubles to list of lists
  @param[in,out] lol List of Lists
  @param[in] vals Array of doubles to be added to list of lists
  @param[in] len number of doubles to be added to the list
  @param[in] name name of list of doubles being added
  @see lol_push
*/
void lol_push_dbl(ListOfLists *lol, double *vals, int len,
		  const char *name);
/** Append a list of integers to list of lists
  @param[in,out] lol List of Lists
  @param[in] vals Array of integers to be added to list of lists
  @param[in] len number of integers to be added to the list
  @param[in] name name of list of integers being added
  @see lol_push
*/
void lol_push_int(ListOfLists *lol, int *vals, int len, 
		  const char *name);

/** Append a *char string to list of lists
  @param[in,out] lol List of Lists
  @param[in] vals *char string to be added to list of lists
  @param[in] len length of string to be added to the list
  @param[in] name name of string being added
  @see lol_push
*/
void lol_push_charvec(ListOfLists *lol, char **vals, int len,
		      const char *name);

/** Append data from a Matrix to list of lists
  A single list is appended to the list of lists (lol).
  The appended list contains the matrix rows and names
  @code
  //Layout of list to append
  matrixList->
	     //Rows with their column number as name
	      1-> 	  (list)
		  //Data at row 1, column 1 through n
		  2.4	     (double)
		  2.5        (double)
		  ...	    
		   n	     (double)
	      2-> 	  (list)
		  ...
		  ...
		  ...
	      ...
	      n 	  (list)
		  ...
		  ...
		  ...
	      row.names-> (list)
		  //Row numbers
	          1	  (char string)
		  2	  (char string)
		  ...
		  n	  (char string)
  @endcode
  @param[in,out] lol List of Lists
  @param[in] mat Matrix to be added to list of lists
  @param[in] name name of the Matrix being added
  @see lol_push
*/
void lol_push_matrix(ListOfLists *lol, Matrix *mat, 
		     const char *name);

/** Append an MSA to list of lists
  @param[in,out] lol List of Lists
  @param[in] msa msa to be added to list of lists
  @param[in] name name of msa being added
  @see lol_push
*/
void lol_push_msa(ListOfLists *lol, MSA *msa, const char *name);

/** Append data from a Tree Model to list of lists.
  A single list is appended to the list of lists (lol).
  The appended list contains up to 13 lists as shown below (depending on available data)
  @code
  //Layout of list to append
  treeModList->
	      alphabet		(char string)
	      backgd 		(double)
	      rate.matrix 	(matrix)
	      subst.mod		(char string)
	      likelihood	(double)
	      alpha		(double)
	      nratecats		(int)
	      rate.consts	(double)
	      rate.weights	(double)
	      selection		(double)
	      tree		(char string)
	      root.leaf		(int)
	      alt.model->	(list)
		     	  model->	(list) //For each alternative model
			  	subst.mod	(char string)
	      		  	backgd		(double)
			  	rate.matrix	(matrix)
			  	selection	(double)
			  	bgc		(double)
			  	defn		(char string)
  @endcode
  @param[in,out] lol List of Lists
  @param[in] tm Tree Model to be added to list of lists
  @param[in] name name of Tree Model being added
  @see lol_push
*/
void lol_push_treeModel(ListOfLists *lol, TreeModel *tm,
			const char *name);

/** Append data from a GFF to list of lists.
    A single list is appended to the list of lists (lol).
    The appended list contains 5 to 9 lists as shown below (depending on available data)
    @code
    //Layout of list to append
     gffList->
              List of sequence names (char string)
	      List of source   	 (char string)
	      List of feature   	 (char string)
	      List of start     	 (int)
	      List of end       	 (int)
	      List of score	    	 (double)
	      List of strand    	 (char)
	      List of frame     	 (int)
	      List of attribute 	 (char)
  @endcode
  @param[in,out] lol List of Lists
  @param[in] gff Feature Set containing data to add to list of lists
  @param[in] name name of gff being added
  @see lol_push
*/
void lol_push_gff(ListOfLists *lol, GFF_Set *gff,
		  const char *name);

/** Append a GFF (via pointer) to list of lists 
  @param[in,out] lol List of Lists
  @param[in] gff Feature Set to add
  @param[in] name Name to refer to Feature Set by
*/
void lol_push_gff_ptr(ListOfLists *lol, GFF_Set *gff, const char *name);

/** Append a MSA (via pointer) to list of lists 
  @param lol[in,out] lol List of Lists
  @param[in] msa Multiple Alignment to add
  @param[in] name Name to refer to MSA by
*/
void lol_push_msa_ptr(ListOfLists *lol, MSA *msa, const char *name);

/** Append an MS to list of lists
  @param[in,out] lol List of Lists
  @param[in] ms ms to be added to list of lists
  @param[in] name name of ms being added
  @see lol_push
*/
void lol_push_ms_ptr(ListOfLists *lol, MS *ms, const char *name);

/** Append a Double Array to list of lists 
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
