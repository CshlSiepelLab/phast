#include <lists.h>
/* Recursive list structure emulating R's lists.  List of basic list types:
   char*, integer, or double currently supported */

#ifndef __LIST_OF_LISTS__
#define __LIST_OF_LISTS__

typedef enum { INT_LIST, DBL_LIST, CHAR_LIST, LIST_LIST } list_element_type;

struct list_of_list_struct {
  List *lst;      //list of lists.  Each element can be be a List * 
                  //containing ints, doubles, or char*.  OR an element can
                  //be a ListOfList*.
  List *lstName;  //list of char* giving name of each element of lst
  List *lstType;  //list of list_element_type giving type of each element in lst
  int isMatrix;      //if isMatrix or isDataFrame is TRUE, all lsts 
                     //should be same length.
  int isDataFrame;
};

typedef struct list_of_list_struct ListOfLists;

ListOfLists *ListOfLists_new(int approx_size);
void ListOfLists_push(ListOfLists *lol, void *data,
		     const char *name, 
		     list_element_type listType);

void ListOfLists_push_list(ListOfLists *lol, List *data, 
			   const char *name, list_element_type listType);

void ListOfLists_push_listOfLists(ListOfLists *lol, 
				  ListOfLists *data, const char *name);

void ListOfLists_push_dbl(ListOfLists *lol, double *vals, int len,
			    const char *name);
void ListOfLists_push_int(ListOfLists *lol, int *vals, int len, 
			  const char *name);

void ListOfLists_push_charvec(ListOfLists *lol, char **vals, int len,
			      const char *name);

void ListOfLists_free(ListOfLists *lol);

#endif
