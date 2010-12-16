#include <lists.h>
#include <tree_model.h>
/* Recursive list structure emulating R's lists.  List of particular list types:
   char*, integer, or double currently supported */

#ifndef __LIST_OF_LISTS__
#define __LIST_OF_LISTS__

typedef enum { INT_LIST, DBL_LIST, CHAR_LIST, MSA_PTR_LIST, GFF_PTR_LIST, 
	       LIST_LIST } 
  list_element_type;

struct list_of_list_struct {
  List *lst;      //list of lists.  Each element can be be a List * 
                  //containing ints, doubles, char*, or void*.  OR an element 
                  //can be a ListOfList*.
  List *lstName;  //list of char* giving name of each element of lst
  List *lstType;  //list of list_element_type giving type of each element in lst
  char *class;    //NULL if this should be treated as a list; otherwise tells
                  //R to coerce this list into this type of element (ie,
                  //"matrix", "data.frame", "tm")
};

typedef struct list_of_list_struct ListOfLists;

ListOfLists *lol_new(int approx_size);

void lol_set_class(ListOfLists *lol, char *class);

void lol_push(ListOfLists *lol, void *data,
	      const char *name, 
	      list_element_type listType);

void lol_push_list(ListOfLists *lol, List *data, 
		   const char *name, list_element_type listType);

void lol_push_lol(ListOfLists *lol, 
		  ListOfLists *data, const char *name);

void lol_push_dbl(ListOfLists *lol, double *vals, int len,
		  const char *name);
void lol_push_int(ListOfLists *lol, int *vals, int len, 
		  const char *name);

void lol_push_charvec(ListOfLists *lol, char **vals, int len,
		      const char *name);

void lol_push_matrix(ListOfLists *lol, Matrix *mat, 
		     const char *name);

void lol_push_msa(ListOfLists *lol, MSA *msa, const char *name);

void lol_push_treeModel(ListOfLists *lol, TreeModel *tm,
			const char *name);

void lol_push_gff(ListOfLists *lol, GFF_Set *gff,
		  const char *name);

void lol_push_gff_ptr(ListOfLists *lol, GFF_Set *gff, const char *name);

void lol_push_msa_ptr(ListOfLists *lol, MSA *msa, const char *name);

void lol_push_wig(ListOfLists *lol, double *scores, MSA *msa);

void lol_push_dbl_array(ListOfLists *lol, void *data, char *name, 
			int ndim, int *dimsize, char ***dimname);

List *lol_find_list(ListOfLists *lol, const char *lstName, 
		    list_element_type lstType);

ListOfLists *lol_find_lol(ListOfLists *lol, const char *name);

void lol_free(ListOfLists *lol);

#endif
