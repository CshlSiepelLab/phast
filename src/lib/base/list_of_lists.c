#include <stdlib.h>
#include <misc.h>
#include <list_of_lists.h>

//make a new list of lists (LOL)
ListOfLists *ListOfLists_new(int approx_size) {
  ListOfLists *lol = malloc(sizeof(ListOfLists));
  lol->lst = lst_new_ptr(approx_size);
  lol->lstName = lst_new_ptr(approx_size);
  lol->lstType = lst_new_int(approx_size);
  lol->isMatrix = lol->isDataFrame = 0;
  return lol;
}

//adds a pointer to list to end of LOL object.
void ListOfLists_push(ListOfLists *lol, void *data, const char *name, 
		     list_element_type listType) {
  char *tempname;
  lst_push_ptr(lol->lst, (void*)data);
  if (name == NULL)
    tempname = NULL;
  else {
    tempname = smalloc((strlen(name)+1)*sizeof(char));
    strcpy(tempname, name);
  }
  lst_push_ptr(lol->lstName, (void*)tempname);
  lst_push_int(lol->lstType, listType);
}


void ListOfLists_push_list(ListOfLists *lol, List *data, const char *name,
			   list_element_type listType) {
  ListOfLists_push(lol, (void*)data, name, listType);
}

void ListOfLists_push_listOfLists(ListOfLists *lol, ListOfLists *data,
				  const char *name) {
  ListOfLists_push(lol, (void*)data, name, LIST_LIST);
}


//add a list of doubles to end of LOL object.  Copies all values.
void ListOfLists_push_dbl(ListOfLists *lol, double *vals, int len,
			  const char *name) {
  List *lst = lst_new_dbl(len);
  int i;
  for (i=0; i<len;  i++) 
    lst_push_dbl(lst, vals[i]);
  ListOfLists_push(lol, lst, name, DBL_LIST);
}



//add a list of ints to end of LOL object.  Copies all values.
void ListOfLists_push_int(ListOfLists *lol, int *vals, int len,
			const char *name) {
  List *lst = lst_new_dbl(len);
  int i;
  for (i=0; i<len;  i++) 
    lst_push_int(lst, vals[i]);
  ListOfLists_push(lol, lst, name, INT_LIST);
}


//add a list of char* to end of LOL object.  Copies all values.
void ListOfLists_push_charvec(ListOfLists *lol, char **vals, int len,
			      const char *name) {
  List *lst = lst_new_ptr(len);
  char *tempstr;
  int i;
  for (i=0; i<len; i++) {
    tempstr = smalloc((strlen(vals[i])+1)*sizeof(char));
    strcpy(tempstr, vals[i]);
    lst_push_ptr(lst, tempstr);
  }
  ListOfLists_push(lol, lst, name, CHAR_LIST);
}


//free a lol and all associated memory (if type is char free the strings)
void ListOfLists_free(ListOfLists *lol) {
  List *currlst;
  int i, j;
  list_element_type currtype;
  char *currname=NULL, *currstr;

  for (i=0; i<lst_size(lol->lst); i++) {
    currlst = (List*)lst_get_ptr(lol->lst, i);
    currtype = lst_get_int(lol->lstType, i);
    currname = (char*)lst_get_ptr(lol->lstName, i);
    if (currname != NULL)
      free(currname);
    if (currtype == LIST_LIST) 
      ListOfLists_free((ListOfLists*)currlst);
    else {
      if (currtype == CHAR_LIST) {
	for (j=0; j<lst_size(currlst); j++) {
	  currstr = (char*)lst_get_ptr(currlst, j);
	  free(currstr);
	}
      }
      lst_free(currlst);
    }
  }
  lst_free(lol->lstType);
  lst_free(lol->lstName);
  free(lol);
}
