#include "misc.h"
#include "list_of_lists.h"
#include <Rdefines.h>
#include <rph_util.h>

//convert from C ListOfLists object to a list in R.
SEXP rph_listOfLists_to_SEXP(ListOfLists *lol) {
  int numprotect=0, listlen,lstSize, col, *intp, i, hasHeader=0;
  double *doublep;
  List *currlist;
  SEXPTYPE lstType;
  SEXP result, header=R_NilValue, tempvec=R_NilValue;

  listlen = lst_size(lol->lst);

  PROTECT(result = allocVector(VECSXP, listlen));
  numprotect++;

  if (lol->lstName != NULL) {
    for (i=0; i<listlen; i++) 
      if (lst_get_ptr(lol->lstName, i) !=NULL) {
	hasHeader=1;
	break;
      }
  }
  if (hasHeader) {
    PROTECT(header = allocVector(STRSXP, listlen));
    numprotect++;
  }
  for (col=0; col < listlen; col++) {
    checkInterrupt();
    lstType = lst_get_int(lol->lstType, col);
    if (lstType == LIST_LIST) {
      SET_VECTOR_ELT(result, col,
	             rph_listOfLists_to_SEXP((ListOfLists*)lst_get_ptr(lol->lst, col)));
    } else if (lstType == GFF_PTR_LIST) {
      SET_VECTOR_ELT(result, col,
		     rph_gff_new_extptr((GFF_Set*)lst_get_ptr(lol->lst, col)));
    } else if (lstType == MSA_PTR_LIST) {
      SET_VECTOR_ELT(result, col, rph_msa_new_extptr((MSA*)lst_get_ptr(lol->lst, col)));
    } else if (lstType == MS_PTR_LIST) {
      SET_VECTOR_ELT(result, col, rph_ms_new_extptr((MS*)lst_get_ptr(lol->lst, col)));
    } else {
      currlist = (List*)lst_get_ptr(lol->lst, col);
      lstSize = lst_size(currlist);
      if (lstType == INT_LIST) {
	SET_VECTOR_ELT(result, col, tempvec = allocVector(INTSXP, lstSize));
	intp = INTEGER_POINTER(tempvec);
	for (i=0; i<lstSize; i++) {
	  intp[i] = lst_get_int(currlist, i);
	}
      } else if (lstType == DBL_LIST) {
	SET_VECTOR_ELT(result, col, tempvec = allocVector(REALSXP, lstSize));
	doublep = NUMERIC_POINTER(tempvec);
	for (i=0; i < lstSize; i++) {
	  doublep[i] = lst_get_dbl(currlist, i);
	}
      } else if (lstType == CHAR_LIST) {
	SET_VECTOR_ELT(result, col, tempvec = allocVector(STRSXP, lstSize));
	for (i=0; i < lstSize; i++) {
	  SET_STRING_ELT(tempvec, i, mkChar((char*)lst_get_ptr(currlist, i)));
	}
      } else {
	die("rph_listOfList_to_R unknown type %i\n", lstType);
      }
    }

    if (lol->lstName != NULL && lst_size(lol->lstName) > col) {
      char *name = lst_get_ptr(lol->lstName, col);
      if ((void*)name != NULL)
	SET_STRING_ELT(header, col, mkChar((char*)lst_get_ptr(lol->lstName, col)));
    }
  }

  if (lol->lstName != NULL)
    SET_NAMES(result, header);
  
  if (lol->class != NULL) {
    SEXP attrName, attrVal;
    PROTECT(attrName = allocVector(STRSXP, 1));
    SET_STRING_ELT(attrName, 0, mkChar("class"));
    PROTECT(attrVal = allocVector(STRSXP, 1));
    SET_STRING_ELT(attrVal, 0, mkChar(lol->class));
    SET_ATTR(result, attrName, attrVal);
    numprotect += 2;
  }

  if (numprotect > 0) UNPROTECT(numprotect);
  return result;
}



