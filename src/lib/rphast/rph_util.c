/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <rph_util.h>

SEXP rph_new_mem_handler() {
  phast_new_mem_handler();
  return R_NilValue;
}


SEXP rph_free_all() {
  phast_free_all();
  return R_NilValue;
}

//return a phast Vector from an R vector
Vector *rph_get_vector(SEXP doubleP) {
  int dim, i;
  Vector *rv;
  double *p;
  PROTECT(doubleP = AS_NUMERIC(doubleP));
  p = NUMERIC_POINTER(doubleP);
  dim = LENGTH(doubleP);
  rv = vec_new(dim);
  for (i=0; i<dim; i++)
    vec_set(rv, i, p[i]);
  UNPROTECT(1);
  return rv;
}


// returns a phast Matrix from an R matrix.  Assumes matrix is SQUARE!
Matrix *rph_get_matrix(SEXP matP) {
  double *doubleP;
  Matrix *m;
  int dim, pos, i, j;
  PROTECT(matP = AS_NUMERIC(matP));
  doubleP = NUMERIC_POINTER(matP);
  dim = (int)sqrt(LENGTH(matP));
  m = mat_new(dim, dim);
  pos=0;
  for (i=0; i<dim; i++)
    for (j=0; j<dim; j++)
      mat_set(m, j, i, doubleP[pos++]);
  UNPROTECT(1);
  return m;
}


SEXP rph_lst_len(SEXP listP) {
  List* l= (List*)EXTPTR_PTR(listP);
  SEXP rv;
  PROTECT(rv = allocVector(INTSXP, 1));
  INTEGER(rv)[0] = (int)lst_size(l);
  UNPROTECT(1);
  return rv;
}


void rph_lst_free(SEXP listP) {
  List *l;
  l = (List*)EXTPTR_PTR(listP);
  lst_free(l);
}


SEXP rph_lst_new_extptr(List *l) {
  SEXP result;
  PROTECT(result=R_MakeExternalPtr((void*)l, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, rph_lst_free, 1);
  lst_protect(l);
  UNPROTECT(1);
  return result;
}


