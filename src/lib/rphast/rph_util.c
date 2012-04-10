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

struct rph_likelihood_struct {
  SEXP functionCall;
  SEXP env;
};


double rph_likelihood_wrapper(Vector *params, void *data) {
  SEXP paramsP, funcVal, fcall, env, fn;
  double *paramsd, rv;
  int i;
  env = ((struct rph_likelihood_struct*)data)->env;
  fn = ((struct rph_likelihood_struct*)data)->functionCall;
  if (!isFunction(fn)) 
    die("rph_likelihood_wrapper: fn is not a function");
  if (!isEnvironment(env)) 
    die("rph_likelihood_wrapper: env is not an environment");
  PROTECT(paramsP = NEW_NUMERIC(params->size));
  paramsd = NUMERIC_POINTER(paramsP);
  for (i=0; i < params->size; i++)
    paramsd[i] = params->data[i];
  PROTECT(fcall = lang2(fn, paramsP));
  PROTECT(funcVal = eval(fcall, env));
  rv = NUMERIC_VALUE(funcVal);
  UNPROTECT(3);
  return -rv;
}

SEXP rph_opt_bfgs(SEXP likelihoodFunctionP, SEXP paramsP, 
		  SEXP lowerP, SEXP upperP, SEXP precisionP, 
		  SEXP logfileP, SEXP envP) {
  
  struct rph_likelihood_struct *data = smalloc(sizeof(struct rph_likelihood_struct));
  Vector *params, *lower=NULL, *upper=NULL;
  int i, numprotect=1, numeval;
  opt_precision_type precision;
  double retval, val;
  ListOfLists *result;
  FILE *logfile=NULL;
 
  if (!isFunction(likelihoodFunctionP)) 
    die("rph_opt_bfgs: likelihoodFunction is not funtion\n");
  PROTECT(paramsP = AS_NUMERIC(paramsP));
  params = vec_new_from_array(NUMERIC_POINTER(paramsP), LENGTH(paramsP));
  if (lowerP != R_NilValue) {
    PROTECT(lowerP = AS_NUMERIC(lowerP)); numprotect++;
    lower = vec_new_from_array(NUMERIC_POINTER(lowerP), LENGTH(lowerP));
    for (i=0; i < lower->size; i++) {
      val = vec_get(lower, i);
      if (!isfinite(val))
	vec_set(lower, i, INFTY * (signbit(val)!=0 ? -1 : 1.0));
    }
  }
  if (upperP != R_NilValue) {
    PROTECT(upperP = AS_NUMERIC(upperP)); numprotect++;
    upper = vec_new_from_array(NUMERIC_POINTER(upperP), LENGTH(upperP));
    for (i=0; i <  upper->size; i++) {
      val = vec_get(upper, i);
      if (!isfinite(val))
	vec_set(upper, i, INFTY * (signbit(val)!=0 ? -1 : 1.0));
    }
  }
  precision = get_precision(CHARACTER_VALUE(precisionP));
  if (precision == OPT_UNKNOWN_PREC) die("unknown precision");

  data->functionCall = likelihoodFunctionP;
  data->env = envP;
  if (logfileP != R_NilValue) 
    logfile = phast_fopen(CHARACTER_VALUE(logfileP), "a");

  opt_bfgs(rph_likelihood_wrapper, params, data, &retval, lower, 
	   upper, logfile, NULL, precision, NULL, &numeval);

  if (logfile != NULL)
    phast_fclose(logfile);

  //need to create a list with likelihood value and parameter estimates
  result = lol_new(3);
  lol_push_dbl(result, &retval, 1, "value");
  lol_push_dbl(result, params->data, params->size, "par");
  lol_push_int(result, &numeval, 1, "neval");
  UNPROTECT(numprotect);
  return rph_listOfLists_to_SEXP(result);
}

