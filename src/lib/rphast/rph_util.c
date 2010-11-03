#include "misc.h"
#include "stringsplus.h"
#include <Rdefines.h>
#undef Matrix
#include <matrix.h>

static void **mem_list=NULL;
static int *mem_do_free=NULL;
static int mem_list_len=0;
static int mem_list_alloc_len=0;
#define MEM_LIST_START_SIZE 10000

void rph_make_mem_list() {
  mem_list = malloc(MEM_LIST_START_SIZE * sizeof(void*));
  mem_do_free = malloc(MEM_LIST_START_SIZE * sizeof(int));
  mem_list_len = 0;
  mem_list_alloc_len = MEM_LIST_START_SIZE;
}

void rph_realloc_mem_list() {
  mem_list_alloc_len += MEM_LIST_START_SIZE;
  mem_list = realloc(mem_list, mem_list_alloc_len*sizeof(void*));
  mem_do_free = realloc(mem_do_free, mem_list_alloc_len * sizeof(int));
}

void rph_add_to_mem_list(void *ptr) {
  int *iptr;
  if (mem_list == NULL)
    rph_make_mem_list();
  else if (mem_list_len == mem_list_alloc_len)
    rph_realloc_mem_list();
  mem_list[mem_list_len] = ptr;
  mem_do_free[mem_list_len] = 1;
  iptr = (int*)ptr;
  iptr[0] = mem_list_len++;
}


void rph_protect_mem(void *ptr) {
  int *iptr = (int*)ptr;
  mem_do_free[iptr[-1]] = 0;
}

SEXP rph_free_all() {
  int i;
  for (i=0; i < mem_list_len; i++)
    if (mem_do_free[i]) 
      free(mem_list[i]);
  free(mem_list);
  free(mem_do_free);
  mem_list = NULL;
  mem_do_free = NULL;
  mem_list_alloc_len = mem_list_len = 0;
  return R_NilValue;
}
  
void *smalloc(size_t size) {
  int *retval = (int*)malloc(size + sizeof(int));
  if (retval == NULL)
    die("ERROR: out of memory\n");
  rph_add_to_mem_list(retval);
  return (void*)(retval + 1);
}

void *srealloc(void *ptr, size_t size) {
  int *iptr, *newptr, oldval;
  if (ptr == NULL) return smalloc(size);
  iptr = ((int*)ptr) - 1;
  if (size == 0) {
    if (mem_do_free[iptr[0]]) {
      free(mem_list[iptr[0]]);
      mem_do_free[iptr[0]] = 0;
    }
    return NULL;
  }
  oldval = iptr[0];
  newptr = (int*)realloc(iptr, size + sizeof(int));
  mem_list[oldval] = newptr;
  return (void*)(newptr + 1);
}

void sfree(void *ptr) {
  int *iptr;
  if (ptr == NULL) return;
  iptr = ((int*)ptr - 1);
  mem_do_free[iptr[0]] = 0;
  free(iptr);
}


void rph_lst_protect(List *l) {
  if (l == NULL) return;
  rph_protect_mem(l);
  rph_protect_mem(l->array);
}


void rph_string_protect(String *s) {
  if (s == NULL) return;
  rph_protect_mem(s);
  if (s->chars != NULL) rph_protect_mem(s->chars);
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
  rph_lst_protect(l);
  UNPROTECT(1);
  return result;
}
