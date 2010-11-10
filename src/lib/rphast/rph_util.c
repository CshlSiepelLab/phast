#include "misc.h"
#include "stringsplus.h"
#include <Rdefines.h>
#undef Matrix
#include <matrix.h>

static void **mem_list=NULL;
static int mem_list_len=0;
static int mem_list_alloc_len=0;
static void ***static_mem_list=NULL;
static int static_mem_list_len=0;
static int static_mem_list_alloc_len=0;
static int *mem_available_list=NULL;
static int mem_available_list_len = 0;
static int mem_available_alloc_len = 0;
#define MEM_LIST_START_SIZE 100000
#define MEM_LIST_INCREASE_SIZE 1000000

/*
  RPHAST memory handler.  Basic idea- all C memory allocations must be done
  via smalloc or srealloc (careful to use copy_charstr rather than strdup,
  which uses malloc).  All frees have to be done with sfree.  Every time 
  memory is allocated it gets added to mem_list.  Each allocated memory is 
  bigger than the requested size by sizeof(int), and smalloc actually returns 
  ((int*)ptr)[1].  At position 0 is stored the index in mem_list where we 
  have a pointer to this memory.

  Many rphast functions use on.exit(freeall.rphast).  This invokes
  rph_free_all which frees all memory in mem_list, and also resets all 
  the mem_list, static_mem_list, mem_available_list, etc.

  If we want to keep objects past a freeall.rphast call, it needs to be
  protected.  This is simply done by setting mem_list[i] to NULL for
  this piece of memory.  Anything that is protected should have a finalizer
  arranged.  There are functions rph_msa_protect, rph_gff_protect, 
  rph_lst_protect, etc.  Note that rph_lst_protect doesn't protect the 
  objects in the list.

  sfree sets mem_list[(int*)ptr[-1]]=NULL, and frees the memory.

  mem_available_list is an array of indices in mem_list whose values no
  longer need to be stored (either because they have been protected or freed).
  When new memory is allocated it recycles indices from this list, if available,
  to prevent mem_list from getting too long (there are a surprising number
  of allocations performed when, for example, reading in an entire chromsome
  alignment).  Without this, the memory handler can occasionally take much more
  memory than the rest of the program!

  There is a complication having to do with static variables within some
  phast functions.  One option would be to protect them all.  But that would
  involve a lot of ugly protect statements in the code, and it's not clear
  that it is appropriate for these values to always be stored between various
  rphast calls, even if it is appropriate to store them throughout a run
  of command-line phast.  So instead there is set_static_var, which stores
  a pointer to the static variable in static_mem_list (which is the same 
  length as mem_list; all non-static variables have NULL in this array).  
  Then, when rph_free_all is called, these pointers get set to NULL.  
  Functions which use static variables need to check whether the this value 
  is NULL to see if it needs to be re-computed.  It is only necessary to use
  set_static_var on the static variables whose values are checked to see if
  they need to be reset.

  Another complication arises in R code.  There are many places in the R code
  where we use as.pointer.obj to convert an object to a pointer, and then
  send that pointer to a C function via .Call.  This is perfectly safe for msa
  and feat objects, since external pointers to these objects are always
  protected and have finalizers attached.  But tm and hmm objects do not
  get automatically protected (at least for now), since as.pointer.tm and
  as.pointer.hmm are internal functions (not visible to user) and we generally
  do not need to keep these pointers around for very long.  So, whenever
  as.pointer.hmm or as.pointer.tm are called, we have to be careful not to
  call any R functions which will invoke freeall.rphast until we are done
  using these pointers.  (Almost all rphast functions besides a few internal
  ones do invoke freeall.rphast.)  Or we have to protect them.  It may be a 
  good idea at some point to automatically protect them, if this becomes a 
  problem, but usually it is easy to create these pointers at the last moment.

  Finally, one other complication arises with protected objects.  If a C 
  function could potentially change the value of an object, and particularly 
  if it can allocate some new memory and have the object point to it- 
  rph_mem_protect needs to be called on the object again.  It is perfectly 
  safe (though perhaps inefficient) to call rph_mem_protect multiple times on 
  the same object.  The efficiency is a factor of the number of pointers in 
  the object.  MSA objects stored as SS can be particularly inefficienct, since
  each tuple is allocated individually.
 */


void rph_make_mem_list() {
  mem_list = malloc(MEM_LIST_START_SIZE * sizeof(void*));
  mem_list_len = 0;
  mem_list_alloc_len = MEM_LIST_START_SIZE;
}

void rph_realloc_mem_list() {
  mem_list_alloc_len += MEM_LIST_INCREASE_SIZE;
  mem_list = realloc(mem_list, mem_list_alloc_len*sizeof(void*));
}


void rph_add_to_mem_list(void *ptr) {
  int *iptr, idx;
  if (mem_list == NULL)
    rph_make_mem_list();
  if (mem_available_list_len != 0) {
    idx = mem_available_list[mem_available_list_len-1];
    mem_available_list_len--;
  } else {
    idx = mem_list_len++;
    if (idx == mem_list_alloc_len) 
      rph_realloc_mem_list();
  }
  mem_list[idx] = ptr;
  iptr = (int*)ptr;
  iptr[0] = idx;
}


void rph_protect_mem(void *ptr) {
  int *iptr = (int*)ptr;
  if (iptr[-1] >= 0) {
    mem_list[iptr[-1]] = NULL;
    iptr[-1] = -1;
  }
}

SEXP rph_free_all() {
  int i;
  if (mem_list == NULL) return R_NilValue;
  for (i=0; i < mem_list_len; i++) {
    if (mem_list[i] != NULL) 
      free(mem_list[i]);
  }
  free(mem_list);
  for (i=0; i < static_mem_list_len; i++)
    *(static_mem_list[i]) = NULL;
  free(static_mem_list);
  mem_list = NULL;
  static_mem_list = NULL;
  mem_list_alloc_len = mem_list_len = 0;
  static_mem_list_len = static_mem_list_alloc_len = 0;
  if (mem_available_list != NULL) {
    free(mem_available_list);
    mem_available_list = NULL;
  }
  mem_available_list_len = mem_available_alloc_len = 0;
  return R_NilValue;
}
  
void *smalloc(size_t size) {
  int *retval = (int*)malloc(size + sizeof(int));
  if (retval == NULL)
    die("ERROR: out of memory\n");
  rph_add_to_mem_list(retval);
  return (void*)(retval + 1);
}


void set_static_var(void **ptr) {
  if (static_mem_list_len == static_mem_list_alloc_len) {
    if (static_mem_list_alloc_len == 0) {
      static_mem_list_alloc_len = MEM_LIST_START_SIZE;
      static_mem_list = malloc(static_mem_list_alloc_len*sizeof(void**));
    }
    else {
      static_mem_list_alloc_len += MEM_LIST_INCREASE_SIZE;
      static_mem_list = realloc(static_mem_list,
				static_mem_list_alloc_len*sizeof(void*));
    }
  }
  static_mem_list[static_mem_list_len++] = ptr;
}

  
void *srealloc(void *ptr, size_t size) {
  int *iptr, *newptr, oldval;
  if (ptr == NULL) return smalloc(size);
  iptr = ((int*)ptr) - 1;
  if (size == 0) {
    if (iptr[0] >= 0 && mem_list[iptr[0]] != NULL) 
      sfree(ptr);
    return NULL;
  }
  oldval = iptr[0];
  newptr = (int*)realloc(iptr, size + sizeof(int));
  if (oldval >= 0) mem_list[oldval] = newptr;
  return (void*)(newptr + 1);
}

void sfree(void *ptr) {
  int *iptr;
  if (ptr == NULL) return;
  iptr = ((int*)ptr - 1);
  if (iptr[0] >= 0) {
    if (iptr[0] == mem_list_len-1)
      mem_list_len--;
    else {
      if (mem_available_list == NULL) {
	mem_available_list = malloc(MEM_LIST_START_SIZE*sizeof(int));
      } else if (mem_available_list_len == mem_available_alloc_len) {
	mem_available_alloc_len += MEM_LIST_INCREASE_SIZE;
	mem_available_list = realloc(mem_available_list, mem_available_alloc_len*sizeof(int));
      }
      mem_list[iptr[0]] = NULL;
      mem_available_list[mem_available_list_len++] = iptr[0];
    }
  }
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
