#include <rph_util.h>

static void **mem_list=NULL;
static int mem_list_len=0;
static int mem_list_alloc_len=0;
static void ***static_mem_list=NULL;
static int static_mem_list_len=0;
static int static_mem_list_alloc_len=0;
static void **mem_available_list=NULL;
static int mem_available_list_len = 0;
static int mem_available_alloc_len = 0;

struct protected_object_struct {
  void *object;
  void (*function)(void*);
};
static struct protected_object_struct *protected_objects = NULL;
static int num_protected_objects=0;
static int protected_objects_alloc_len=0;


#define MEM_LIST_START_SIZE 100000
#define MEM_LIST_INCREASE_SIZE 1000000

#define USE_RPHAST_MEMORY_HANDLER

/*
  RPHAST memory handler.  Basic idea (see next section for usage rules):

  1. All C memory allocation  happens via smalloc.  smalloc(n) allocates
     an object x of size n+sizeof(void*), and returns ((void**)x)[1].  At
     x[0] is stored a void* pointer to mem_list[i], where mem_list is a 
     void ** and mem_list[i] points to x.

  2. All memory frees should happen via sfree.  This frees x and sets
     mem_list[i] to NULL, and adds mem_list[i] to mem_available_list.

  3. Whenever R functions call C functions via .Call.rphast, rph_free_all
     is automatically called once the C function has returned.  rph_free_all
     frees all memory pointed to by mem_list.

  4. Any memory that should not be freed by rph_free_all needs to
     be "protected".  protection is done by setting mem_list[i] and x[0] to
     NULL without freeing x.  All protected memory should have finalizers
     arranged to prevent memory leaks.

  5. Any C function which takes an external pointer object and potentially
     assigns new memory to any part of that object needs to re-protect the
     object.  Protection needs to be arranged before the change is made,
     because the code could be interrupted mid-change.  However the
     protection needs to occur after the new memory is allocated.  Therefore
     objects can be registered for protection via functions like 
     rph_msa_register_protect, rph_tm_register_protect, etc.  Then
     rph_free_all first protects registered objects before freeing all 
     remaining memory.

  6. mem_available_list is a void** with pointers to positions in 
     mem_list where the memory has been freed or protected.  Whenever
     mem_available_list has non-zero length, smalloc will set x[0] to 
     mem_available_list[mem_available_list_len--] rather than unnecessarily
     increasing the length of mem_list.  Without this, mem_list can quickly
     grow to an unreasonable size due to loops with repeated smalloc/sfrees.
   
  7. Static variables within C functions require special handling.  They 
     could be protected, but an easier option is to use set_static_var(&x)
     after smalloc.  This adds the address of x to static_mem_list.  When
     rph_free_all is called, x is freed and reset to NULL.  Functions which
     use static variables can then check if the static variable is NULL, it
     needs to be re-computed.

  In order for this to work, some rules need to be followed writing new
    code:

  1. All allocations should be through smalloc, and memory should be freed
     using sfree.  Careful not to use strdup which calls malloc, instead
     use copy_charstr.  If any memory is allocated/freed with smalloc/free
     or malloc/sfree, RPHAST will probably crash.

  2. When new external pointer objects are created in C and passed to R,
     they need to be protected and have a finalizer registered.  This is
     done automatically by functions rph_msa_new_extptr, rph_tm_new_extptr,
     rph_gff_new_extptr, rph_hmm_new_extptr, etc.

  3. Any C function which takes an external pointer object and potentially
     allocates new memory which becomes associated with that object must
     arrange for the object to be re-protected.  This can be done with
     functions rph_msa_register_protect, rph_tm_register_protect, etc,
     which should be called *before* any change is made to the objet.

  4. Any C function which uses a static variable x should call 
     set_static_var(&x).  This causes x to be re-set to NULL once it is 
     freed.  Then the function can check if x is NULL and recompute the
     value as necessary.
 */


#ifdef USE_RPHAST_MEMORY_HANDLER
void rph_make_mem_list() {
  mem_list = malloc(MEM_LIST_START_SIZE * sizeof(void*));
  mem_list_len = 0;
  mem_list_alloc_len = MEM_LIST_START_SIZE;
}

void rph_realloc_mem_list() {
  int i;
  void *old_mem_list = mem_list, **ptr;
  if (mem_available_list_len != 0) die("error");  //we shouldn't be reallocating if there are available spots, so we don't have to worry about these pointers being destroyed
  mem_list_alloc_len += MEM_LIST_INCREASE_SIZE;
  mem_list = realloc(mem_list, mem_list_alloc_len*sizeof(void*));
  if (mem_list != old_mem_list) {
    for (i=0; i < mem_list_len; i++) {
      if (mem_list[i] != NULL) {
	ptr = (void**)mem_list[i];
	ptr[0] = &mem_list[i];
      }
    }
  }
}


void rph_add_to_mem_list(void **ptr) {
  int idx;
  if (mem_list == NULL)
    rph_make_mem_list();
  if (mem_available_list_len != 0) {
    *((void**)mem_available_list[mem_available_list_len-1]) = (void*)ptr;
    ptr[0] = mem_available_list[mem_available_list_len-1];
    mem_available_list_len--;
  } else {
    if (mem_list_len == mem_list_alloc_len)
      rph_realloc_mem_list();
    idx = mem_list_len++;
    mem_list[idx] = (void*)ptr;
    ptr[0] = &mem_list[idx];
  }
}

void rph_add_to_mem_available_list(void *val) {
  if (mem_available_list == NULL) {
    mem_available_list = malloc(MEM_LIST_START_SIZE*sizeof(void*));
    mem_available_alloc_len = MEM_LIST_START_SIZE;
  } else if (mem_available_list_len == mem_available_alloc_len) {
    mem_available_alloc_len += MEM_LIST_INCREASE_SIZE;
    mem_available_list = realloc(mem_available_list, mem_available_alloc_len*sizeof(void*));
  }
  mem_available_list[mem_available_list_len++] = val;
}

#endif



void rph_mem_protect(void *ptr0) {
#ifdef USE_RPHAST_MEMORY_HANDLER
  void **ptr = ((void**)ptr0)-1;
  if (ptr[0] != NULL) {
    rph_add_to_mem_available_list(ptr[0]);
    *(void**)ptr[0] = NULL;
    ptr[0] = NULL;
  }
#endif
}


void rph_register_protected_object(void *ptr, void (*function)(void*)) {
  int idx;
  if (protected_objects == NULL) {
    protected_objects = malloc(100*sizeof(struct protected_object_struct));
    protected_objects_alloc_len = 100;
  }
  idx = num_protected_objects++;
  if (idx == protected_objects_alloc_len) {
    protected_objects_alloc_len += 1000;
    protected_objects = realloc(protected_objects, 
				protected_objects_alloc_len*
				sizeof(struct protected_object_struct));
  }
  protected_objects[idx].object = ptr;
  protected_objects[idx].function = function;
}

//this is inefficient but the protected_objects list is generally length
//1 or 2, and currently should not exceed 4.
void rph_unregister_protected(void *ptr) {
  int i;
  for (i=0; i < num_protected_objects; i++) {
    if (protected_objects[i].object == ptr) {
      protected_objects[i].object = NULL;
      break;
    }
  }
}


SEXP rph_free_all() {
  int i;
#ifdef USE_RPHAST_MEMORY_HANDLER
  if (protected_objects != NULL) {
    for (i=0; i < num_protected_objects; i++) {
      if (protected_objects[i].object != NULL) 
	(*protected_objects[i].function)(protected_objects[i].object);
    }
    if (protected_objects != NULL) {
      free(protected_objects);
      protected_objects = NULL;
      num_protected_objects = 0;
      protected_objects_alloc_len = 0;
    }
  }
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
#endif
  return R_NilValue;
}


#ifdef USE_RPHAST_MEMORY_HANDLER
void *smalloc(size_t size) {
  void **retval = (void**)malloc(size + sizeof(void*));
  if (retval == NULL)
    die("ERROR: out of memory\n");
  rph_add_to_mem_list(retval);
  return (void*)(retval+1);
}
#else
void *smalloc(size_t size) {
  void *retval = malloc(size);
  if (retval == NULL) die("ERROR: out of memory\n");
  return retval;
}
#endif

void set_static_var(void **ptr) {
#ifdef USE_RPHAST_MEMORY_HANDLER
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
#endif
}


#ifdef USE_RPHAST_MEMORY_HANDLER  
void *srealloc(void *ptr0, size_t size) {
  void *newptr, **ptr;
  if (ptr0 == NULL) return smalloc(size);
  ptr = (void**)ptr0-1;
  if (size == 0) {
    sfree(ptr0);
    return NULL;
  }
  newptr = realloc((void*)ptr, size + sizeof(void*));
  ptr = (void**)newptr;
  if (ptr[0] != NULL)
    *(void**)ptr[0] = newptr;
  return (void*)(ptr+1);
}
#else
void *srealloc(void *ptr, size_t size) {
  void *retval = realloc(ptr, size);
  if (retval == NULL && ptr != NULL && size != 0)
    die("FATAL ERROR: out of memory.\n");
  return retval;
}
#endif

#ifdef USE_RPHAST_MEMORY_HANDLER

void sfree(void *ptr0) {
  void **ptr;
  if (ptr0 == NULL) return;
  ptr = (void**)ptr0-1;
  if (ptr[0] != NULL) {
    rph_add_to_mem_available_list(ptr[0]);
    *(void**)ptr[0] = NULL;
  }
  free((void*)ptr);
}
#else
void sfree(void *ptr) {
  free(ptr);
}
#endif


void rph_lst_protect(List *l) {
  if (l == NULL) return;
  rph_mem_protect(l);
  rph_mem_protect(l->array);
}


void rph_str_protect(String *s) {
  if (s == NULL) return;
  rph_mem_protect(s);
  if (s->chars != NULL) rph_mem_protect(s->chars);
}

void rph_vec_protect(Vector *v) {
  if (v == NULL) return;
  rph_mem_protect(v);
  if (v->data != NULL) rph_mem_protect(v->data);
}


void rph_zvec_protect(Zvector *v) {
  if (v == NULL) return;
  rph_mem_protect(v);
  if (v->data != NULL) rph_mem_protect(v->data);
}


void rph_mat_protect(Matrix *m) {
  int i;
  if (m == NULL) return;
  rph_mem_protect(m);
  if (m->data != NULL) {
    for (i=0; i < m->nrows; i++)
      rph_mem_protect(m->data[i]);
    rph_mem_protect(m->data);
  }
}


void rph_zmat_protect(Zmatrix *m) {
  int i;
  if (m == NULL) return;
  rph_mem_protect(m);
  if (m->data != NULL) {
    for (i=0; i < m->nrows; i++) 
      rph_mem_protect(m->data[i]);
    rph_mem_protect(m->data);
  }
}


void rph_mm_protect(MarkovMatrix *mm) {
  if (mm == NULL) return;
  rph_mem_protect(mm);
  rph_mat_protect(mm->matrix);
  rph_zmat_protect(mm->evec_matrix_z);
  rph_zmat_protect(mm->evec_matrix_inv_z);
  rph_zvec_protect(mm->evals_z);
  rph_mat_protect(mm->evec_matrix_r);
  rph_mat_protect(mm->evec_matrix_inv_r);
  rph_vec_protect(mm->evals_r);
  rph_mem_protect(mm->states);
}


void rph_gp_protect(GapPatternMap *gpm) {
  int i;
  if (gpm == NULL) return;
  rph_mem_protect(gpm);
  rph_mem_protect(gpm->gapcat_to_cat);
  rph_mem_protect(gpm->gapcat_to_pattern);
  if (gpm->cat_x_pattern_to_gapcat != NULL) {
    rph_mem_protect(gpm->cat_x_pattern_to_gapcat);
    for (i=0; i < gpm->ncats; i++) {
      if (gpm->cat_x_pattern_to_gapcat[i] != NULL)
	rph_mem_protect(gpm->cat_x_pattern_to_gapcat[i]);
    }
  }
  if (gpm->node_to_branch != NULL)
    rph_mem_protect(gpm->node_to_branch);
  if (gpm->pattern_to_node != NULL)
    rph_mem_protect(gpm->pattern_to_node);
  if (gpm->pattern != NULL) {
    rph_mem_protect(gpm->pattern);
    for (i=0; i < gpm->ngap_patterns; i++)
      rph_mem_protect(gpm->pattern[i]);
  }
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
