/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <memory_handler.h>

typedef struct mem_list_type MemList;

struct protected_object_struct {
  void *object;
  void (*function)(void*);
};

struct mem_list_type {
  void **mem_list;
  int mem_list_len;
  int mem_list_alloc_len;
  void ***static_mem_list;
  int static_mem_list_len;
  int static_mem_list_alloc_len;
  void **mem_available_list;
  int mem_available_list_len;
  int mem_available_alloc_len;
  struct protected_object_struct *protected_objects;
  int num_protected_objects;
  int protected_objects_alloc_len;
  FILE **open_files;
  int num_open_files;
  int open_files_alloc_len;
};

static MemList *memlist=NULL;
#ifdef USE_PHAST_MEMORY_HANDLER
static MemList *big_memlist = NULL;
static int num_memlist=0;
#endif

#define FILE_LIST_SIZE 100
#define MEM_LIST_START_SIZE 100000
#define MEM_LIST_INCREASE_SIZE 1000000

#ifdef RPHAST
#undef malloc 
#define malloc(x) (void*)Calloc((x),char)
#undef realloc
#define realloc(x,n) (void*)Realloc((x),(n),char)
#undef free
#define free(x) Free((x))
#endif


/*
  RPHAST memory handler.  Basic idea (see next section for usage rules):

  1. All C memory allocation  happens via smalloc.  smalloc(n) allocates
     an object x of size n+sizeof(void*), and returns ((void**)x)[1].  At
     x[0] is stored a void* pointer to mem_list[i], where mem_list is a 
     void** and mem_list[i] points to x.

  2. All memory frees should happen via sfree.  This frees x and sets
     mem_list[i] to NULL, and adds mem_list[i] to mem_available_list.

  3. Whenever R functions call C functions via .Call.rphast, phast_free_all
     is automatically called once the C function has returned.  phast_free_all
     frees all memory pointed to by mem_list.

  4. Any memory associated with external pointer objects that should not be 
     freed by phast_free_all needs to be "protected".  protection is done by 
     setting mem_list[i] and x[0] to NULL without freeing x.  All protected 
     memory should have finalizers arranged to prevent memory leaks.

  5. Any C function which takes an external pointer object and potentially
     assigns new memory to any part of that object needs to re-protect the
     object.  Protection needs to be arranged before the change is made,
     because the code could be interrupted mid-change.  However the
     protection needs to occur after the new memory is allocated.  Therefore
     objects can be registered for protection via functions like 
     msa_register_protect, tm_register_protect, etc.  Then
     phast_free_all first protects registered objects before freeing all 
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
     phast_free_all is called, x is freed and reset to NULL.  Functions which
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
     done automatically by functions phast_msa_new_extptr, phast_tm_new_extptr,
     phast_gff_new_extptr, phast_hmm_new_extptr, etc.

  3. Any C function which takes an external pointer object and potentially
     allocates new memory which becomes associated with that object must
     arrange for the object to be re-protected.  This can be done with
     functions msa_register_protect, tm_register_protect, etc,
     which should be called *before* any change is made to the objet.

  4. Any C function which uses a static variable x should call 
     set_static_var(&x).  This causes x to be re-set to NULL once it is 
     freed.  Then the function can check if x is NULL and recompute the
     value as necessary.
 */


void phast_init_memlist() {
  memlist->mem_list=NULL;
  memlist->mem_list_len=0;
  memlist->mem_list_alloc_len=0;
  memlist->static_mem_list=NULL;
  memlist->static_mem_list_len=0;
  memlist->static_mem_list_alloc_len=0;
  memlist->mem_available_list=NULL;
  memlist->mem_available_list_len = 0;
  memlist->mem_available_alloc_len = 0;
  memlist->protected_objects = NULL;
  memlist->num_protected_objects=0;
  memlist->protected_objects_alloc_len=0;
  memlist->open_files = NULL;
  memlist->num_open_files = 0;
  memlist->open_files_alloc_len = 0;
}

void phast_new_mem_handler() {
#ifdef USE_PHAST_MEMORY_HANDLER
  if (big_memlist == NULL) 
    big_memlist = malloc((num_memlist+1)*sizeof(MemList));
  else big_memlist = realloc(big_memlist, (num_memlist+1)*sizeof(MemList));
  memlist = &big_memlist[num_memlist];
  phast_init_memlist();
  num_memlist++;
#else
  phast_warning("warning: phast_new_memory_handler: memory handler not turned on.  To use, recompile with -DUSE_PHAST_MEMORY_HANDLER.\n");
#endif
}


int phast_num_mem_handlers() {
#ifdef USE_PHAST_MEMORY_HANDLER
  return num_memlist;
#else
  return -1;
#endif
}

void phast_make_mem_list() {
#ifdef USE_PHAST_MEMORY_HANDLER
  memlist->mem_list = malloc(MEM_LIST_START_SIZE * sizeof(void*));
  memlist->mem_list_alloc_len = MEM_LIST_START_SIZE;
#endif
}

void phast_realloc_mem_list() {
  int i;
  void *old_mem_list, **ptr;
  old_mem_list= memlist->mem_list;
  if (memlist->mem_available_list_len != 0) die("error");  //we shouldn't be reallocating if there are available spots, so we don't have to worry about these pointers being destroyed
  memlist->mem_list_alloc_len += MEM_LIST_INCREASE_SIZE;
  memlist->mem_list = realloc(memlist->mem_list, memlist->mem_list_alloc_len*sizeof(void*));
  if (memlist->mem_list != old_mem_list) {
    for (i=0; i < memlist->mem_list_len; i++) {
      if (memlist->mem_list[i] != NULL) {
	ptr = (void**)memlist->mem_list[i];
	ptr[0] = &memlist->mem_list[i];
      }
    }
  }
}


void phast_add_to_mem_list(void **ptr) {
  int idx;
  if (memlist->mem_list == NULL)
    phast_make_mem_list();
  if (memlist->mem_available_list_len != 0) {
    *((void**)memlist->mem_available_list[memlist->mem_available_list_len-1]) = (void*)ptr;
    ptr[0] = memlist->mem_available_list[memlist->mem_available_list_len-1];
    memlist->mem_available_list_len--;
  } else {
    if (memlist->mem_list_len == memlist->mem_list_alloc_len)
      phast_realloc_mem_list();
    idx = memlist->mem_list_len++;
    memlist->mem_list[idx] = (void*)ptr;
    ptr[0] = &memlist->mem_list[idx];
  }
}

void phast_add_to_mem_available_list(void *val) {
  if (memlist->mem_available_list == NULL) {
    memlist->mem_available_list = malloc(MEM_LIST_START_SIZE*sizeof(void*));
    memlist->mem_available_alloc_len = MEM_LIST_START_SIZE;
  } else if (memlist->mem_available_list_len == memlist->mem_available_alloc_len) {
    memlist->mem_available_alloc_len += MEM_LIST_INCREASE_SIZE;
    memlist->mem_available_list = realloc(memlist->mem_available_list, memlist->mem_available_alloc_len*sizeof(void*));
  }
  memlist->mem_available_list[memlist->mem_available_list_len++] = val;
}



void phast_mem_protect(void *ptr0) {
#ifdef USE_PHAST_MEMORY_HANDLER
  void **ptr = ((void**)ptr0)-1;
  if (ptr[0] != NULL) {
    phast_add_to_mem_available_list(ptr[0]);
    *(void**)ptr[0] = NULL;
    ptr[0] = NULL;
  }
#endif
}


void phast_register_protected_object(void *ptr, void (*function)(void*)) {
#ifdef USE_PHAST_MEMORY_HANDLER
  int idx;
  if (memlist->protected_objects == NULL) {
    memlist->protected_objects = malloc(100*sizeof(struct protected_object_struct));
    memlist->protected_objects_alloc_len = 100;
  }
  idx = memlist->num_protected_objects++;
  if (idx == memlist->protected_objects_alloc_len) {
    memlist->protected_objects_alloc_len += 1000;
    memlist->protected_objects = realloc(memlist->protected_objects, 
					memlist->protected_objects_alloc_len*
					sizeof(struct protected_object_struct));
  }
  memlist->protected_objects[idx].object = ptr;
  memlist->protected_objects[idx].function = function;
#endif
}

//this is inefficient but the protected_objects list is generally length
//1 or 2, and currently should not exceed 4.
void phast_unregister_protected(void *ptr) {
  int i;
  if (memlist == NULL) return;
  for (i=0; i < memlist->num_protected_objects; i++) {
    if (memlist->protected_objects[i].object == ptr) {
      memlist->protected_objects[i].object = NULL;
      break;
    }
  }
}


void phast_free_all() {
#ifdef USE_PHAST_MEMORY_HANDLER
  int i;
  if (memlist->protected_objects != NULL) {
    for (i=0; i < memlist->num_protected_objects; i++) {
      if (memlist->protected_objects[i].object != NULL) 
	(*memlist->protected_objects[i].function)(memlist->protected_objects[i].object);
    }
    if (memlist->protected_objects != NULL) {
      free(memlist->protected_objects);
      memlist->protected_objects = NULL;
      memlist->num_protected_objects = 0;
      memlist->protected_objects_alloc_len = 0;
    }
  }
  if (memlist->mem_list != NULL) {
    for (i=0; i < memlist->mem_list_len; i++) {
      if (memlist->mem_list[i] != NULL) 
	free(memlist->mem_list[i]);
    }
    free(memlist->mem_list);
    for (i=0; i < memlist->static_mem_list_len; i++)
      *(memlist->static_mem_list[i]) = NULL;
    free(memlist->static_mem_list);
    memlist->mem_list = NULL;
    memlist->static_mem_list = NULL;
    memlist->mem_list_alloc_len = memlist->mem_list_len = 0;
    memlist->static_mem_list_len = memlist->static_mem_list_alloc_len = 0;
    if (memlist->mem_available_list != NULL) {
      free(memlist->mem_available_list);
      memlist->mem_available_list = NULL;
    }
    memlist->mem_available_list_len = memlist->mem_available_alloc_len = 0;
    for (i=0; i < memlist->num_open_files; i++) {
      if (memlist->open_files[i] != NULL) 
	fclose(memlist->open_files[i]);
    }
    if (memlist->open_files != NULL)
      free(memlist->open_files);
    memlist->num_open_files = 0;
    memlist->open_files_alloc_len = 0;
  }
  num_memlist--;
  if (num_memlist > 0)
    memlist = &big_memlist[num_memlist-1];
  else {
    memlist = NULL;
    free(big_memlist);
    big_memlist = NULL;
  }
#endif
}


void set_static_var(void **ptr) {
#ifdef USE_PHAST_MEMORY_HANDLER
  if (memlist->static_mem_list_len == memlist->static_mem_list_alloc_len) {
    if (memlist->static_mem_list_alloc_len == 0) {
      memlist->static_mem_list_alloc_len = MEM_LIST_START_SIZE;
      memlist->static_mem_list = malloc(memlist->static_mem_list_alloc_len*sizeof(void**));
    }
    else {
      memlist->static_mem_list_alloc_len += MEM_LIST_INCREASE_SIZE;
      memlist->static_mem_list = realloc(memlist->static_mem_list,
					 memlist->static_mem_list_alloc_len*sizeof(void*));
    }
  }
  memlist->static_mem_list[memlist->static_mem_list_len++] = ptr;
#endif
}



/* This is not very efficient but there should never be a huge number
   of open files
 */
void register_open_file(FILE *F) {
#ifdef USE_PHAST_MEMORY_HANDLER
  int i;
  for (i=0; i < memlist->num_open_files; i++)
    if (memlist->open_files[i]==NULL) {
      memlist->open_files[i] = F;
      return;
    }
  if (memlist->open_files_alloc_len <= memlist->num_open_files) {
    if (memlist->open_files == NULL) 
      memlist->open_files = malloc(FILE_LIST_SIZE*sizeof(FILE*));
    else memlist->open_files = realloc(memlist->open_files, 
				       (memlist->open_files_alloc_len + FILE_LIST_SIZE)*sizeof(FILE*));
  }
  memlist->open_files[memlist->num_open_files++] = F;
#endif
}

void unregister_open_file(FILE *F) {
#ifdef USE_PHAST_MEMORY_HANDLER
  int i;
  for (i=0; i < memlist->num_open_files; i++)
    if (memlist->open_files[i] == F) {
      memlist->open_files[i] = NULL;
      return;
    }
  die("Memory handler error: could not un-register file\n");
#endif
}

/* without mem handler versions of smalloc, srealloc, sfree defined in misc.h */
#ifdef USE_PHAST_MEMORY_HANDLER  
void *smalloc(size_t size) {
  void **retval = (void**)malloc(size + sizeof(void*));
  if (retval == NULL)
    die("ERROR: out of memory\n");
  phast_add_to_mem_list(retval);
  return (void*)(retval+1);
}

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

void sfree(void *ptr0) {
  void **ptr;
  if (ptr0 == NULL) return;
  ptr = (void**)ptr0-1;
  if (ptr[0] != NULL) {
    phast_add_to_mem_available_list(ptr[0]);
    *(void**)ptr[0] = NULL;
  }
#ifdef RPHAST
  Free(ptr);
#else
  free((void*)ptr);
#endif
}
#endif


void lst_protect(List *l) {
  if (l == NULL) return;
  phast_mem_protect(l);
  phast_mem_protect(l->array);
}


void str_protect(String *s) {
  if (s == NULL) return;
  phast_mem_protect(s);
  if (s->chars != NULL) phast_mem_protect(s->chars);
}

void vec_protect(Vector *v) {
  if (v == NULL) return;
  phast_mem_protect(v);
  if (v->data != NULL) phast_mem_protect(v->data);
}


void zvec_protect(Zvector *v) {
  if (v == NULL) return;
  phast_mem_protect(v);
  if (v->data != NULL) phast_mem_protect(v->data);
}


void mat_protect(Matrix *m) {
  int i;
  if (m == NULL) return;
  phast_mem_protect(m);
  if (m->data != NULL) {
    for (i=0; i < m->nrows; i++)
      phast_mem_protect(m->data[i]);
    phast_mem_protect(m->data);
  }
}


void zmat_protect(Zmatrix *m) {
  int i;
  if (m == NULL) return;
  phast_mem_protect(m);
  if (m->data != NULL) {
    for (i=0; i < m->nrows; i++) 
      phast_mem_protect(m->data[i]);
    phast_mem_protect(m->data);
  }
}


void mm_protect(MarkovMatrix *mm) {
  if (mm == NULL) return;
  phast_mem_protect(mm);
  mat_protect(mm->matrix);
  zmat_protect(mm->evec_matrix_z);
  zmat_protect(mm->evec_matrix_inv_z);
  zvec_protect(mm->evals_z);
  mat_protect(mm->evec_matrix_r);
  mat_protect(mm->evec_matrix_inv_r);
  vec_protect(mm->evals_r);
  phast_mem_protect(mm->states);
}


void gp_protect(GapPatternMap *gpm) {
  int i;
  if (gpm == NULL) return;
  phast_mem_protect(gpm);
  phast_mem_protect(gpm->gapcat_to_cat);
  phast_mem_protect(gpm->gapcat_to_pattern);
  if (gpm->cat_x_pattern_to_gapcat != NULL) {
    phast_mem_protect(gpm->cat_x_pattern_to_gapcat);
    for (i=0; i < gpm->ncats; i++) {
      if (gpm->cat_x_pattern_to_gapcat[i] != NULL)
	phast_mem_protect(gpm->cat_x_pattern_to_gapcat[i]);
    }
  }
  if (gpm->node_to_branch != NULL)
    phast_mem_protect(gpm->node_to_branch);
  if (gpm->pattern_to_node != NULL)
    phast_mem_protect(gpm->pattern_to_node);
  if (gpm->pattern != NULL) {
    phast_mem_protect(gpm->pattern);
    for (i=0; i < gpm->ngap_patterns; i++)
      phast_mem_protect(gpm->pattern[i]);
  }
}


void tm_rmp_protect(TreeModel *tm) {
  int nparams = tm_get_nparams(tm);
  int i;
  for (i=0; i < nparams; i++) {
    if (tm->rate_matrix_param_row[i] != NULL) {
      lst_protect(tm->rate_matrix_param_row[i]);
      lst_protect(tm->rate_matrix_param_col[i]);
    }
  }
  phast_mem_protect(tm->rate_matrix_param_row);
  phast_mem_protect(tm->rate_matrix_param_col);
}

void tm_altmod_protect(AltSubstMod *am) {
  int i;
  phast_mem_protect(am);
  if (am->backgd_freqs != NULL)
    vec_protect(am->backgd_freqs);
  if (am->rate_matrix != NULL)
    mm_protect(am->rate_matrix);
  if (am->param_list != NULL) {
    lst_protect(am->param_list);
    for (i=0; i < lst_size(am->param_list); i++)
      str_protect(lst_get_ptr(am->param_list, i));
  }
  str_protect(am->defString);
  if (am->noopt_arg != NULL)
    str_protect(am->noopt_arg);
}


void tm_protect(TreeModel *tm) {
  int i, j;
  phast_mem_protect(tm);
  tree_protect(tm->tree);
  vec_protect(tm->backgd_freqs);
  mm_protect(tm->rate_matrix);
  if (tm->msa_seq_idx != NULL) phast_mem_protect(tm->msa_seq_idx);
  if (tm->P != NULL) {
    for (i=0; i < tm->tree->nnodes; i++) {
      for (j=0; j < tm->nratecats; j++)
	mm_protect(tm->P[i][j]);
      phast_mem_protect(tm->P[i]);
    }
    phast_mem_protect(tm->P);
  }
  if (tm->rK != NULL) phast_mem_protect(tm->rK);
  if (tm->freqK != NULL) phast_mem_protect(tm->freqK);
  if (tm->rate_matrix_param_row != NULL) 
    tm_rmp_protect(tm);
  if (tm->in_subtree != NULL) phast_mem_protect(tm->in_subtree);
  if (tm->ignore_branch != NULL) phast_mem_protect(tm->ignore_branch);
  if (tm->alt_subst_mods != NULL) {
    lst_protect(tm->alt_subst_mods);
    for (i=0; i < lst_size(tm->alt_subst_mods); i++)
      tm_altmod_protect(lst_get_ptr(tm->alt_subst_mods, i));
  }
  if (tm->alt_subst_mods_ptr != NULL) {
    for (i=0; i < tm->tree->nnodes; i++)
      phast_mem_protect(tm->alt_subst_mods_ptr[i]);
    phast_mem_protect(tm->alt_subst_mods_ptr);
  }
  if (tm->all_params != NULL) vec_protect(tm->all_params);
  if (tm->param_map != NULL) phast_mem_protect(tm->param_map);
  if (tm->bound_arg != NULL) {
    lst_protect(tm->bound_arg);
    for (i=0; i < lst_size(tm->bound_arg); i++)
      str_protect(lst_get_ptr(tm->bound_arg, i));
  }
  if (tm->noopt_arg != NULL)
    str_protect(tm->noopt_arg);
  if (tm->iupac_inv_map != NULL) {
    for (i=0; i < 256; i++)
      if (tm->iupac_inv_map[i] != NULL) phast_mem_protect(tm->iupac_inv_map[i]);
    phast_mem_protect(tm->iupac_inv_map);
  }
}

void tm_register_protect(TreeModel *tm) {
  phast_register_protected_object(tm, (void (*)(void*))tm_protect);
}


void gff_feat_protect(GFF_Feature *feat) {
  phast_mem_protect(feat);
  str_protect(feat->seqname);
  str_protect(feat->source);
  str_protect(feat->feature);
  str_protect(feat->attribute);
}


void gff_protect(GFF_Set *gff) {
  int i;
  phast_mem_protect(gff);
  if (gff->features != NULL) {
    lst_protect(gff->features);
    for (i=0; i < lst_size(gff->features); i++) {
      gff_feat_protect(lst_get_ptr(gff->features, i));
    }
  }
  str_protect(gff->gff_version);
  str_protect(gff->source);
  str_protect(gff->source_version);
  str_protect(gff->date);
  if (gff->groups != NULL) {
    lst_protect(gff->groups);
    for (i=0; i < lst_size(gff->groups); i++) {
      GFF_FeatureGroup *g = lst_get_ptr(gff->groups, i);
      phast_mem_protect(g);
      str_protect(g->name);
      lst_protect(g->features);
    }
  }
  str_protect(gff->group_tag);
}


void gff_register_protect(GFF_Set *gff) {
  phast_register_protected_object(gff, (void (*)(void *))gff_protect);
}


void cm_protect(CategoryMap *cm) {
  int i;
  if (cm == NULL) return;
  phast_mem_protect(cm);
  for (i=0; i <= cm->ncats; i++) {
    int len=0;
    if (cm->ranges[i] != NULL) {
      len = cm->ranges[i]->end_cat_no - cm->ranges[i]->start_cat_no;
      lst_protect(cm->ranges[i]->feature_types);
      phast_mem_protect(cm->ranges[i]);
      lst_protect(cm->conditioned_on[i]);
      i += len;
    }
  }
  phast_mem_protect(cm->ranges);
  phast_mem_protect(cm->conditioned_on);
  phast_mem_protect(cm->labelling_precedence);
  phast_mem_protect(cm->fill_precedence);
  if (cm->unspooler != NULL) 
    phast_mem_protect(cm->unspooler);
}


void cm_register_protect(CategoryMap *cm) {
  phast_register_protected_object(cm, (void (*)(void *))cm_protect);
}


void msa_protect_ss(MSA_SS *ss) {
  int i;
  phast_mem_protect(ss);
  if (ss->col_tuples != NULL) {
    phast_mem_protect(ss->col_tuples);
    for (i=0; i < ss->alloc_ntuples; i++) 
      if (ss->col_tuples[i] != NULL)
	phast_mem_protect(ss->col_tuples[i]);
  }
  if (ss->tuple_idx != NULL)
    phast_mem_protect(ss->tuple_idx);
  if (ss->counts != NULL)
    phast_mem_protect(ss->counts);
  if (ss->cat_counts != NULL) {
    phast_mem_protect(ss->cat_counts);
    for (i=0; i < ss->msa->ncats; i++)
      phast_mem_protect(ss->cat_counts[i]);
  }
}

void msa_protect(MSA *msa) {
  int i;
  phast_mem_protect(msa);
  if (msa->alphabet != NULL)
    phast_mem_protect(msa->alphabet);
  if (msa->names != NULL) {
    phast_mem_protect(msa->names);
    for (i=0; i < msa->nseqs; i++)
      phast_mem_protect(msa->names[i]);
  }
  if (msa->seqs != NULL) {
    phast_mem_protect(msa->seqs);
    for (i=0; i < msa->nseqs; i++)
      phast_mem_protect(msa->seqs[i]);
  }
  if (msa->categories != NULL)
    phast_mem_protect(msa->categories);
  if (msa->ss != NULL) {
    if (msa != msa->ss->msa) {
      die("ss pointer mismatch\n");
    }
    msa_protect_ss(msa->ss);
  }
  if (msa->is_informative != NULL)
    phast_mem_protect(msa->is_informative);
}


void ms_protect(MS *ms) {
  int i;
  phast_mem_protect(ms);

  if (ms->names != NULL) {
    phast_mem_protect(ms->names);
    for (i=0; i < ms->nseqs; i++)
      phast_mem_protect(ms->names[i]);
  }
  if (ms->seqs != NULL) {
    phast_mem_protect(ms->seqs);
    for (i=0; i < ms->nseqs; i++)
      phast_mem_protect(ms->seqs[i]);
  }
  if (ms->idx_offsets != NULL) {
    phast_mem_protect(ms->idx_offsets);
  }
  if (ms->alphabet != NULL) {
    phast_mem_protect(ms->alphabet);
  }

}

void msa_register_protect(MSA *msa) {
  phast_register_protected_object(msa, (void (*)(void*))msa_protect);
}

void ms_register_protect(MS *ms) {
  phast_register_protected_object(ms, (void (*)(void*))ms_protect);
}


void hmm_protect(HMM *hmm) {
  int i;
  if (hmm == NULL) return;
  phast_mem_protect(hmm);
  mm_protect(hmm->transition_matrix);
  mat_protect(hmm->transition_score_matrix);
  vec_protect(hmm->begin_transitions);
  vec_protect(hmm->end_transitions);
  vec_protect(hmm->eq_freqs);
  vec_protect(hmm->begin_transition_scores);
  vec_protect(hmm->end_transition_scores);
  for (i=0; i < hmm->nstates; i++) {
    lst_protect(hmm->predecessors[i]);
    lst_protect(hmm->successors[i]);
  }
  phast_mem_protect(hmm->predecessors);
  phast_mem_protect(hmm->successors);
  lst_protect(hmm->begin_successors);
  lst_protect(hmm->end_predecessors);
}


void hmm_register_protect(HMM *hmm) {
  phast_register_protected_object(hmm, (void (*)(void *))hmm_protect);
}


void phmm_protect(PhyloHmm *p) {
  int i;
  if (p == NULL) return;
  phast_mem_protect(p);
  cm_protect(p->cm);
  hmm_protect(p->hmm);
  hmm_protect(p->functional_hmm);
  hmm_protect(p->autocorr_hmm);
  if (p->mods != NULL) {
    phast_mem_protect(p->mods);
    for (i=0; i < p->nmods; i++)
      tm_protect(p->mods[i]);
  }
  gp_protect(p->gpm);
  if (p->state_to_mod != NULL)
    phast_mem_protect(p->state_to_mod);
  if (p->state_to_cat != NULL)
    phast_mem_protect(p->state_to_cat);
  if (p->reverse_compl != NULL)
    phast_mem_protect(p->reverse_compl);
  if (p->cat_to_states != NULL) {
    phast_mem_protect(p->cat_to_states);
    for (i=0; i <= p->cm->ncats; i++) 
      lst_protect(p->cat_to_states[i]);
  }
  if (p->state_to_pattern != NULL)
    phast_mem_protect(p->state_to_pattern);
  if (p->emissions != NULL) {
    for (i=0; i < p->hmm->nstates; i++)
      if (p->state_pos[p->state_to_mod[i]] == i ||
          p->state_neg[p->state_to_mod[i]] == i || 
          p->state_to_pattern[i] >= 0)
        phast_mem_protect(p->emissions[i]);
    phast_mem_protect(p->emissions);
    phast_mem_protect(p->state_pos);
    phast_mem_protect(p->state_neg);
  }
  if (p->forward != NULL) {
    for (i=0; i < p->hmm->nstates; i++) 
      phast_mem_protect(p->forward[i]);
    phast_mem_protect(p->forward);
  }
  //topology isn't freed by phmm_free so shouldn't be protected?
  //  if (p->topology != NULL)
  //    tree_protect(p->topology);

  if (p->alpha != NULL) 
    phast_mem_protect(p->alpha);
  if (p->beta != NULL)
    phast_mem_protect(p->beta);
  if (p->tau != NULL)
    phast_mem_protect(p->tau);
  if (p->T != NULL) {
    phast_mem_protect(p->T);
    phast_mem_protect(p->t);
    for (i=0; i < p->functional_hmm->nstates; i++) {
      phast_mem_protect(p->T[i]);
      phast_mem_protect(p->t[i]);
    }
  }
  if (p->em_data != NULL) {
    phast_mem_protect(p->em_data);
    mat_protect(p->em_data->H);
    //msa_protect(p->em_data->msa);  //assume this is protected elsewhere/
  }

}


void phmm_register_protect(PhyloHmm *phmm) {
  phast_register_protected_object(phmm, (void (*)(void *))phmm_protect);
}


void tree_protect(TreeNode *tr) {
  TreeNode *n;
  int i;
  if (tr == NULL) return;
  phast_mem_protect(tr);
  for (i=0; i < tr->nnodes; i++) {
    n = (TreeNode*)lst_get_ptr(tr->nodes, i);
    phast_mem_protect(n);
    if (n->nodes != NULL) 
      lst_protect(n->nodes);
    if (n->preorder != NULL)
      lst_protect(n->preorder);
    if (n->postorder != NULL)
      lst_protect(n->postorder);
    if (n->inorder != NULL)
      lst_protect(n->inorder);
    if (n->label != NULL)
      phast_mem_protect(n->label);
  }
}

