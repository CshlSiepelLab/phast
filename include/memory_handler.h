/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** \file memory_handler.h
    Memory handler for phast, not used for command-line phast, but useful
    for cleaning memory when linking phast to other platforms such as 
    R or python.

    To use, the macro USE_PHAST_MEMORY_HANDLER needs to be defined.
    When the memory handler is in use, memory cannot be allocated until 
    phast_new_mem_handler()  has been called.  The memory handler keeps 
    track of all allocated memory and frees it all upon a call to 
    phast_free_all().  Once phast_free_all() is called, a new memory 
    handler needs to be created before any allocation can occur.
    Objects can be protected from being freed by phast_free_all() via 
    protect functions, in which case the user needs to arrange to free
    them manually.

    If phast_new_mem_handler() is called multiple times without calling
    phast_free_all(), this creates a stack of memory handlers.  When 
    memory is allocated or freed, or registered for protection, it is
    handled by the memory handler on the top of the stack (which is the 
    one most recently created).  phast_free_all() will then free everything
    associated with the memory handler on top of the stack, and then pop
    that memory handler off the stack.  This is useful in RPHAST when 
    R can potentially call C functions, which can invoke R functions, which
    can invoke C functions, etc.
 */

#ifndef MEMORY_HANDLER_H
#define MEMORY_HANDLER_H

#include <stdio.h>
#include <stdlib.h>
#include "misc.h"
#include "lists.h"
#include "tree_model.h"
#include "gff.h"
#include "msa.h"
#include "ms.h"
#include "matrix.h"
#include "complex_matrix.h"
#include "vector.h"
#include "complex_vector.h"
#include "hmm.h"
#include "phylo_hmm.h"
#include "stringsplus.h"
#include "category_map.h"
#include "trees.h"
#include "sufficient_stats.h"

/** Start new memory handler and push it on top of memory-handler stack.
 */
void phast_new_mem_handler();

/** Get the size of the memory handler stack.

    @return the size of the memory handler stack, or -1 phast was compiled
    without memory handler support.  The size of the memory handler stack 
    is equal to the number of times phast_new_mem_handler() has been called 
    minus the number of times phast_free_all()  has been called.
 */
int phast_num_mem_handlers();

/** First protects any memory which has been registered for protection.
    Then free all memory which is not protected, that has been allocated
    since the last call to phast_new_mem_handler().  Pop the current
    memory handler off the stack.  If the stack is empty then
    phast_new_mem_handler() will have to be called before any memory is
    allocated again.
 */
void phast_free_all();

/** Register an object to be protected when phast_free_all() is called on 
    the current memory handler.  All memory associated with a protected
    object will not be freed by phast_free_all() and will need to be 
    manually freed.
    @param ptr The object to be protected
    @param function The function to call to invoke protection (usually "obj_protect", where obj denotes the type of object.  Examples are msa_protect, gff_protect, tm_protect, mat_protect, etc.
 */
void phast_register_protected_object(void *ptr, void (*function)(void*));

/** Un-register an object for protection.  This may need to be called before
    freeing an object which has been registered for protection.
 */
void phast_unregister_protected(void *ptr);

/** Allocate memory.
    @param size the size of memory (in bytes) to allocate.
    @return A pointer to newly allocated memory
    @note When memory handler is being used, the amount of memory which is
    allocated is increased by sizeof(void*), but the returned value is a
    pointer to the last size bytes of the allocated memory.
 */
void *smalloc(size_t size);

/**  Re-allocate memory.  Works the same as realloc, but is adapted for memory
     handler and reports an error on failure.
     @param ptr0 Pointer to memory to re-allocate.
     @param size New size.
 */
void *srealloc(void *ptr0, size_t size);

/** Free memory
    @param Pointer to object to free, as returned by smalloc or srealloc.
 */
void sfree(void *ptr0);

/** Register a static variable.  This is necessary when memory handler is
    in use.  This re-sets the value of the variable to NULL after it is
    freed by phast_free_all().
    @param ptr A pointer to a static variable
    @usage set_static_var((void**)&static_var);
 */
void set_static_var(void **ptr);

/**\name Object protection functions 
   \{ */

/** Protect memory from being freed by phast_free_all()
    @param A pointer to the memory to protect
 */
void phast_mem_protect(void *ptr0);

/** Protect a list from being freed by phast_free_all()
    @param l An list object
    @note If l is a pointer list, does not protect objects pointed to by the list.
 */
void lst_protect(List *l);

/** Protect a string object from being freed by phast_free_all()
    @param s A string object
 */
void str_protect(String *s);

/** Protect a vector from being freed by phast_free_all()
    @param v A vector
 */
void vec_protect(Vector *v);

/** Protect a complex vector from being freed by phast_free_all()
    @param v A complex vector
 */
void zvec_protect(Zvector *v);

/** Protect a matrix from being freed by phast_free_all()
    @param m A matrix
 */
void mat_protect(Matrix *m);

/** Protect a complex matrix from being freed by phast_free_all()
    @param m A complex matrix
 */
void zmat_protect(Zmatrix *m);

/** Protect a Markov matrix from being freed by phast_free_all()
    @param m A Markov matrix
 */

void mm_protect(MarkovMatrix *mm);

/** Protect a GapPatternmap object from being freed by phast_free_all()
    @param m A GapPatternMap object
 */
void gp_protect(GapPatternMap *gpm);

/** Protect the rate_matrix_param_row element of a tree model from
    being freed by phast_free_all()
    @param tm An tree model object
 */
void tm_rmp_protect(TreeModel *tm);

/** Protect an alternate substitution model object from being freed 
    by phast_free_all()
    @param am An alternate substition model object
 */
void tm_altmod_protect(AltSubstMod *am);

/** Protect a tree model object from being freed by phast_free_all()
    @param tm A tree model object */
void tm_protect(TreeModel *tm);

/** Protect a GFF_Feature object from being freed by phast_free_all()
    @param tm A GFF_Feature object */
void gff_feat_protect(GFF_Feature *feat);

/** Protect a GFF object from being freed by phast_free_all()
    @param tm A GFF object */
void gff_protect(GFF_Set *gff);

/** Protect a category map object from being freed by phast_free_all()
    @param tm A category map object */
void cm_protect(CategoryMap *cm);

/** Protect a sufficient statistics object from being freed by phast_free_all()
    @param tm A sufficient statistics object */
void msa_protect_ss(MSA_SS *ss);

/** Protect an MSA object from being freed by phast_free_all()
    @param tm An MSA object */
void msa_protect(MSA *msa);

/** Protect an HMM object from being freed by phast_free_all()
    @param tm An HMM object */
void hmm_protect(HMM *hmm);

/** Protect a phylo-HMM object from being freed by phast_free_all()
    @param tm A phylo-HMM object */
void phmm_protect(PhyloHmm *p);

/** Protect a tree object from being freed by phast_free_all()
    @param tm A tree object */
void tree_protect(TreeNode *tr);

/** \name Functions to register an object for protection
    \{ */

/** Register a tree model object for protection 
    @param A tree model object */
void tm_register_protect(TreeModel *tm);

/** Register a GFF object for protection 
    @param A GFF object */
void gff_register_protect(GFF_Set *gff);

/** Register a Category Map object for protection 
    @param A Category Map object */
void cm_register_protect(CategoryMap *cm);

/** Register an MSA object for protection 
    @param An MSA object */
void msa_register_protect(MSA *msa);

/** Register an MS object for protection
    @param An MS object */
void ms_register_protect(MS *ms);

/** Register an HMM object for protection 
    @param An HMM object */
void hmm_register_protect(HMM *hmm);

/** Register a phylo-HMM object for protection 
    @param A phylo-HMM object */
void phmm_register_protect(PhyloHmm *phmm);

/** \} */

#endif
