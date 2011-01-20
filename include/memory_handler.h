/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#ifndef MEMORY_HANDLER_H
#define MEMORY_HANDLER_H

#include <stdio.h>
#include <stdlib.h>
#include "misc.h"
#include "lists.h"
#include "tree_model.h"
#include "gff.h"
#include "msa.h"
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

void phast_new_mem_handler();
void phast_register_protected_object(void *ptr, void (*function)(void*));
void phast_unregister_protected(void *ptr);
void phast_free_all();
void *smalloc(size_t size);
void *srealloc(void *ptr0, size_t size);
void sfree(void *ptr0);
void set_static_var(void **ptr);

void phast_mem_protect(void *ptr0);
void lst_protect(List *l);
void str_protect(String *s);
void vec_protect(Vector *v);
void zvec_protect(Zvector *v);
void mat_protect(Matrix *m);
void zmat_protect(Zmatrix *m);
void mm_protect(MarkovMatrix *mm);
void gp_protect(GapPatternMap *gpm);
void tm_rmp_protect(TreeModel *tm);
void tm_altmod_protect(AltSubstMod *am);
void tm_protect(TreeModel *tm);
void tm_register_protect(TreeModel *tm);
void gff_feat_protect(GFF_Feature *feat);
void gff_protect(GFF_Set *gff);
void gff_register_protect(GFF_Set *gff);
void cm_protect(CategoryMap *cm);
void cm_register_protect(CategoryMap *cm);
void msa_protect_ss(MSA_SS *ss);
void msa_protect(MSA *msa);
void msa_register_protect(MSA *msa);
void hmm_protect(HMM *hmm);
void hmm_register_protect(HMM *hmm);
void phmm_protect(PhyloHmm *p);
void phmm_register_protect(PhyloHmm *phmm);
void tree_protect(TreeNode *tr);

#endif
