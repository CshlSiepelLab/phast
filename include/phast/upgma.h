/* PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file upgma.h
    Simple UPGMA tree inference  
    @ingroup phylo
*/

#ifndef UPGMA_H
#define UPGMA_H

#include <stdio.h>
#include <phast/matrix.h>
#include <phast/trees.h>
#include <phast/tree_model.h>
#include <phast/misc.h>

/* for use with min-heap in fast upgma algorithm */
typedef struct {
  int i, j;
  double val;  // raw distance D(i, j)
} UPGMAHeapNode;

void upgma_find_min(Matrix *D, Vector *active, int *u, int *v);

void upgma_updateD(Matrix *D, int u, int v, int w, Vector *active,
                   Vector *sizes, Vector *heights);

TreeNode* upgma_infer_tree(Matrix *initD, char **names, Matrix *dt_dD);

void upgma_set_dt_dD(TreeNode *tree, Matrix* dt_dD);

UPGMAHeapNode* upgma_heap_node(int i, int j, Matrix *D);

TreeNode* upgma_fast_infer(Matrix *initD, char **names, Matrix *dt_dD);

#endif
