/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <float.h>
#include "phast/stacks.h"
#include "phast/trees.h"
#include "phast/misc.h"
#include "phast/stringsplus.h"
#include "phast/nj.h"
#include "phast/upgma.h"
#include "phast/heap.h"

void upgma_find_min(Matrix *D, Vector *active, int *u, int *v) {
  int i, j, n = D->nrows;
  double min = INFINITY;

  for (i = 0; i < n; i++) {
    if (vec_get(active, i) == FALSE) continue;
    for (j = i+1; j < n; j++) {
      if (vec_get(active, j) == FALSE) continue;
      double d = mat_get(D, i, j);
      if (d < min) {
        min = d;
        *u = i;
        *v = j;
      }
    }
  }

  if (min == INFINITY)
    die("ERROR in upgma_find_min: fewer than two active taxa\n");
}

void upgma_updateD(Matrix *D, int u, int v, int w, Vector *active, Vector *sizes,
                   Vector *heights) {
  double size_u = vec_get(sizes, u);
  double size_v = vec_get(sizes, v);
  double size_w = size_u + size_v;

  for (int k = 0; k < w; k++) {
    if (vec_get(active, k) == FALSE) continue;
    if (k == u || k == v) continue;

    double duk = (u < k ? mat_get(D, u, k) : mat_get(D, k, u));
    double dvk = (v < k ? mat_get(D, v, k) : mat_get(D, k, v));
    double dnew = (size_u * duk + size_v * dvk) / size_w;

    mat_set(D, k, w, dnew);
    if (signbit(mat_get(D, k, w)))
      mat_set(D, k, w, 0);
  }

  double hw = mat_get(D, u, v) / 2.0;
  mat_set(D, u, w, hw - vec_get(heights, u));
  mat_set(D, v, w, hw - vec_get(heights, v));
  vec_set(heights, w, hw);
  vec_set(sizes, w, size_w);

  /* we can't let the distances go negative in this implementation
     because it will mess up the likelihood calculation */
  if (signbit(mat_get(D, u, w))) /* covers -0 case */
    mat_set(D, u, w, 0);
  if (signbit(mat_get(D, v, w)))
    mat_set(D, v, w, 0);
}

/* version of nj_infer_tree simplified to use the UPGMA algorithm and
   return an ultrametric tree. If dt_dD is non-NULL,
   will be populated with Jacobian for 2n-3 branch lengths
   vs. n-choose-2 pairwise distances  */
TreeNode* upgma_infer_tree(Matrix *initD, char **names, Matrix *dt_dD) {
  int n = initD->nrows;
  int N = 2*n - 2;
  int i, j, u = -1, v = -1, w;
  Matrix *D;
  Vector *active, *sizes, *heights;
  List *nodes;
  TreeNode *node_u, *node_v, *node_w, *root;
  double hw;

  if (initD->nrows != initD->ncols || n < 2)
    die("ERROR upgma_infer_tree: bad distance matrix\n");

  D = mat_new(N, N); mat_zero(D);
  active = vec_new(N); vec_set_all(active, FALSE);
  sizes = vec_new(N); vec_zero(sizes);  /* FIXME.  Use list of ints */
  heights = vec_new(N);
  nodes = lst_new_ptr(N);
  tr_reset_id();

  for (i = 0; i < n; i++) {
    node_u = tr_new_node();
    strcat(node_u->name, names[i]);
    lst_push_ptr(nodes, node_u);
    vec_set(active, i, TRUE);
    vec_set(sizes, i, 1.0); /* FIXME */
    for (j = i+1; j < n; j++)
      mat_set(D, i, j, mat_get(initD, i, j));
  }
  
  /* main loop, over internal nodes w */
  for (w = n; w < N; w++) {
    upgma_find_min(D, active, &u, &v);  // find closest pair
    upgma_updateD(D, u, v, w, active, sizes, heights);

    node_w = tr_new_node();
    lst_push_ptr(nodes, node_w);
    node_u = lst_get_ptr(nodes, u);
    node_v = lst_get_ptr(nodes, v);
    tr_add_child(node_w, node_u);
    tr_add_child(node_w, node_v);
    node_u->dparent = mat_get(D, u, w);
    node_v->dparent = mat_get(D, v, w);

    vec_set(active, u, FALSE);
    vec_set(active, v, FALSE);
    vec_set(active, w, TRUE);
  }

  /* there should be exactly two active nodes left. Join them under
     a root. */
  node_u = NULL; node_v = NULL;
  root = tr_new_node();
  for (i = 0; i < N; i++) {
    if (vec_get(active, i) == TRUE) {
      if (node_u == NULL) {
        u = i;
        node_u = lst_get_ptr(nodes, i);
      }
      else if (node_v == NULL) {
        v = i;
        node_v = lst_get_ptr(nodes, i);
      }
      else 
        die("ERROR upgma_infer_tree: more than two nodes left at root\n");
    }
  }
  tr_add_child(root, node_u);
  tr_add_child(root, node_v);

  hw = mat_get(D, u, v) / 2.0;
  node_u->dparent = hw - vec_get(heights, u);
  node_v->dparent = hw - vec_get(heights, v);
  vec_set(heights, root->id, hw);
 
  root->nnodes = N+1;
  tr_set_nnodes(root);

  /* set jacobian; can be done in postprocessing */
  if (dt_dD != NULL)
    upgma_set_dt_dD(root, dt_dD);
  
  lst_free(nodes);
  vec_free(active);
  vec_free(sizes);
  vec_free(heights);
  mat_free(D);

  return root;
}

void upgma_set_dt_dD(TreeNode *tree, Matrix* dt_dD) {
  int i, j, k;
  Matrix *H;
  int nnodes = tree->nnodes, nleaves = (nnodes+2)/2, ndist = nleaves * (nleaves-1) / 2;
    
  /* initialize lists for leaves beneath each node */
  List **leaf_lst = smalloc(nnodes * sizeof(void*));
  for (i = 0; i < nnodes; i++) 
    leaf_lst[i] = lst_new_ptr(nnodes);

  /* populate lists of leaves */
  tr_list_leaves(tree, leaf_lst);
  assert(lst_size(leaf_lst[tree->id]) == nleaves);

  /* now compute node height derivatives */
  H = mat_new(nnodes, ndist); 
  mat_zero(H);
  for (i = 0; i < nnodes; i++) {
    TreeNode *n = lst_get_ptr(tree->nodes, i), *ll, *rl;
    List *lleaves, *rleaves;
    double weight;
    
    if (n->lchild == NULL || n->rchild == NULL)
      continue;

    lleaves = leaf_lst[n->lchild->id];
    rleaves = leaf_lst[n->rchild->id];
    weight = 1.0 / (2.0 * lst_size(lleaves) * lst_size(rleaves));

    for (j = 0; j < lst_size(lleaves); j++) {
      ll = lst_get_ptr(lleaves, j);
      for (k = 0; k < lst_size(rleaves); k++) {
        rl = lst_get_ptr(rleaves, k);
        mat_set(H, i, nj_i_j_to_dist(ll->id, rl->id, nleaves), weight);
      }
    }
  }

  /* finally convert height Jacobian H to branch length Jacobian */
  mat_zero(dt_dD);
  for (i = 0; i < nnodes; i++) {
    TreeNode *n = lst_get_ptr(tree->nodes, i);
    if (n == tree || n == tree->rchild) continue; /* deal with unrooted tree */

    for (j = 0; j < dt_dD->ncols; j++) {
      double val = mat_get(H, n->parent->id, j) - mat_get(H, i, j);
      mat_set(dt_dD, i, j, val);
    }
  }
  
  for (i = 0; i < nnodes; i++) 
    lst_free(leaf_lst[i]);
  free(leaf_lst);
  mat_free(H);
}

UPGMAHeapNode* upgma_heap_node(int i, int j, Matrix *D) {
  UPGMAHeapNode *hn = malloc(sizeof(UPGMAHeapNode));
  hn->i = i;
  hn->j = j;
  hn->val = (i < j) ? mat_get(D, i, j) : mat_get(D, j, i);
  return hn;
}

TreeNode* upgma_fast_infer(Matrix *initD, char **names, Matrix *dt_dD) {
  int n = initD->nrows;
  int N = 2*n - 2;
  int i, j, u = -1, v = -1, w;
  Matrix *D;
  Vector *active, *sizes, *heights;
  List *nodes;
  TreeNode *node_u, *node_v, *node_w, *root;
  HeapNode *heap = NULL;
  UPGMAHeapNode *hn, *newhn;
  double hw;

  if (initD->nrows != initD->ncols || n < 2)
    die("ERROR upgma_fast_infer: bad distance matrix\n");

  D = mat_new(N, N); mat_zero(D);
  active = vec_new(N); vec_set_all(active, FALSE);
  sizes = vec_new(N); vec_zero(sizes);
  heights = vec_new(N);
  nodes = lst_new_ptr(N);
  tr_reset_id();

  /* Initialize leaf nodes and heap */
  for (i = 0; i < n; i++) {
    node_u = tr_new_node();
    strcat(node_u->name, names[i]);
    lst_push_ptr(nodes, node_u);
    vec_set(active, i, TRUE);
    vec_set(sizes, i, 1.0);

    for (j = i+1; j < n; j++) {
      double d = mat_get(initD, i, j);
      mat_set(D, i, j, d);
      hn = upgma_heap_node(i, j, D);
      heap = hp_insert(heap, hn->val, hn);
    }
  }
  
  /* main loop, over internal nodes w */
  for (w = n; w < N; w++) {
    /* Extract minimum from heap */
    while (TRUE) {
      heap = hp_delete_min(heap, (void**)&hn);
      if (vec_get(active, hn->i) == TRUE && vec_get(active, hn->j) == TRUE) break;
      free(hn);
    }

    /* join u and v; w is the new node */
    u = hn->i;
    v = hn->j;
    upgma_updateD(D, u, v, w, active, sizes, heights);
    node_w = tr_new_node();
    lst_push_ptr(nodes, node_w);

    /* attach child nodes to parent and set branch lengths */
    node_u = lst_get_ptr(nodes, u);
    node_v = lst_get_ptr(nodes, v);
    tr_add_child(node_w, node_u);
    tr_add_child(node_w, node_v);
    node_u->dparent = mat_get(D, u, w);
    node_v->dparent = mat_get(D, v, w);

    /* Mark status */
    vec_set(active, u, FALSE);
    vec_set(active, v, FALSE);
    vec_set(active, w, TRUE);

    /* Insert new distances into heap */
    for (i = 0; i < w; i++) {
      if (vec_get(active, i)) {
        newhn = upgma_heap_node(i, w, D);
        heap = hp_insert(heap, newhn->val, newhn);
      }
    }

    free(hn);
  }

  /* Final join */
  node_u = node_v = NULL;
  root = tr_new_node();
  for (i = 0; i < N; i++) {
    if (vec_get(active, i) == TRUE) {
      if (node_u == NULL) {
        u = i;
        node_u = lst_get_ptr(nodes, i);
      }
      else if (node_v == NULL) {
        v = i;
        node_v = lst_get_ptr(nodes, i);
      }
      else 
        die("ERROR upgma_fast_infer: more than two nodes left at root\n");
    }
  }
  tr_add_child(root, node_u);
  tr_add_child(root, node_v);

  hw = mat_get(D, u, v) / 2.0;
  node_u->dparent = hw - vec_get(heights, u);
  node_v->dparent = hw - vec_get(heights, v);
  vec_set(heights, root->id, hw);

  root->nnodes = N + 1;
  tr_set_nnodes(root);

  if (dt_dD != NULL)
    upgma_set_dt_dD(root, dt_dD);  // Postprocess

  lst_free(nodes);
  vec_free(active);
  vec_free(sizes);
  vec_free(heights);
  mat_free(D);

  return root;
}
