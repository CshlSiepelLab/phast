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
#include "phast/tree_likelihoods.h"
#include "phast/eigen.h"
#include "phast/sufficient_stats.h"
#include "phast/markov_matrix.h"
#include "phast/sparse_matrix.h"
#include "phast/lists.h"
#include "phast/mvn.h"
#include "phast/multi_mvn.h"
#include "phast/heap.h"
#include "phast/crispr.h"
#include "phast/upgma.h"

/* uncomment to dump gradients to a file called "grads_log.txt" */
//#define DUMPGRAD 1


/* Reset Q matrix based on distance matrix.  Assume upper triangular
   square Q and D.  Only touches active rows and columns of Q and D up
   to maxidx. As a side-effect set u and v to the indices of the
   closest neighbors.  Also update sums to sum of distances from each
   node */
void nj_resetQ(Matrix *Q, Matrix *D, Vector *active, Vector *sums, int *u,
               int *v, int maxidx) {
  int i, j, n = 0;
  double min = INFINITY;
  
  if ((D->nrows != D->ncols) || (D->nrows != Q->nrows) ||
      (D->nrows != Q->ncols))
    die("ERROR nj_setQ: dimension mismatch\n");

  *u = *v = 0;
  
  /* update row sums */
  vec_zero(sums);
  for (i = 0; i < maxidx; i++) {
    if (vec_get(active, i) == TRUE) {
      n++;
      for (j = i+1; j < maxidx; j++) {
        if (vec_get(active, j) == TRUE) {
          sums->data[i] += mat_get(D, i, j);
          sums->data[j] += mat_get(D, i, j);
        }
      }
    }
  }
  
  /* now reset Q */
  for (i = 0; i < maxidx; i++) {
    if (vec_get(active, i) == TRUE) {
      for (j = i+1; j < maxidx; j++) {
        if (vec_get(active, j) == TRUE) {
          double qij = (n-2) * mat_get(D, i, j) - vec_get(sums, i) -
            vec_get(sums, j);
          mat_set(Q, i, j, qij);
          if (qij < min) {
            min = qij;
            *u = i;
            *v = j;
          }
        }
      }
    }
  }

  if (min == INFINITY)
    die("ERROR in nj_resetQ: fewer than two active taxa\n");
}

/* Update distance matrix after operation that joins neighbors u and v
   and adds new node w.  Update active list accordingly. Assumes u < v
   < w.  Also assumes all nodes > w are inactive. Assumes sums are
   precomputed */
void nj_updateD(Matrix *D, int u, int v, int w, Vector *active, Vector *sums) {
  int k;
  int n = vec_sum(active);

  if (D->nrows != D->ncols)
    die("ERROR nj_updateD: dimension mismatch\n");
  if (v <= u || w <= v)
    die("ERROR nj_updateD: indices out of order\n");
  if (n <= 2)
    die("ERROR nj_updateD: too few active nodes\n");
  
  mat_set(D, u, w, 0.5 * mat_get(D, u, v) +
	  1.0/(2.0*(n-2)) * (vec_get(sums, u) - vec_get(sums, v)));

  mat_set(D, v, w, mat_get(D, u, v) - mat_get(D, u, w));

  /* we can't let the distances go negative in this implementation
     because it will mess up the likelihood calculation */
  if (signbit(mat_get(D, u, w))) /* covers -0 case */
    mat_set(D, u, w, 0);
  if (signbit(mat_get(D, v, w)))
    mat_set(D, v, w, 0);
  
  for (k = 0; k < w; k++) {
    if (vec_get(active, k) == TRUE && k != u && k != v) {
      double du, dv;

      /* needed because of upper triangular constraint */
      du = (u < k ? mat_get(D, u, k) : mat_get(D, k, u));
      dv = (v < k ? mat_get(D, v, k) : mat_get(D, k, v));
      
      mat_set(D, k, w, 0.5 * (du + dv - mat_get(D, u, v)));

      if (mat_get(D, k, w) < 0)
        mat_set(D, k, w, 0);
    }
  }
}


/* Main function to infer the tree from a starting distance matrix.
   Does not alter the provided distance matrix.  If dt_dD is non-NULL,
   will be populated with Jacobian for 2n-3 branch lengths
   vs. n-choose-2 pairwise distances */
TreeNode* nj_infer_tree(Matrix *initD, char **names, Matrix *dt_dD) {
    int n = initD->nrows;
    int N = 2*n - 2;   /* number of nodes in unrooted tree */
    int i, j, u = -1, v = -1, w;
    Matrix *D, *Q;
    Vector *sums, *active;
    List *nodes;  /* could just be an array */
    TreeNode *node_u, *node_v, *node_w, *root;
    int npairs = n * (n-1) / 2, Npairs = N * (N-1) / 2;
    double *Jk = NULL, *Jnext = NULL;
    
    if (initD->nrows != initD->ncols || n < 3)
      die("ERROR nj_infer_tree: bad distance matrix\n");

    if (dt_dD != NULL && (dt_dD->nrows != N || dt_dD->ncols != npairs))
      die("ERROR nj_infer_tree: bad dimension in dt_dD\n");
    
    /* create a larger distance matrix of dimension N x N to
       accommodate internal nodes; also set up list of active nodes
       and starting tree nodes */
    D = mat_new(N, N); mat_zero(D);
    active = vec_new(N); vec_set_all(active, FALSE);
    sums = vec_new(N); vec_zero(sums);
    nodes = lst_new_ptr(N);
    tr_reset_id();
    
    for (i = 0; i < n; i++) {
      node_u = tr_new_node();
      strcat(node_u->name, names[i]);
      lst_push_ptr(nodes, node_u);
      vec_set(active, i, TRUE);
      for (j = i+1; j < n; j++)
        mat_set(D, i, j, mat_get(initD, i, j));
    }
   
    /* set up Q */
    Q = mat_new(N, N); mat_zero(Q);

    /* set up backprop data */
    if (dt_dD != NULL) {
      Jk = malloc(Npairs * npairs * sizeof(double));
      Jnext = malloc(Npairs * npairs * sizeof(double));
      nj_backprop_init(Jk, n);
      mat_zero(dt_dD);
    }
    
    /* main loop, over internal nodes w */
    for (w = n; w < N; w++) {   
      nj_resetQ(Q, D, active, sums, &u, &v, w);
      
      nj_updateD(D, u, v, w, active, sums);                    
      node_w = tr_new_node();
      lst_push_ptr(nodes, node_w);

      /* attach child nodes to parent and set branch lengths */
      node_u = lst_get_ptr(nodes, u);
      node_v = lst_get_ptr(nodes, v);
      tr_add_child(node_w, node_u);
      tr_add_child(node_w, node_v);
      node_u->dparent = mat_get(D, u, w);
      node_v->dparent = mat_get(D, v, w);

      if (dt_dD != NULL) {
        nj_backprop_set_dt_dD(Jk, dt_dD, n, u, v, node_u->id, node_v->id, active);
        nj_backprop(Jk, Jnext, n, u, v, w, active);
      }

      /* this has to be done after the backprop calls */
      vec_set(active, u, FALSE);
      vec_set(active, v, FALSE);
      vec_set(active, w, TRUE);

      if (dt_dD != NULL) {
        free(Jk);
        Jk = Jnext;
        Jnext = malloc(Npairs * npairs * sizeof(double));
      }
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
          die("ERROR nj_infer_tree: more than two nodes left at root\n");
      }
    }
    tr_add_child(root, node_u);
    tr_add_child(root, node_v);
    node_u->dparent = mat_get(D, u, v) / 2;
    node_v->dparent = mat_get(D, u, v) / 2;

    if (dt_dD != NULL) 
      nj_backprop_set_dt_dD(Jk, dt_dD, n, u, v, node_u->id, node_v->id, active);
    
    /* finish set up of tree */
    root->nnodes = N+1;
    tr_reset_nnodes(root);

    assert(root->id == root->nnodes - 1); /* important for indexing */
    
    lst_free(nodes);
    vec_free(active);
    vec_free(sums);
    mat_free(D);
    mat_free(Q);

    if (dt_dD != NULL) {
      free(Jk);
      free(Jnext);
    }
    
    return root;
}

/* Faster version of function to infer the tree from a starting
   distance matrix. Uses a min-heap for efficient lookup of minimum Q
   values, with lazy evaluation to avoid unnecessary computations.
   Does not alter the provided distance matrix. If dt_dD is non-NULL,
   will be populated with Jacobian for 2n-3 branch lengths
   vs. n-choose-2 pairwise distances */
TreeNode* nj_fast_infer(Matrix *initD, char **names, Matrix *dt_dD) {
  int n = initD->nrows, orign = n;
  int N = 2*n - 2;   /* number of nodes in unrooted tree */
  int i, j, u = -1, v = -1, w;
  Matrix *D;
  Vector *sums, *active;
  List *nodes;  
  TreeNode *node_u, *node_v, *node_w, *root;
  HeapNode *heap = NULL;
  NJHeapNode *hn, *newhn;
  int rev[N];
  int npairs = n * (n-1) / 2, Npairs = N * (N-1) / 2;
  //  double *Jk = NULL, *Jnext = NULL;
  static SparseMatrix *Jk = NULL, *Jnext = NULL;
    
  if (initD->nrows != initD->ncols || n < 3)
    die("ERROR nj_fast_infer: bad distance matrix\n");

  if (dt_dD != NULL && (dt_dD->nrows != N || dt_dD->ncols != npairs))
    die("ERROR nj_fast_infer: bad dimension in dt_dD\n");
    
  /* initialize revision numbers for all nodes */
  for (i = 0; i < N; i++) rev[i] = 0;

  /* create a larger distance matrix of dimension N x N to
     accommodate internal nodes; also set up list of active nodes
     and starting tree nodes */
  D = mat_new(N, N); mat_zero(D);
  active = vec_new(N); vec_set_all(active, FALSE);
  sums = vec_new(N); vec_zero(sums);
  nodes = lst_new_ptr(N);
  tr_reset_id();
    
  for (i = 0; i < n; i++) {
    node_u = tr_new_node();
    strcat(node_u->name, names[i]);
    lst_push_ptr(nodes, node_u);
    vec_set(active, i, TRUE);
    for (j = i+1; j < n; j++) {
      double d = mat_get(initD, i, j);
      mat_set(D, i, j, d);
      sums->data[i] += d;
      sums->data[j] += d;
    }
  }

  /* also set up the heap for Q values */
  for (i = 0; i < n; i++) {
    for (j = i+1; j < n; j++) {
      hn = nj_heap_computeQ(i, j, n, D, sums, rev);
      heap = hp_insert(heap, hn->val, hn);
    }
  }

  /* set up backprop data */
  if (dt_dD != NULL) {
    if (Jk == NULL) { /* first call */
      Jk = spmat_new(Npairs, npairs, 100);
      Jnext = spmat_new(Npairs, npairs, 100);
    }
    else
      assert(Npairs == Jk->nrows && npairs == Jk->ncols);

    nj_backprop_init_sparse(Jk, n);
    mat_zero(dt_dD);
  }
    
  /* main loop, over internal nodes w */
  for (w = n; w < N; w++) {

    /* get the minimum Q value from the heap; use lazy evaluation */
    while (TRUE) {
      heap = hp_delete_min(heap, (void**)&hn);

      if (vec_get(active, hn->i) == FALSE || vec_get(active, hn->j) == FALSE) {
        free(hn);
        continue;
      }
      else if (hn->rev_i == rev[hn->i] && hn->rev_j == rev[hn->j]) 
        break; /* valid and active */
      else {
        /* active but stale; recompute */
        newhn = nj_heap_computeQ(hn->i, hn->j, n, D, sums, rev);
        heap = hp_insert(heap, newhn->val, newhn);
        free(hn);
      }
    }
      
    /* join u and v; w is the new node */
    u = hn->i;
    v = hn->j;
    nj_updateD(D, u, v, w, active, sums);
    node_w = tr_new_node();
    lst_push_ptr(nodes, node_w);

    /* attach child nodes to parent and set branch lengths */
    node_u = lst_get_ptr(nodes, u);
    node_v = lst_get_ptr(nodes, v);
    tr_add_child(node_w, node_u);
    tr_add_child(node_w, node_v);
    node_u->dparent = mat_get(D, u, w);
    node_v->dparent = mat_get(D, v, w);
    
    /* update row sums and revision numbers */
    vec_set(sums, w, 0);  
    for (i = 0; i < w; i++) {
      if (vec_get(active, i) == TRUE && i != u && i != v) {   
        double du = (u < i ? mat_get(D, u, i) : mat_get(D, i, u)); /* upper triangular */
        double dv = (v < i ? mat_get(D, v, i) : mat_get(D, i, v));
        sums->data[i] += (mat_get(D, i, w) - du - dv); /* can be updated */
        sums->data[w] += mat_get(D, i, w); /* have to compute from scratch */
      }
      rev[i]++;
    }
    rev[w]++;
    
    if (dt_dD != NULL) {
            nj_backprop_set_dt_dD_sparse(Jk, dt_dD, orign, u, v, node_u->id, node_v->id, active);
            nj_backprop_sparse(Jk, Jnext, orign, hn->i, hn->j, w, active);
    }

    /* this has to be done after the backprop calls */
    vec_set(active, u, FALSE);
    vec_set(active, v, FALSE);
    vec_set(active, w, TRUE);
    n--;  /* one fewer active nodes */

    if (dt_dD != NULL) {
      /* swap pointers to avoid deep copy */
      SparseMatrix *tmp = Jk; 
      Jk = Jnext; 
      Jnext = tmp;
    }
      
    /* finally, add new Q values to the heap */
    for (i = 0; i < w; i++) {
      if (vec_get(active, i) == TRUE) {
        newhn = nj_heap_computeQ(i, w, n, D, sums, rev); /* add to heap */
        heap = hp_insert(heap, newhn->val, newhn);
      }
    }

    free(hn);
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
        die("ERROR nj_fast_infer: more than two nodes left at root\n");
    }
  }
  tr_add_child(root, node_u);
  tr_add_child(root, node_v);
  node_u->dparent = mat_get(D, u, v) / 2;
  node_v->dparent = mat_get(D, u, v) / 2;

  if (dt_dD != NULL) 
    nj_backprop_set_dt_dD_sparse(Jk, dt_dD, orign, u, v, node_u->id, node_v->id, active);
  //nj_backprop_set_dt_dD(Jk, dt_dD, orign, u, v, node_u->id, node_v->id, active);
  
  /* finish set up of tree */
  root->nnodes = N+1;
  tr_reset_nnodes(root);

  assert(root->id == root->nnodes - 1); /* important for indexing */

  /* drain heap */
  while (heap != NULL) {
    heap = hp_delete_min(heap, (void**)&hn);
    free(hn);
  }
  
  hp_free(heap);
  lst_free(nodes);
  vec_free(active);
  vec_free(sums);
  mat_free(D);
    
  return root;
}

NJHeapNode* nj_heap_computeQ(int i, int j, int n, Matrix *D, Vector *sums, int *rev) {
    NJHeapNode *hn = malloc(sizeof(NJHeapNode));
    hn->i = i;
    hn->j = j;
    hn->val = (n - 2) * mat_get(D, i, j) - vec_get(sums, i) -
      vec_get(sums, j);
    hn->rev_i = rev[i];
    hn->rev_j = rev[j];
    return hn;
}

/* compute pairwise distance between two DNA seqs using the
   Jukes-Cantor model */
double nj_compute_JC_dist(MSA *msa, int i, int j) {
  int k, diff = 0, n = 0;
  double d;
  for (k = 0; k < msa->length; k++) {
    if (msa->seqs[i][k] == GAP_CHAR || msa->seqs[j][k] == GAP_CHAR ||
	msa->is_missing[(int)msa->seqs[i][k]] ||
	msa->is_missing[(int)msa->seqs[j][k]])
      continue;
    n++;
    if (msa->seqs[i][k] != msa->seqs[j][k])
      diff++;
  }
  if ((double)diff/n >= 0.75)
    /* in this case, there are too many differences for the correction
       to work.  We basically want the distance to be "really long" but
       not so long that it completely skews the tree or makes it
       impossible for the variational algorithm to recover.  */
    d = 3.0;
  else
    d = -0.75 * log(1 - 4.0/3 * diff/n);   /* Jukes-Cantor correction */

  assert(isfinite(d));
  return d;
}

/* based on a multiple alignment, build and return a distance matrix
   using the Jukes-Cantor model.  Assume DNA alphabet */
Matrix *nj_compute_JC_matr(MSA *msa) {
  int i, j;
  Matrix *retval = mat_new(msa->nseqs, msa->nseqs);

  mat_zero(retval);
  
  for (i = 0; i < msa->nseqs; i++) 
    for (j = i+1; j < msa->nseqs; j++) 
      mat_set(retval, i, j,
              nj_compute_JC_dist(msa, i, j));

  return retval;  
}

/* compute a distance matrix from a tree, defining each pairwise
   distance as the edge length between the corresponding taxa */ 
Matrix *nj_tree_to_distances(TreeNode *tree, char **names, int n) {
  TreeNode *n1, *n2;
  List *leaves = lst_new_ptr(tree->nnodes);
  int i, j, ii, jj;
  Matrix *D;
  double dist;
  int *seq_idx;
  unsigned int all_zeroes = TRUE;
  
  assert(tree->nodes != NULL);  /* assume list of nodes exists */
  
  for (i = 0; i < tree->nnodes; i++) {
    n1 = lst_get_ptr(tree->nodes, i);
    if (all_zeroes == TRUE && n1->dparent > 0) all_zeroes = FALSE;
    if (n1->lchild == NULL && n1->rchild == NULL)
      lst_push_ptr(leaves, n1);
  }
  
  if (lst_size(leaves) != n)
    die("ERROR in nj_tree_to_distances: number of names must match number of leaves in tree.\n");

  /* if the input tree had no branch lengths defined, all will have
     values of zero, which will be a problem.  In this case,
     just initialize them all to a small constant value */
  if (all_zeroes == TRUE) {
    for (i = 0; i < tree->nnodes; i++) {
      n1 = lst_get_ptr(tree->nodes, i);
      if (n1->parent != NULL)
        n1->dparent = 0.1;
    }
  }
  
  D = mat_new(n, n);
  mat_zero(D);

  seq_idx = nj_build_seq_idx(leaves, names);

  /* O(n^2) operation but seems plenty fast in practice */
  for (i = 0; i < lst_size(leaves); i++) {
    n1 = lst_get_ptr(leaves, i);
    ii = seq_idx[n1->id];   /* convert to seq indices for matrix */
    for (j = i+1; j < lst_size(leaves); j++) {
      n2 = lst_get_ptr(leaves, j);
      jj = seq_idx[n2->id];
      dist = nj_distance_on_tree(tree, n1, n2);
      if (ii < jj)
        mat_set(D, ii, jj, dist);
      else
        mat_set(D, jj, ii, dist);
    }
  }

  lst_free(leaves);
  sfree(seq_idx);
  
  return(D);
}

double nj_distance_on_tree(TreeNode *root, TreeNode *n1, TreeNode *n2) {
  double dist[root->nnodes];
  int id;
  TreeNode *n;
  double totd1, totd2;
  
  /* initialize distance from n1 to each ancestor to be -1 */
  for (id = 0; id < root->nnodes; id++)
    dist[id] = -1;
  
  /* find distance to each ancestor of n1 */
  for (n = n1, totd1 = 0; n->parent != NULL; n = n->parent) {
    totd1 += n->dparent;
    dist[n->parent->id] = totd1;
  }
  dist[root->id] = totd1;

  /* now trace ancestry of n2 until an ancestor of n1 is found */
  for (n = n2, totd2 = 0; dist[n->id] == -1 && n->parent != NULL; n = n->parent) 
    totd2 += n->dparent;

  if (n->parent == NULL && dist[n->id] == -1)
    die("ERROR in nj_distance_on_tree: got to root without finding LCA\n");
  
  /* at this point, it must be true that n is the LCA of n1 and n2 */ 
  return totd2 + dist[n->id];
  
}

/* convert an nd-dimensional vector to an nxn upper triangular
   distance matrix.  Assumes each taxon is represented as a point in
   d-dimensional space.  Wrapper for versions that assume either
   Euclidean or hyperbolic geometry */
void nj_points_to_distances(Vector *points, CovarData *data) {
  if (data->hyperbolic)
    nj_points_to_distances_hyperbolic(points, data);
  else
    nj_points_to_distances_euclidean(points, data);
}
  
/* convert an nd-dimensional vector to an nxn upper triangular
   distance matrix.  Assumes each taxon is represented as a point in
   d-dimensional space and computes Euclidean distances between these
   points */ 
void nj_points_to_distances_euclidean(Vector *points, CovarData *data) {
  int i, j, k, vidx1, vidx2, n, d;
  double sum;
  Matrix *D = data->dist;

  n = D->nrows;
  d = points->size / n;
  
  if (points->size != n*d || D->nrows != D->ncols) 
    die("ERROR nj_points_to_distances_euclidean: bad dimensions\n");

  mat_zero(D);
  for (i = 0; i < n; i++) {
    vidx1 = i*d;
    for (j = i+1; j < n; j++) {
      vidx2 = j*d;
      sum = 0;
      for (k = 0; k < d; k++) {
        double diff = vec_get(points, vidx1 + k) -
          vec_get(points, vidx2 + k);
        sum += diff*diff;
      }
      mat_set(D, i, j, sqrt(sum) / data->pointscale);
    }
  }
}

/* convert an nd-dimensional vector to an nxn upper triangular
   distance matrix.  Assumes each taxon is represented as a point in
   d-dimensional space and computes hyperbolic distances between these
   points */ 
void nj_points_to_distances_hyperbolic(Vector *points, CovarData *data) {
  int i, j, k, vidx1, vidx2, n, d;
  double lor_inner, ss1, ss2, x0_1, x0_2, Dij, u;
  Matrix *D = data->dist;
  double alpha = 1.0 / sqrt(data->negcurvature);   /* curvature radius */

  n = D->nrows;
  d = points->size / n;
  
  if (points->size != n*d || D->nrows != D->ncols) {
    die("ERROR nj_points_to_distances_hyperbolic: bad dimensions\n");
  }

  mat_zero(D);
  for (i = 0; i < n; i++) {
    vidx1 = i*d;
    for (j = i+1; j < n; j++) {
      vidx2 = j*d;
      lor_inner = 0;
      ss1 = 1;
      ss2 = 1;
      for (k = 0; k < d; k++) {
        double xi = vec_get(points, vidx1 + k);
        double xj = vec_get(points, vidx2 + k);
        lor_inner += xi * xj ;
        ss1 += xi * xi;
        ss2 += xj * xj;
      }
      x0_1 = sqrt(ss1); /* the 0th dimension for each point is determined by the
                           others, to stay on the hyperboloid */
      x0_2 = sqrt(ss2);

      lor_inner -= x0_1 * x0_2;  /* last term of Lorentz inner product */

      u = -lor_inner;
      Dij = (alpha / data->pointscale) * acosh_stable(u);

      mat_set(D, i, j, Dij);      
      assert(isfinite(Dij) && Dij >= 0);
    }
  }
}

/* compute the gradient of the log likelihood for a tree model with
   respect to the free parameters of the MVN averaging distribution,
   starting from a given MVN sample (points). Returns log likelihood
   of current model, which is computed as a by-product.  */
double nj_compute_model_grad(TreeModel *mod, multi_MVN *mmvn, 
                             Vector *points, Vector *points_std,
                             Vector *grad, CovarData *data) {
  int n = data->nseqs; /* number of taxa */
  int d = mmvn->n * mmvn->d / n; /* dimensionality; have to accommodate diagonal case */
  int dim = n*d; /* full dimension of point vector */
  int i, j, k;
  double porig, ll_base, loglambda_grad;
  Vector *dL_dx = vec_new(dim);
  
  if (grad->size != dim + data->params->size)
    die("ERROR in nj_compute_model_grad: bad gradient dimension.\n");
  
  if (data->type == DIST && data->Lapl_pinv_evals == NULL)
    die("ERROR in nj_compute_model_grad: Laplacian pseudoinverse and eigendecomposition required in DIST case.\n");
  else if (data->type == LOWR && data->R == NULL)
    die("ERROR in nj_compute_model_grad: low-rank matrix R required in LOWR case.\n");
  
  /* obtain gradient with respect to points, dL/dx */
  ll_base = nj_dL_dx_smartest(points, dL_dx, mod, data);

  /* TEMPORARY: compare with numerical version */
  /* printf("analytical (%f):\n", ll_base); */
  /* vec_print(dL_dx, stdout); */
  /* ll_base = nj_dL_dx_dumb(points, dL_dx, mod, data); */
  /* printf("numerical (%f):\n", ll_base); */
  /* vec_print(dL_dx, stdout); */
  /* exit(0); */

  if (!isfinite(ll_base)) /* can happen with crispr model; force calling code to deal with it */
    return ll_base;
        
  /* now derive partial derivatives wrt free parameters from dL/dx */
  vec_zero(grad);
  loglambda_grad = 0;
  for (i = 0; i < n; i++) {  
    for (k = 0; k < d; k++) {
      int pidx = i*d + k;
      porig = vec_get(points, pidx);

      /* the partial derivative wrt the mean parameter is equal to the
         derivative with respect to the point, because the point is just a
         translation of a 0-mean MVN variable via the reparameterization
         trick */
      vec_set(grad, pidx, vec_get(dL_dx, pidx));

      /* the partial derivative wrt the variance parameter is more
         complicated because of the reparameterization trick */      
      if (data->type == CONST || data->type == DIST)
        loglambda_grad += 0.5 * vec_get(dL_dx, pidx) * (porig - mmvn_get_mu_el(mmvn, pidx));
      /* assumes log parameterization of scale factor lambda */
      
      else if (data->type == DIAG) 
        /* in the DIAG case, the partial derivative wrt the
           corresponding variance parameter can be computed directly
           based on a single point and coordinate */
        vec_set(grad, (i+n)*d + k, 0.5 * vec_get(dL_dx, pidx) * (porig - mmvn_get_mu_el(mmvn, pidx)));
      
    }
  }
  if (data->type == CONST || data->type == DIST) /* in this case, need to update the final
                                                    gradient component corresponding to the
                                                    lambda parameter */
    vec_set(grad, dim, loglambda_grad); 

  else if (data->type == LOWR) { /* in this case have to sum across
                                    dimensions because there is a
                                    many-to-many relationship with the
                                    variance parameters */
    for (i = 0; i < n; i++) 
      for (j = 0; j < data->lowrank; j++) 
        for (k = 0; k < d; k++) 
          vec_set(grad, dim + i*data->lowrank + j, vec_get(grad, dim + i*data->lowrank + j) +
                  vec_get(grad, i*d + k) * vec_get(points_std, j*d + k));
          /* update for parameter corresponding to element R[i, j],
             which has index in grad of dim + (i*data->lowrank) + j.
             The gradient is a dot product of dL/dx and dx/dp, where L
             is the log likelihood, x is the point, and p is the
             variance parameter.  In this case, dxj/dp is simply zj, the
             corresponding standardized variable.
          */
  }

  vec_free(dL_dx);
  return ll_base;
}  

/* version of nj_compute_model_grad that calculates all derivatives
   numerically, to check correctness of analytical calculations.  The
   check knows nothing about the details of the parameterization; it
   just uses a brute force numerical calculation */
double nj_compute_model_grad_check(TreeModel *mod, multi_MVN *mmvn, 
                                   Vector *points, Vector *points_std,
                                   Vector *grad, CovarData *data) {
  int n = data->msa->nseqs; /* number of taxa */
  int d = mmvn->n * mmvn->d / n; /* dimensionality; have to accommodate diagonal case */
  int dim = n*d; /* full dimension of point vector */
  int i, j;
  double porig, ll_base, ll, deriv;
  TreeNode *tree, *orig_tree;   /* has to be rebuilt repeatedly; restore at end */
  Vector *points_tweak = vec_new(points->size);
  Vector *sigmapar = data->params;
  Vector *dL_dx = vec_new(points->size);
  Matrix *D = data->dist;
  
  if (grad->size != dim + data->params->size)
    die("ERROR in nj_compute_model_grad_check: bad gradient dimension.\n");
  
  /* set up tree model and get baseline log likelihood */
  nj_points_to_distances(points, data);    
  tree = nj_inf(D, data->names, NULL, data);
  orig_tree = tr_create_copy(tree);   /* restore at the end */
  nj_reset_tree_model(mod, tree);
  ll_base = nj_compute_log_likelihood(mod, data, NULL);

  if (!isfinite(ll_base)) /* can happen with crispr model; force
                             calling code to deal with it */
    return ll_base;
  
  /* Now perturb each point and propagate perturbation through distance
     calculation, neighbor-joining reconstruction, and likelihood
     calculation on tree */
  for (i = 0; i < dim; i++) {
    double mu_orig, dxi_dmui;

    porig = vec_get(points, i);
    vec_set(points, i, porig + DERIV_EPS);

    nj_points_to_distances(points, data); 
    tree = nj_inf(D, data->names, NULL, data);
    nj_reset_tree_model(mod, tree);      
    ll = nj_compute_log_likelihood(mod, data, NULL);

    if (!isfinite(ll)) /* can happen with crispr model; force
                          calling code to deal with it */
      return ll;

    deriv = (ll - ll_base) / DERIV_EPS; 
    vec_set(dL_dx, i, deriv); /* derive of log likelihood wrt dim i of
                                 point; need to save this for later */
    
    /* the mean is straightforward but check it anyway to be sure */
    mu_orig = mmvn_get_mu_el(mmvn, i);
    mmvn_set_mu_el(mmvn, i, mu_orig + DERIV_EPS);
    vec_copy(points_tweak, points_std);
    mmvn_map_std(mmvn, points_tweak);
    dxi_dmui = (vec_get(points_tweak, i) - porig) / DERIV_EPS;
    /* should be 1! */
    vec_set(grad, i, deriv * dxi_dmui);
    mmvn_set_mu_el(mmvn, i, mu_orig);
            
    vec_set(points, i, porig); /* restore orig */
  }

  /* for the variance parameters, we have to consider potential
     changes to entire vector of points */
  for (i = 0; i < sigmapar->size; i++) {
    double dL_dp = 0, origp = vec_get(sigmapar, i);    
    vec_set(sigmapar, i, origp + DERIV_EPS);
    nj_update_covariance(mmvn, data);
    vec_copy(points_tweak, points_std);
    mmvn_map_std(mmvn, points_tweak);
    vec_minus_eq(points_tweak, points);
    vec_scale(points_tweak, 1.0/DERIV_EPS); /* now contains dx / dp
                                               where p is the variance
                                               parameter */

    /* dL/dp is a dot product of dL/dx and dx/dp */
    for (j = 0; j < points_tweak->size; j++)
      dL_dp += vec_get(dL_dx, j) * vec_get(points_tweak, j);
    
    vec_set(grad, i + dim, dL_dp);
    vec_set(sigmapar, i, origp); /* restore orig */
  }
  
  nj_update_covariance(mmvn, data); /* make sure to leave it in original state */

  nj_reset_tree_model(mod, orig_tree);
  vec_free(points_tweak);
  vec_free(dL_dx);
  return ll_base;
}  

/* optimize variational model by stochastic gradient ascent using the
   Adam algorithm.  Takes initial tree model and alignment and
   distance matrix, dimensionality of Euclidean space to work in.
   Note: alters distance matrix */
void nj_variational_inf(TreeModel *mod, multi_MVN *mmvn,
                        int nminibatch, double learnrate, int nbatches_conv,
                        int min_nbatches, CovarData *data, FILE *logf) {

  Vector *points, *grad, *kldgrad, *avegrad, *m, *m_prev, *v, *v_prev,
    *best_mu, *best_sigmapar, *rescaledgrad, *sparsitygrad = NULL, 
    *points_std;
  Vector *sigmapar = data->params;
  int n = data->nseqs, i, j, t, stop = FALSE, bestt = -1, graddim,
    dim = data->dim, fulld = n*dim;
  double ll, avell, kld, bestelb = -INFTY, bestll = -INFTY, bestkld = -INFTY,
    running_tot = 0, last_running_tot = -INFTY, trace, logdet, penalty = 0,
    bestpenalty = 0, ave_nf_logdet, best_nf_logdet = -INFTY;
  FILE *gradf = NULL;

  /* for nuisance parameters; these are parameters that are optimized
     by stochastic gradient descent but are not fully sampled via the
     variational distribution */
  int n_nuisance_params = nj_get_num_nuisance_params(mod, data);
  Vector *nuis_grad = NULL, *ave_nuis_grad = NULL, *m_nuis = NULL,
    *v_nuis = NULL, *m_nuis_prev = NULL, *v_nuis_prev = NULL,
    *best_nuis_params = NULL;
  
#ifdef DUMPGRAD
  gradf = phast_fopen("grads_log.txt", "w");
#endif
  
  if (mmvn->d * mmvn->n != dim * n)
    die("ERROR in nj_variational_inf: bad dimensions\n");

  points = vec_new(fulld);
  graddim = fulld + data->params->size;
  grad = vec_new(graddim);  
  kldgrad = vec_new(graddim);
  avegrad = vec_new(graddim);
  rescaledgrad = vec_new(graddim);
  m = vec_new(graddim);
  m_prev = vec_new(graddim);
  v = vec_new(graddim);
  v_prev = vec_new(graddim);

  if (n_nuisance_params > 0) {
    nuis_grad = vec_new(n_nuisance_params);
    ave_nuis_grad = vec_new(n_nuisance_params);
    m_nuis = vec_new(n_nuisance_params);
    v_nuis = vec_new(n_nuisance_params);
    m_nuis_prev = vec_new(n_nuisance_params);
    v_nuis_prev = vec_new(n_nuisance_params);
    best_nuis_params = vec_new(n_nuisance_params);
  }

  if (data->type == LOWR) /* in this case, the underlying standard
                             normal MVN is of the lower dimension */
    points_std = vec_new(data->lowrank * dim);
  else
    points_std = vec_new(fulld);

  if (data->type == LOWR && data->sparsity != -1)
    sparsitygrad = vec_new(graddim);
  
  best_mu = vec_new(fulld);
  mmvn_save_mu(mmvn, best_mu);
  best_sigmapar = vec_new(sigmapar->size);
  vec_copy(best_sigmapar, sigmapar);

  /* set up log file */
  if (logf != NULL) {
    fprintf(logf, "# varPHAST logfile\n");
    fprintf(logf, "state\tll\tkld\tnfld\telb\t");
    if (data->type == LOWR && data->sparsity != -1)
      fprintf(logf, "penalty\t");
    for (j = 0; j < fulld; j++)
      fprintf(logf, "mu.%d\t", j);
    for (j = 0; j < sigmapar->size; j++)
      fprintf(logf, "sigma.%d\t", j);
    for (j = 0; j < n_nuisance_params; j++)
      fprintf(logf, "%s\t", nj_get_nuisance_param_name(mod, data, j));
    fprintf(logf, "\n");
  }
  if (gradf != NULL) {
    fprintf(gradf, "state\t");
    for (j = 0; j < fulld; j++)
      fprintf(gradf, "g_mu.%d\t", j);
    for (j = fulld; j < grad->size; j++)
      fprintf(gradf, "g_sig.%d\t", j-fulld);
    fprintf(gradf, "\n");
  }
  
  /* initialize moments for Adam algorithm */
  vec_zero(m);  vec_zero(m_prev);
  vec_zero(v);  vec_zero(v_prev);
  if (n_nuisance_params > 0) {
    vec_zero(m_nuis);  vec_zero(m_nuis_prev);
    vec_zero(v_nuis);  vec_zero(v_nuis_prev);
  }
  t = 0;
  
  do {
    /* we can precompute the KLD because it does not depend on the data under this model */
    /* (see equation 7, Doersch arXiv 2016) */
    trace = mmvn_trace(mmvn);  /* we'll reuse this */
    logdet = mmvn_log_det(mmvn);
    
    kld = 0.5 * (trace + mmvn_mu2(mmvn) - fulld - logdet);

    kld *= data->kld_upweight/(data->pointscale*data->pointscale); 
    
    /* we can also precompute the contribution of the KLD to the gradient */
    /* Note KLD is subtracted rather than added, so compute the gradient of -KLD */
    for (j = 0; j < grad->size; j++) {
      double gj = 0.0;

      if (j < n*dim)  /* partial deriv wrt mu_j is just mu_j */
        gj = -1.0*mmvn_get_mu_el(mmvn, j);
      else {            /* partial deriv wrt sigma_j is more
                           complicated because of the trace and log
                           determinant */
        if (data->type == CONST || data->type == DIST)
          gj = 0.5 * (fulld - trace);
        else if (data->type == DIAG) 
          gj = 0.5 * (1.0 - mat_get(mmvn->mvn->sigma, j-fulld, j-fulld)); 
        else 
          continue; /* LOWR case is messy; handle below */
      }
      vec_set(kldgrad, j, gj);
    }
    
    if (data->type == LOWR) {
      nj_set_kld_grad_LOWR(kldgrad, mmvn);

      if (data->sparsity != -1) { /* can also precompute */
        nj_set_sparsity_penalty_LOWR(sparsitygrad, mmvn, data);
        penalty = data->penalty;
      }
    }

    vec_scale(kldgrad, data->kld_upweight/(data->pointscale*data->pointscale));
    
    /* now sample a minibatch from the MVN averaging distribution and
       compute log likelihoods and gradients */
    vec_zero(avegrad);
    if (n_nuisance_params > 0)
      vec_zero(ave_nuis_grad);
    avell = 0;
    ave_nf_logdet = 0;
    for (i = 0; i < nminibatch; i++) {
      int bail = 0;
      double nf_logdet = 0; /* log det of Jacobian for normalizing flows;
                               will only be set if needed */
      do {
        /* do this in a way that keeps track of the original standard
           normal variable (points_std) for use in computing
           gradients.  Also use antithetic sampling to reduce
           variance */
        nj_sample_points(mmvn, data, points, points_std, &nf_logdet);
 
        ll = nj_compute_model_grad(mod, mmvn, points, points_std, grad, data);

        if (++bail > 10 && !isfinite(ll)) {
          fprintf(stderr, "WARNING: repeatedly sampling zero-probability trees. Prohibiting zero-length branches.\n");
          data->no_zero_br = TRUE;
          assert(bail < 15); /* prohibit infinite loop */
        }
      } while (!isfinite(ll));  /* in certain cases under the
                                   irreversible CRISPR model, trees
                                   can have likelihoods of zero; we'll
                                   try to just keep sampling if that
                                   happens and if necessary constrain
                                   the tree structure */
      
      avell += ll;
      vec_plus_eq(avegrad, grad);
      ave_nf_logdet += nf_logdet;

      if (n_nuisance_params > 0) {
        nj_update_nuis_grad(mod, data, nuis_grad);
        vec_plus_eq(ave_nuis_grad, nuis_grad);
      }
    }

    /* divide by nminibatch to get expected gradient */
    vec_scale(avegrad, 1.0/nminibatch);
    avell /= nminibatch;
    vec_plus_eq(avegrad, kldgrad);
    ave_nf_logdet /= nminibatch;
    if (data->type == LOWR && data->sparsity != -1)
      vec_plus_eq(avegrad, sparsitygrad);

    if (n_nuisance_params > 0) 
      vec_scale(ave_nuis_grad, 1.0/nminibatch);
    
    /* store parameters if best yet */
    if (avell - kld + ave_nf_logdet - penalty > bestelb) {
      bestelb = avell - kld + ave_nf_logdet - penalty;
      bestll = avell;  /* not necessarily best ll but ll corresponding to bestelb */
      bestkld = kld;  /* same comment */
      bestpenalty = penalty; 
      best_nf_logdet = ave_nf_logdet;
      bestt = t;
      mmvn_save_mu(mmvn, best_mu);
      vec_copy(best_sigmapar, sigmapar);
      if (n_nuisance_params > 0)
        nj_save_nuis_params(best_nuis_params, mod, data);
    }

    /* rescale gradient by approximate inverse Fisher information to
       put on similar scales; seems to help with optimization */
    if (data->natural_grad == TRUE)
      nj_rescale_grad(avegrad, rescaledgrad, mmvn, data);
    else
      vec_copy(rescaledgrad, avegrad);
    /* we won't do this with nuisance params */
    
    /* Adam updates; see Kingma & Ba, arxiv 2014 */
    t++;
    for (j = 0; j < rescaledgrad->size; j++) {   
      double mhatj, vhatj, g = vec_get(rescaledgrad, j);
      
      vec_set(m, j, ADAM_BETA1 * vec_get(m_prev, j) + (1.0 - ADAM_BETA1) * g);
      vec_set(v, j, ADAM_BETA2 * vec_get(v_prev, j) + (1.0 - ADAM_BETA2) * pow(g,2));
      mhatj = vec_get(m, j) / (1.0 - pow(ADAM_BETA1, t));
      vhatj = vec_get(v, j) / (1.0 - pow(ADAM_BETA2, t));

      /* update mu or sigma, depending on parameter index */
      if (j < fulld)
        mmvn_set_mu_el(mmvn, j, mmvn_get_mu_el(mmvn, j) +
                               learnrate * mhatj / (sqrt(vhatj) + ADAM_EPS)); 
      else 
        vec_set(sigmapar, j-fulld, vec_get(sigmapar, j-fulld) +
                learnrate * mhatj / (sqrt(vhatj) + ADAM_EPS)); 
    }
    nj_update_covariance(mmvn, data);
    
    vec_copy(m_prev, m);
    vec_copy(v_prev, v);

    /* same thing for nuisance params, if necessary */
    for (j = 0; j < n_nuisance_params; j++) {   
      double mhatj_nuis, vhatj_nuis, g = vec_get(ave_nuis_grad, j);
      vec_set(m_nuis, j, ADAM_BETA1 * vec_get(m_nuis_prev, j) + (1.0 - ADAM_BETA1) * g);
      vec_set(v_nuis, j, ADAM_BETA2 * vec_get(v_nuis_prev, j) + (1.0 - ADAM_BETA2) * pow(g,2));
      mhatj_nuis = vec_get(m_nuis, j) / (1.0 - pow(ADAM_BETA1, t));
      vhatj_nuis = vec_get(v_nuis, j) / (1.0 - pow(ADAM_BETA2, t));
      nj_nuis_param_pluseq(mod, data, j, learnrate * mhatj_nuis / (sqrt(vhatj_nuis) + ADAM_EPS));
    }
    if (n_nuisance_params > 0) {
      vec_copy(m_nuis_prev, m_nuis);
      vec_copy(v_nuis_prev, v_nuis);
    }
    
    /* report to log file */
    if (logf != NULL) {
      fprintf(logf, "%d\t%f\t%f\t%f\t%f\t", t, avell, kld, ave_nf_logdet, avell - kld + ave_nf_logdet - penalty);
      if (data->type == LOWR && data->sparsity != -1)
        fprintf(logf, "%f\t", data->penalty);
      mmvn_print(mmvn, logf, TRUE, FALSE);
      for (j = 0; j < sigmapar->size; j++)
        fprintf(logf, "%f\t", vec_get(sigmapar, j));
      for (j = 0; j < n_nuisance_params; j++)
        fprintf(logf, "%f\t", nj_nuis_param_get(mod, data, j)); 
      fprintf(logf, "\n");
    }
    if (gradf != NULL) {
      //      vec_print(avegrad, gradf);
      //vec_print(kldgrad, gradf);
      fprintf(gradf, "%d\t", t);
      for (j = 0; j < rescaledgrad->size; j++)
        fprintf(gradf, "%f\t", vec_get(rescaledgrad, j));
      fprintf(gradf, "\n");
    }
    
    /* check total elb every nbatches_conv to decide whether to stop */
    running_tot += avell - kld + ave_nf_logdet - penalty;
    if (t % nbatches_conv == 0) {
      if (logf != NULL)
        fprintf(logf, "# Average ELBO for last %d: %f\n", nbatches_conv, running_tot/nbatches_conv);
      if (t >= min_nbatches && 1.001*running_tot <= last_running_tot*0.999)
        /* sometimes get stuck increasingly asymptotically; stop if increase not more than about 0.1% */
        stop = TRUE;
      else {
        last_running_tot = running_tot;
        running_tot = 0;
      }
    }    
  } while(stop == FALSE);

  mmvn_set_mu(mmvn, best_mu);
  vec_copy(sigmapar, best_sigmapar);
  nj_update_covariance(mmvn, data);
  if (n_nuisance_params > 0)
    nj_update_nuis_params(best_nuis_params, mod, data);
  
  if (logf != NULL) {
    fprintf(logf, "# Reverting to parameters from iteration %d; ELB: %.2f, LNL: %.2f, KLD: %.2f, RFLD: %.2f",
            bestt+1, bestelb, bestll, bestkld, best_nf_logdet);
    if (data->type == LOWR && data->sparsity != -1)
      fprintf(logf, ", penalty: %f", bestpenalty);
    fprintf(logf, "\n");
  }
  
  vec_free(grad);
  vec_free(avegrad);
  vec_free(rescaledgrad);
  vec_free(kldgrad);
  if (data->type == LOWR && data->sparsity != -1)
    vec_free(sparsitygrad);
  vec_free(points);
  vec_free(points_std);
  vec_free(m);
  vec_free(m_prev);
  vec_free(v);
  vec_free(v_prev);
  vec_free(best_mu);
  vec_free(best_sigmapar);

  if (n_nuisance_params > 0) {
    vec_free(nuis_grad);
    vec_free(ave_nuis_grad);
    vec_free(m_nuis);
    vec_free(v_nuis);
    vec_free(m_nuis_prev);
    vec_free(v_nuis_prev);
    vec_free(best_nuis_params);
  }
    
}

/* sample a list of trees from the approximate posterior distribution
   and return as a new list.  If logdens is non-null, return
   corresponding vector of log densities for the samples */
List *nj_var_sample(int nsamples, multi_MVN *mmvn, CovarData *data, char** names,
                    Vector *logdens) {
  List *retval = lst_new_ptr(nsamples);
  int i;
  TreeNode *tree;
  Vector *points = vec_new(mmvn->d * mmvn->n);
  
  for (i = 0; i < nsamples; i++) {
    nj_sample_points(mmvn, data, points, NULL, NULL);

    if (logdens != NULL) 
      vec_set(logdens, i, mmvn_log_dens(mmvn, points));
     
    nj_points_to_distances(points, data);
    tree = nj_inf(data->dist, names, NULL, data);
    lst_push_ptr(retval, tree);
  }
  
  vec_free(points);
  return(retval);
}

/* return a single tree representing the approximate posterior mean */
TreeNode *nj_mean(Vector *mu, char **names, CovarData *data) {
  TreeNode *tree;
  
  if (data->nseqs * data->dim != mu->size)
    die("ERROR in nj_mean: bad dimensions\n");

  nj_points_to_distances(mu, data);  
  tree = nj_inf(data->dist, names, NULL, data);
  
  return(tree);
}

/* fully reset a tree model for use in likelihood calculations with a
   new tree */
void nj_reset_tree_model(TreeModel *mod, TreeNode *newtree) {
  if (mod->tree != NULL)
    tr_free(mod->tree);

  mod->tree = newtree;
  
  /* the size of the tree won't change, so we don't need to do a full
     call to tm_reset_tree; we also don't need to call
     tm_set_subst_matrices because it will be called each time the
     likelihood function is invoked */
}

/* generate an approximate multivariate normal distribution from a
   distance matrix, for use in initializing the variational inference
   algorithm.  */
void nj_estimate_mmvn_from_distances(CovarData *data, multi_MVN *mmvn) {
  if (data->hyperbolic)
    nj_estimate_mmvn_from_distances_hyperbolic(data, mmvn);
  else
    nj_estimate_mmvn_from_distances_euclidean(data, mmvn);  
}

/* generate an approximate multivariate normal distribution from a distance matrix, for
   use in initializing the variational inference algorithm. Uses multidimensional scaling  */
void nj_estimate_mmvn_from_distances_euclidean(CovarData *data, multi_MVN *mmvn) {
  Matrix *D = data->dist;
  int n = D->nrows;
  Matrix *Dsq, *G, *revec_real;
  Vector *eval_real;
  int i, j, d;
  Vector *mu_full = vec_new(data->dim * n);
  
  if (D->nrows != D->ncols || mmvn->d * mmvn->n != data->dim * n)
    die("ERROR in nj_estimate_points_from_distances: bad dimensions\n");

  /* build matrix of squared distances; note that D is upper
     triangular but Dsq must be symmetric */
  Dsq = mat_new(n, n);
  for (i = 0; i < n; i++) {
    mat_set(Dsq, i, i, 0);
    for (j = i + 1; j < n; j++) {
      double d2 = mat_get(D, i, j) * mat_get(D, i, j);
      mat_set(Dsq, i, j, d2);
      mat_set(Dsq, j, i, d2);
    }
  }

  /* double center */
  G = mat_new(n, n);
  mat_double_center(G, Dsq, FALSE);
  
  /* find eigendecomposition of G */
  eval_real = vec_new(n); vec_zero(eval_real);
  revec_real = mat_new(n, n); mat_zero(revec_real);
  if (mat_diagonalize_sym(G, eval_real, revec_real) != 0)
    die("ERROR in nj_estimate_mmvn_from_distances_euclidean: diagonalization failed.\n");
  
  /* create a vector of points based on the first 'dim' eigenvalues */
  for (d = 0; d < data->dim; d++) {
    double eval = vec_get(eval_real, n-1-d);
    if (eval < 0) eval = 0;
    double sqeval = sqrt(eval);
    /* product of evalsqrt and corresponding column of revec will define
       the dth component of each point */    
    for (i = 0; i < n; i++) {
      vec_set(mu_full, i*data->dim + d,
              sqeval * mat_get(revec_real, i, n-1-d));
    }
  }
  
  /* rescale */
  vec_scale(mu_full, data->pointscale);
  mmvn_set_mu(mmvn, mu_full);

  /* covariance parameters should already be initialized */
  nj_update_covariance(mmvn, data);
  
  mat_free(Dsq);
  mat_free(G);
  vec_free(eval_real);
  mat_free(revec_real);
  vec_free(mu_full);
}

/* generate an approximate mu and sigma from a distance matrix, for
   use in initializing the variational inference algorithm. In this
   version, use the 'hydra' algorithm to solve the problem
   approximately in hyperbolic space (Keller-Ressel & Nargang,
   arXiv:1903.08977, 2019) */
void nj_estimate_mmvn_from_distances_hyperbolic(CovarData *data, multi_MVN *mmvn) {
  Matrix *D = data->dist;
  int n = D->nrows;
  Matrix *A, *revec_real;
  Vector *eval_real;
  int i, j, d;
  Vector *mu_full = vec_new(data->dim*n);
    
  if (D->nrows != D->ncols || mmvn->d * mmvn->n != data->dim * n)
    die("ERROR in nj_estimate_points_from_distances_hyperbolic: bad dimensions\n");
  
  /* build matrix A of transformed distances; note that D is upper
     triangular but A must be symmetric */
  A = mat_new(n, n);
  for (i = 0; i < n; i++) {
    mat_set(A, i, i, 1);
    for (j = i + 1; j < n; j++) {
      double a = cosh(sqrt(data->negcurvature) * mat_get(D, i, j) * data->pointscale);
      mat_set(A, i, j, a);
      mat_set(A, j, i, a);
    }
  }
  
  /* find eigendecomposition of A */
  eval_real = vec_new(n);
  revec_real = mat_new(n, n);
  if (mat_diagonalize_sym(A, eval_real, revec_real) != 0)
    die("ERROR in nj_estimate_mmvn_from_distances_hyperbolic: diagonalization failed.\n");

  /* create a vector of points based on the first 'dim' eigenvalues */
  for (d = 0; d < data->dim; d++) {
    /* product of evalsqrt and corresponding column of revec will define
       the dth component of each point */
    double ev = -vec_get(eval_real, d);
    assert(isfinite(ev));
    if (ev < 0) ev = 1e-6;
    for (i = 0; i < n; i++) {
      vec_set(mu_full, i*data->dim + d, sqrt(ev) * mat_get(revec_real, i, d)); 
      assert(isfinite(vec_get(mu_full, i*data->dim + d)));
    }
  }

  mmvn_set_mu(mmvn, mu_full); 
  
  /* covariance parameters should already be initialized */
  nj_update_covariance(mmvn, data);
  
  mat_free(A);
  vec_free(eval_real);
  mat_free(revec_real);
  vec_free(mu_full);
}

/* ensure a distance matrix is square, upper triangular, has zeroes on
   main diagonal, and has all non-negative entries.  */
void nj_test_D(Matrix *D) {
  int n = D->nrows;
  int i, j;
  if (n != D->ncols)
    die("ERROR in nj_test_D: bad dimensions in distance matrix D\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j <= i && mat_get(D, i, j) != 0)
	die("ERROR in nj_test_D: distance matrix must be upper triangular and have zeroes on main diagonal.\n");
      else if (mat_get(D, i, j) < 0 || !isfinite(mat_get(D, i, j)))
	die("ERROR in nj_test_D: entries in distance matrix must be nonnegative and finite\n");
    }
  }
}

/* Compute and return the log likelihood of a tree model with respect
   to an alignment.  This is a pared-down version of
   tl_compute_log_likelihood.  It assumes a 0th order model,
   sufficient stats already available, and no rate variation.  It also
   makes use of scaling factors to avoid underflow with large trees.  If
   branchgrad is non-null, it will be populated with the gradient of
   the log likelihood with respect to the individual branches of the
   tree, in post-order.   */
double nj_compute_log_likelihood(TreeModel *mod, CovarData *data, Vector *branchgrad) {

  int i, j, k, nodeidx, rcat = 0, tupleidx;
  int nstates = mod->rate_matrix->size;
  TreeNode *n;
  double total_prob = 0;
  List *traversal;
  double **pL = NULL, **pLbar = NULL;
  Vector *lscale, *lscale_o; /* inside and outside versions */
  double scaling_threshold = DBL_MIN * 1.0e10;  /* need some padding */
  double lscaling_threshold = log(scaling_threshold), ll = 0;
  double tmp[nstates];
  Matrix **grad_mat, **grad_mat_kappa;
  unsigned int rescale;
  MSA *msa = data->msa;
  
  if (msa->ss->tuple_size != 1)
    die("ERROR nj_compute_log_likelihood: need tuple size 1, got %i\n",
	msa->ss->tuple_size);
  if (mod->order != 0)
    die("ERROR nj_compute_log_likelihood: got mod->order of %i, expected 0\n",
	mod->order);
  if (!mod->allow_gaps)
    die("ERROR nj_compute_log_likelihood: need mod->allow_gaps to be TRUE\n");
  if (mod->nratecats > 1)
    die("ERROR nj_compute_log_likelihood: no rate variation allowed\n");
  
  pL = smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++)
    pL[j] = smalloc((mod->tree->nnodes+1) * sizeof(double));

  /* we also need to keep track of the log scale of every node for
     underflow purposes */
  lscale = vec_new(mod->tree->nnodes+1); 
  lscale_o = vec_new(mod->tree->nnodes+1);
  
  if (branchgrad != NULL) {
    if (branchgrad->size != mod->tree->nnodes-1) /* rooted */
      die("ERROR in nj_compute_log_likelihood: size of branchgrad must be 2n-2\n");
    vec_zero(branchgrad);
    pLbar = smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++)
      pLbar[j] = smalloc((mod->tree->nnodes+1) * sizeof(double));
    if (mod->subst_mod == HKY85)
      data->deriv_hky_kappa = 0.0;
    grad_mat = malloc(mod->tree->nnodes * sizeof(void*));
    for (j = 0; j < mod->tree->nnodes; j++)
      grad_mat[j] = mat_new(nstates, nstates);
    if (mod->subst_mod == HKY85) {
      grad_mat_kappa = malloc(mod->tree->nnodes * sizeof(void*));
    for (j = 0; j < mod->tree->nnodes; j++)
      grad_mat_kappa[j] = mat_new(nstates, nstates);
    }
  }

  tm_set_subst_matrices(mod);  /* just call this in all cases; we'll be tweaking the model a lot */
  
  /* get sequence index if not already there */
  if (mod->msa_seq_idx == NULL)
    tm_build_seq_idx(mod, msa);

  for (tupleidx = 0; tupleidx < msa->ss->ntuples; tupleidx++) {
    /* reset scale */
    vec_zero(lscale); vec_zero(lscale_o);
    
    traversal = tr_postorder(mod->tree);
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);

      if (n->lchild == NULL) {
        /* leaf: base case of recursion */
        int state = mod->rate_matrix->
          inv_states[(int)ss_get_char_tuple(msa, tupleidx,
                                            mod->msa_seq_idx[n->id], 0)];
        for (i = 0; i < nstates; i++) {
          if (state < 0 || i == state)
            pL[i][n->id] = 1;
          else
            pL[i][n->id] = 0;
        }
      }
      else {
        /* general recursive case */
        MarkovMatrix *lsubst_mat = mod->P[n->lchild->id][rcat];
        MarkovMatrix *rsubst_mat = mod->P[n->rchild->id][rcat];
        
        rescale = FALSE;
        for (i = 0; i < nstates; i++) {
          double totl = 0, totr = 0;
          for (j = 0; j < nstates; j++)
            totl += pL[j][n->lchild->id] *
              mm_get(lsubst_mat, i, j);
          
          for (k = 0; k < nstates; k++)
            totr += pL[k][n->rchild->id] *
              mm_get(rsubst_mat, i, k);

          pL[i][n->id] = totl * totr;
          if (totl > 0 && totr > 0 && pL[i][n->id] < scaling_threshold)
            rescale = TRUE;
        }

        /* deal with nodewise scaling */
        vec_set(lscale, n->id, vec_get(lscale, n->lchild->id) +
                vec_get(lscale, n->rchild->id));
        if (rescale == TRUE) { /* have to rescale for all states */
          vec_set(lscale, n->id, vec_get(lscale, n->id) + lscaling_threshold);
          for (i = 0; i < nstates; i++) 
            pL[i][n->id] /= scaling_threshold;
        }
      }
    }
  
    /* termination */
    total_prob = 0;
    for (i = 0; i < nstates; i++)
      total_prob += vec_get(mod->backgd_freqs, i) *
        pL[i][mod->tree->id] * mod->freqK[rcat];
    
    ll += (log(total_prob) + vec_get(lscale, mod->tree->id)) * msa->ss->counts[tupleidx];
    assert(isfinite(ll));

    /* to compute gradients efficiently, need to make a second pass
       across the tree to compute "outside" probabilities */
    if (branchgrad != NULL) {
      double expon;
      traversal = tr_preorder(mod->tree);

      for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
        n = lst_get_ptr(traversal, nodeidx);

        if (n->parent == NULL) { /* base case */
          for (i = 0; i < nstates; i++)
            pLbar[i][n->id] = vec_get(mod->backgd_freqs, i);
        }
        else {            /* recursive case */
          TreeNode *sibling = (n == n->parent->lchild ?
                               n->parent->rchild : n->parent->lchild);
          MarkovMatrix *par_subst_mat = mod->P[n->id][rcat];
          MarkovMatrix *sib_subst_mat = mod->P[sibling->id][rcat];
          rescale = FALSE;
          
          for (j = 0; j < nstates; j++) { /* parent state */
            tmp[j] = 0;
            for (k = 0; k < nstates; k++)  /* sibling state */
              tmp[j] += pLbar[j][n->parent->id] *
                pL[k][sibling->id] * mm_get(sib_subst_mat, j, k);
          }
          
          for (i = 0; i < nstates; i++) { /* child state */
            pLbar[i][n->id] = 0;
            for (j = 0; j < nstates; j++) { /* parent state */
              pLbar[i][n->id] +=
                tmp[j] * mm_get(par_subst_mat, j, i);
              if (tmp[j] > 0 && mm_get(par_subst_mat, j, i) > 0 &&
                  pLbar[i][n->id] < scaling_threshold)
                rescale = TRUE;
            }
          }
          vec_set(lscale_o, n->id, vec_get(lscale_o, n->parent->id) +
                  vec_get(lscale, sibling->id));
          if (rescale == TRUE) { /* rescale for all states */
            vec_set(lscale_o, n->id, vec_get(lscale_o, n->id) + lscaling_threshold);
            for (i = 0; i < nstates; i++)
              pLbar[i][n->id] /= scaling_threshold;     
          }
        }
      }

      /* TEMPORARY: check inside/outside */
      /* for (nodeidx = 0; nodeidx < lst_size(mod->tree->nodes); nodeidx++) { */
      /*   double pr = 0; */
      /*   n = lst_get_ptr(mod->tree->nodes, nodeidx); */
      /*   assert(vec_get(lscale, n->id) == 0 && vec_get(lscale_o, n->id) == 0); */
      /*   for (j = 0; j < nstates; j++) */
      /*     pr += pL[j][n->id] * pLbar[j][n->id]; */
      /*   printf("Tuple %d, node %d: %f (%f)\n", tupleidx, nodeidx, log(pr), log(total_prob)); */
      /* } */

      
      /* now compute branchwise derivatives in a final pass */
      traversal = mod->tree->nodes;
      for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
        TreeNode *par, *sib;
        double base_prob = total_prob, deriv;
        
        n = lst_get_ptr(traversal, nodeidx);
        par = n->parent;
	
        if (par == NULL) 
          continue;
       
        sib = (n == n->parent->lchild ?
               n->parent->rchild : n->parent->lchild);

        /* this part is just a constant to propagate through to the
           derivative */
        for (i = 0; i < nstates; i++) {  /* parent */
          tmp[i] = 0;
          for (k = 0; k < nstates; k++)  /* sibling */
            tmp[i] += pL[k][sib->id] * mm_get(mod->P[sib->id][rcat], i, k);
        }

        if (n != mod->tree->rchild) { /* skip branch to right of root because unrooted */
          /* calculate derivative analytically */
          deriv = 0;
          /* only do this first time through */
          if (tupleidx == 0) {
            if (mod->subst_mod == JC69)
              tm_grad_JC69(mod, grad_mat[n->id], n->dparent);
            else if (mod->subst_mod == HKY85) 
              tm_grad_HKY_dt(mod, grad_mat[n->id], data->hky_kappa, n->dparent); 
            else
              die("ERROR in nj_compute_log_likelihood: only JC69 and HKY85 substitution models are supported.\n");
          }
          
          for (i = 0; i < nstates; i++)   
            for (j = 0; j < nstates; j++)    
              deriv +=  tmp[i] * pLbar[i][par->id] * pL[j][n->id] * mat_get(grad_mat[n->id], i, j);

          /* adjust for all relevant scale terms; do everything in log space */
          expon = -vec_get(lscale, mod->tree->id)
            + vec_get(lscale, sib->id) + vec_get(lscale_o, par->id)
            + vec_get(lscale, n->id) - log(base_prob);
          /* note division by base_prob because we need deriv of log P */

          /* avoid overflow */
          if (expon > 700.0) expon = 700.0;
          if (expon < -745.0) expon = -745.0;
          
          deriv *= exp(expon);
          assert(isfinite(deriv));
                  
          vec_set(branchgrad, n->id, vec_get(branchgrad, n->id) +
                  deriv * msa->ss->counts[tupleidx]);
        }

        /* in this case, we need a partial derivative for kappa also;
           it has to be aggregated across all branches */
        if (mod->subst_mod == HKY85) {
          double this_deriv_kappa = 0;
          if (tupleidx == 0)
            tm_grad_HKY_dkappa(mod, grad_mat_kappa[n->id], data->hky_kappa, n->dparent);
          for (i = 0; i < nstates; i++) 
            for (j = 0; j < nstates; j++) 
              this_deriv_kappa += tmp[i] * pLbar[i][par->id] * pL[j][n->id] *
                mat_get(grad_mat_kappa[n->id], i, j);

          /* adjust for all relevant scale terms */
          this_deriv_kappa *= exp(expon);        
          data->deriv_hky_kappa += (this_deriv_kappa * msa->ss->counts[tupleidx]);        
        }
      }
    }
  }
  
  for (j = 0; j < nstates; j++)
    sfree(pL[j]);
  sfree(pL);
  
  if (branchgrad != NULL) {
    for (j = 0; j < nstates; j++)
      sfree(pLbar[j]);
    sfree(pLbar);
    for (j = 0; j < mod->tree->nnodes; j++)      
      mat_free(grad_mat[j]);
    free(grad_mat);
    if (mod->subst_mod == HKY85) {
      for (j = 0; j < mod->tree->nnodes; j++)      
        mat_free(grad_mat_kappa[j]);
      free(grad_mat_kappa);
    }
  }

  vec_free(lscale);
  vec_free(lscale_o);
  
  return ll;
}

/* Build index of leaf ids to sequence indices based on a name
   list. */
int *nj_build_seq_idx(List *leaves, char **names) {
  int i;  
  int *retval = smalloc(lst_size(leaves)*2 * sizeof(int));
  for (i = 0; i < lst_size(leaves)*2; i++) retval[i] = -1;
  for (i = 0; i < lst_size(leaves); i++) {
    TreeNode *n = lst_get_ptr(leaves, i);
    retval[n->id] = nj_get_seq_idx(names, n->name, lst_size(leaves));
    if (retval[n->id] < 0)
      die("ERROR: leaf '%s' not found in name list.\n", n->name);
  }
  return retval;
}

/* Return index of given sequence name or -1 if not found. */
int nj_get_seq_idx(char **names, char *name, int n) {
  int i, retval = -1;
  for (i = 0; retval < 0 && i < n; i++) 
    if (!strcmp(name, names[i]))
      retval = i;
  return retval;
}

/* subsample from a set of trees by importance sampling, using ratio
   of likelihoods to sampling density as weights.  Warning: tree
   objects in returned list will be shared with those in primary list
   and may repeat */
List *nj_importance_sample(int nsamples, List *trees, Vector *logdens,
                           TreeModel *mod, CovarData *data, FILE *logf) {
  List *retval = lst_new_ptr(nsamples);
  Vector *weights = vec_new(lst_size(trees));
  Vector *lls = vec_new(lst_size(trees));
  double ll, maxll = -INFTY, maxweight = -INFTY, sampll = 0;
  int i, count = 0;

  if (logdens->size != lst_size(trees))
    die("ERROR in nj_importance_sample: bad input.\n");
  
  /* calculate importance weights from likelihoods */
  for (i = 0; i < lst_size(trees); i++) {
    TreeNode *t = lst_get_ptr(trees, i);    
    mod->tree = t;
    ll = nj_compute_log_likelihood(mod, data, NULL);
    if (!isfinite(ll))  /* can happen with crispr model */
      ll = -INFTY;
    if (ll > maxll) maxll = ll;
    vec_set(lls, i, ll);
    ll -= vec_get(logdens, i);
    vec_set(weights, i, ll);
    if (ll > maxweight) maxweight = ll;
  }
  assert(maxweight > -INFTY);

  /* exponentiate and renormalize */
  for (i = 0; i < lst_size(trees); i++) {
    ll = vec_get(weights,i);
    vec_set(weights, i, exp(ll-maxweight));  /* avoids underflow */
  }
  pv_normalize(weights);

  for (i = 0; i < lst_size(trees); i++) 
    if (vec_get(weights, i) > 1.0/(lst_size(trees)*10))
      count++;

  if (count < nsamples)
    fprintf(stderr, "Warning: only %d trees eligible for importance sampling...\n", count);
    
  /* now draw nsamples */
  for (i = 0; i < nsamples; i++) {
    int j = pv_draw_idx(weights);
    assert(j >= 0 && j < lst_size(trees));
    lst_push_ptr(retval, lst_get_ptr(trees, j));
    sampll += vec_get(lls, j);
  }

  if (logf != NULL)
    fprintf(logf, "# Importance sampling from %d eligible trees of %d; avelnl: %f, maxlnl: %f\n",
            count, lst_size(trees), sampll/nsamples, maxll);

  vec_free(weights);
  vec_free(lls);
  
  return(retval);
}

/* update covariance matrix based on the parameters and auxiliary
   data */
void nj_update_covariance(multi_MVN *mmvn, CovarData *data) {
  int i, j;
  Vector *sigma_params = data->params;
  
  /* Note: variance parameters now stored as log values and must
     be exponentiated */
  mat_zero(mmvn->mvn->sigma);
  if (data->type == CONST) {
    mat_set_identity(mmvn->mvn->sigma);
    data->lambda = (VARFLOOR + exp(vec_get(sigma_params, 0)));
    if (!isfinite(data->lambda)) data->lambda = VARFLOOR;
    mat_scale(mmvn->mvn->sigma, data->lambda);
  }
  else if (data->type == DIAG) {
    for (i = 0; i < sigma_params->size; i++) 
      mat_set(mmvn->mvn->sigma, i, i, VARFLOOR + exp(vec_get(sigma_params, i)));
  }
  else if (data->type == DIST) {
    data->lambda = VARFLOOR + exp(vec_get(sigma_params, 0));
    mat_copy(mmvn->mvn->sigma, data->Lapl_pinv);
    mat_scale(mmvn->mvn->sigma, data->lambda);

    if (mmvn->mvn->evecs == NULL) {
      mmvn->mvn->evecs = mat_new(mmvn->n, mmvn->n);
      mmvn->mvn->evals = vec_new(mmvn->n);
    }
    mat_copy(mmvn->mvn->evecs, data->Lapl_pinv_evecs); /* can simply derive eigendecomposition 
                                                          from Lapl_pinv */
    vec_copy(mmvn->mvn->evals, data->Lapl_pinv_evals);
    vec_scale(mmvn->mvn->evals, data->lambda);
  }
  else {
    assert(data->type == LOWR);
    for (i = 0; i < data->R->nrows; i++)
      for (j = 0; j < data->R->ncols; j++)
        mat_set(data->R, i, j, vec_get(data->params,
                                       i*data->R->ncols + j));
                                 /* note not log in this case */

    /* update the mvn accordingly */
    mvn_reset_LOWR(mmvn->mvn, data->R);
  }
}

/* create a new CovarData object appropriate for the choice of
   parameterization */
CovarData *nj_new_covar_data(enum covar_type covar_param, Matrix *dist, int dim,
                             MSA *msa, CrisprMutModel* crispr_mod, char **names,
                             unsigned int natural_grad, double kld_upweight,
                             int rank, double sparsity, unsigned int hyperbolic,
                             double negcurvature, unsigned int ultrametric,
                             unsigned int radial_flow, unsigned int planar_flow) {
  static int seeded = 0;
  
  CovarData *retval = smalloc(sizeof(CovarData));
  retval->type = covar_param;
  retval->msa = msa;
  retval->crispr_mod = crispr_mod;
  retval->names = names;
  retval->lambda = LAMBDA_INIT;
  retval->mvn_type = MVN_DIAG;
  retval->dist = dist;
  retval->nseqs = dist->nrows;
  retval->dim = dim;
  retval->natural_grad = natural_grad;
  retval->kld_upweight = kld_upweight;
  retval->Lapl_pinv = NULL;
  retval->Lapl_pinv_evals = NULL;
  retval->Lapl_pinv_evecs = NULL;
  retval->lowrank = -1;
  retval->R = NULL;
  retval->sparsity = sparsity;
  retval->hyperbolic = hyperbolic;
  retval->negcurvature = negcurvature;
  retval->ultrametric = ultrametric;
  retval->no_zero_br = FALSE;

  if (radial_flow == TRUE) {
    retval->rf = rf_new(retval->nseqs, dim);
    rf_rescale(retval->rf, POINTSPAN_EUC/sqrt(2));
  }
  else
    retval->rf = NULL;

  if (planar_flow == TRUE) 
    retval->pf = pf_new(retval->nseqs, dim);
  else
    retval->pf = NULL;
    
  nj_set_pointscale(retval);
  retval->lambda *= retval->pointscale * retval->pointscale;
  
  if (covar_param == CONST) {
    /* store constant */
    retval->params = vec_new(1);
    vec_set(retval->params, 0, log(max(retval->lambda-VARFLOOR, VARFLOOR)));  /* use lambda for scale; log parameterization */
  }
  else if (covar_param == DIAG) {
    int i, j;
    double x = 0, x2 = 0, N = retval->nseqs * (retval->nseqs-1)/2.0;
    retval->params = vec_new(retval->dim * retval->nseqs);
    /* initialize sigma parameters to (log of) identity scaled by
       1/n of the variance across pairwise distances */
    for (i = 0; i < dist->nrows; i++) {
      for (j = i+1; j < dist->ncols; j++) {
        x += mat_get(dist, i, j);
        x2 += mat_get(dist, i, j) * mat_get(dist, i, j);
      }
    }
    vec_set_all(retval->params, log(1.0/(N*N) * (x2/N - x*x/(N*N))) + 2*log(retval->pointscale));
    /* note scale by (log of) pointscale^2 */
  }  
  else if (covar_param == DIST) {
    retval->mvn_type = MVN_GEN;
    retval->params = vec_new(1);
    vec_set(retval->params, 0, log(max(retval->lambda-VARFLOOR, VARFLOOR)));
    retval->Lapl_pinv = mat_new(dist->nrows, dist->ncols);
    retval->Lapl_pinv_evals = vec_new(dist->nrows);
    retval->Lapl_pinv_evecs = mat_new(dist->nrows, dist->nrows);
    nj_laplacian_pinv(retval);  /* set up the Laplacian pseudoinverse */    
  }
  else if (covar_param == LOWR) {
    double sdev;
    int i, j;
    
    retval->lowrank = rank;
    retval->mvn_type = MVN_LOWR;
    retval->params = vec_new(retval->lowrank * retval->nseqs);
    retval->R = mat_new(retval->nseqs, retval->lowrank);

    /* initialization is tricky; we want variances on the order of
       0.01 and expected covariances of 0 but we want to avoid
       orthogonality; initialize randomly with appropriate distrib */
    if (!seeded) {
      srandom((unsigned int)time(NULL));
      seeded = 1;
    }
    sdev = sqrt(LAMBDA_INIT / retval->lowrank); /* yields expected variance
                                                   of LAMBDA_INIT and expected
                                                   covariances of 0 */
    for (i = 0; i < retval->nseqs; i++) {
      for (j = 0; j < retval->lowrank; j++) {
        double draw = norm_draw(0, sdev);
        mat_set(retval->R, i, j, draw);
        vec_set(retval->params, i*retval->lowrank + j, draw);
      }
    }
  }
  else
    die("ERROR in nj_new_covar_data: unrecognized type.\n");

  return (retval);
}

void nj_dump_covar_data(CovarData *data, FILE *F) {
  fprintf(F, "CovarData\nnseqs: %d\ndim: %d\nlambda: %f\n", data->nseqs, data->dim, data->lambda);
  fprintf(F, "distance matrix:\n");
  mat_print(data->dist, F);
  fprintf(F, "Free parameters: ");
  vec_print(data->params, F);
  if (data->type == DIST) {
    fprintf(F, "Laplacian pseudoinverse:\n");
    mat_print(data->Lapl_pinv, F);
    fprintf(F, "Eigenvalues:\n");
    vec_print(data->Lapl_pinv_evals, F);
    fprintf(F, "Eigenvectors:\n");
    mat_print(data->Lapl_pinv_evecs, F);
  }
  else if (data->type == LOWR) {
    fprintf(F, "Low-rank matrix R (rank %d):\n", data->lowrank);
    mat_print(data->R, F);
  }
}

/* define Laplacian pseudoinverse from distance matrix, for use with
   DIST parameterization of covariance.  Also compute
   eigendecomposition for use in gradient calculations. Store
   everything in CovarData object */
void nj_laplacian_pinv(CovarData *data) {
  int i, dim = data->dist->nrows;
  double epsilon, trace;
  
  /* define Laplacian pseudoinverse as double centered version of
     distance matrix */
  mat_double_center(data->Lapl_pinv, data->dist, TRUE);
  
  /* this matrix is only defined up to a translation because it is
     based on pairwise distances.  For it to define a valid covariance
     matrix (up to a positive scale constant) we need to ensure that
     it is positive definite.  We can do this by adding epsilon * I to
     it, where epsilon is equal to the smallest eigenvalue plus a
     small margin.  This will preserve the eigenvectors but shift all
     eigenvalues upward by epsilon */
  if (mat_diagonalize_sym(data->Lapl_pinv, data->Lapl_pinv_evals, data->Lapl_pinv_evecs) != 0)
    die("ERROR in nj_laplacian_pinv: diagonalization failed.\n");    

  epsilon = vec_get(data->Lapl_pinv_evals, 0) + 1e-6;
  /* mat_diagonalize_sym guarantees eigenvalues are in ascending order */

  for (i = 0; i < dim; i++) {
    mat_set(data->Lapl_pinv, i, i, mat_get(data->Lapl_pinv, i, i) + epsilon);
    vec_set(data->Lapl_pinv_evals, i, vec_get(data->Lapl_pinv_evals, i) + epsilon);
  }

  /* finally, rescale so average diagonal element is one (in space of
     final sigma), putting matrix on the same scale as the identity
     matrix used for the CONST parameterization */
  trace = vec_sum(data->Lapl_pinv_evals);
  mat_scale(data->Lapl_pinv, data->nseqs / (trace * data->dim));
  vec_scale(data->Lapl_pinv_evals, data->nseqs / (trace * data->dim));
}

/* wrapper for nj_points_to_distances functions */
void nj_mmvn_to_distances(multi_MVN *mmvn, CovarData *data) {
  Vector *full_mu;
  if (mmvn->type != MVN_GEN && mmvn->type != MVN_LOWR) 
    full_mu = mmvn->mvn->mu;
  else {
    full_mu = vec_new(mmvn->d * mmvn->n);
    mmvn_save_mu(mmvn, full_mu);
  }

  nj_points_to_distances(full_mu, data);

  if (mmvn->type == MVN_GEN || mmvn->type == MVN_LOWR)
    vec_free(full_mu);
}

/* compute partial derivatives of KLD wrt variance parameters in LOWR
   case */
void nj_set_kld_grad_LOWR(Vector *kldgrad, multi_MVN *mmvn) {
  int i, j;
  int offset = mmvn->d * mmvn->n;
  Matrix *Rgrad = mat_new(mmvn->mvn->lowR->nrows, mmvn->mvn->lowR->ncols);

  /* calculate partial derivatives using matrix operations, making use
     of precomputed R^T x R */
  mat_mult(Rgrad, mmvn->mvn->lowR, mmvn->mvn->lowR_invRtR);
  mat_minus_eq(Rgrad, mmvn->mvn->lowR);
  mat_scale(Rgrad, mmvn->d);  /* note: computing negative gradient; that is what we need */

  /* populate vector from matrix */
  for (i = 0; i < mmvn->mvn->lowR->nrows; i++) 
    for (j = 0; j < mmvn->mvn->lowR->ncols; j++) 
      vec_set(kldgrad, offset + i*mmvn->mvn->lowR->ncols + j, mat_get(Rgrad, i, j));

  mat_free(Rgrad);
}

/* rescale gradients by approximate inverse Fisher information for approx
   natural gradient scale */
void nj_rescale_grad(Vector *grad, Vector *rsgrad, multi_MVN *mmvn, CovarData *data) {
  int i, j, fulld = mmvn->n * mmvn->d;
  for (i = 0; i < grad->size; i++) {
    double g = vec_get(grad, i);
    
    if (i < fulld) { /* mean gradients */
      if (data->type == CONST || data->type == DIAG) 
        g *= mat_get(mmvn->mvn->sigma, i, i); /* these are all the same in the CONST case */
      else { /* DIST or LOWR */ /* CHECK: this code untested */
        /* scale by dot product of corresponding row of sigma with the original gradient */
        double dotp = 0.0;
        int sigmarow = i / mmvn->d;  /* project down to sigma */
        int d = i % mmvn->d; /* corresponding dimension */
        for (j = 0; j < mmvn->mvn->sigma->ncols; j++)
          dotp += mat_get(mmvn->mvn->sigma, sigmarow, j) * vec_get(grad, j*mmvn->d + d);
        g *= dotp;
      }
    }
    else { /* variance gradients */
      if (data->type == CONST)
        g *= 2/(mmvn->n * mmvn->d); /* assumes variance is exp(parameter) */
      else if (data->type == DIAG)
        g *= 2; /* assumes variance is exp(parameter) */
      else if (data->type == DIST)
        g *= 2/(mmvn->n-1); /* assumes variance is exp(parameter) */
      else
        break; /* handle LOWR case below */
    }
    
    vec_set(rsgrad, i, g);
  }

  if (data->type == LOWR) { /* CHECK: this code as yet untested */
    /* in this case the rescaled variance gradients can be obtained by
       matrix multiplication with sigma */

    /* first coerce the relevant gradient components into an n x k matrix */
    Matrix *Rgrad = mat_new(mmvn->n, data->lowrank),
      *rsRgrad = mat_new(mmvn->n, data->lowrank);
    for (i = 0; i < mmvn->n; i++)
      for (j = 0; j < data->lowrank; j++)
        mat_set(Rgrad, i, j, vec_get(grad, fulld + i*data->lowrank + j));

    /* multiply on the left by sigma */
    mat_mult(rsRgrad, mmvn->mvn->sigma, Rgrad);

    /* finally extract the rescaled values */
    for (i = 0; i < mmvn->n; i++)
      for (j = 0; j < data->lowrank; j++)
        vec_set(rsgrad, fulld + i*data->lowrank + j, mat_get(rsRgrad, i, j));
    
    mat_free(Rgrad);
    mat_free(rsRgrad);
  }
}

 /* Compute sparsity penalty and its gradient in LOWR case */
void nj_set_sparsity_penalty_LOWR(Vector *grad, multi_MVN *mmvn,
                                  CovarData *data) {  
  double penalty = 0, scale = 0;
  Matrix *sigprime, *Rgrad = mat_new(data->R->nrows, data->R->ncols);
  int i, j, n = mmvn->n, d = mmvn->d;
  
  assert(data->type == LOWR && data->sparsity != -1);

  for (i = 0; i < n; i++) {
    for (j = i+1; j < n; j++)
      penalty += 2 * pow(mat_get(mmvn->mvn->sigma, i, j), 2);
    scale += pow(mat_get(mmvn->mvn->sigma, i, i), 2);
  }

  /* rescale the penalty factor */
  scale *= (n-1);   
  /* If all nondiagonal entries were like the diagonal ones, this
     would be the expected value of the raw penalty.  We will rescale
     it so this expected value is one. */
    
  data->penalty = data->sparsity/scale * penalty; 

  /* calculate gradient in matrix form */
  /* first make copy of sigma with 0s on main diagonal */
  sigprime = mat_create_copy(mmvn->mvn->sigma);
  for (i = 0; i < sigprime->nrows; i++)
    mat_set(sigprime, i, i, 0);
  
  /* multiply by 4R * sparsity for gradient */
  mat_mult(Rgrad, sigprime, data->R);
  mat_scale(Rgrad, 4.0 * data->sparsity/scale);
    
  /* finally add entries to corresponding gradient components */
  vec_zero(grad);
  for (i = 0; i < data->R->nrows; i++)
    for (j = 0; j < data->R->ncols; j++)
      vec_set(grad, n*d + i*data->R->ncols + j,
              vec_get(grad,  n*d + i*data->R->ncols + j) - 
              mat_get(Rgrad, i, j));
  /* note have to subtract because makes negative contribution */

  mat_free(sigprime);
  mat_free(Rgrad);
}

 /* Same as above but use an L1 (LASSO) penalty instead of an L2
    (ridge) penalty */
void nj_set_LASSO_penalty_LOWR(Vector *grad, multi_MVN *mmvn,
                               CovarData *data) {  
  double penalty = 0, scale = 0;
  Matrix *sign = mat_new(mmvn->mvn->sigma->nrows, mmvn->mvn->sigma->ncols),
    *Rgrad = mat_new(data->R->nrows, data->R->ncols);
  int i, j, n = mmvn->n, d = mmvn->d;
  
  assert(data->type == LOWR && data->sparsity != -1);

  for (i = 0; i < n; i++) {
    for (j = i+1; j < n; j++)
      penalty += 2 * fabs(mat_get(mmvn->mvn->sigma, i, j));
    scale += fabs(mat_get(mmvn->mvn->sigma, i, i));
  }

  /* rescale the penalty factor */
  scale *= (n-1);   
  /* If all nondiagonal entries were like the diagonal ones, this
     would be the expected value of the raw penalty.  We will rescale
     it so this expected value is one. */
    
  data->penalty = data->sparsity/scale * penalty; 

  /* calculate gradient in matrix form */

  /* first derive a sign matrix from sigma */
  mat_zero(sign);
  for (i = 0; i < sign->nrows; i++) {
    for (j = i+1; j < sign->ncols; j++) {
      if (mat_get(mmvn->mvn->sigma, i, j) > 0) {
        mat_set(sign, i, j, 1);
        mat_set(sign, j, i, 1);
      }
      else if (mat_get(mmvn->mvn->sigma, i, j) < 0) {
        mat_set(sign, i, j, -1);
        mat_set(sign, j, i, -1);
      }
      /* otherwise leave 0 */
    }
  }
  
  /* multiply by 2R * sparsity for gradient */
  mat_mult(Rgrad, sign, data->R);
  mat_scale(Rgrad, 2.0 * data->sparsity/scale);
    
  /* finally add entries to corresponding gradient components */
  vec_zero(grad);
  for (i = 0; i < data->R->nrows; i++)
    for (j = 0; j < data->R->ncols; j++)
      vec_set(grad, n*d + i*data->R->ncols + j,
              vec_get(grad,  n*d + i*data->R->ncols + j) - 
              mat_get(Rgrad, i, j));
  /* note have to subtract because makes negative contribution */

  mat_free(sign);
  mat_free(Rgrad);
}

/* set scale factor for geometry depending on starting distance matrix */
void nj_set_pointscale(CovarData *data) {
  /* find median pairwise distance */
  double medianD = mat_median_upper_triang(data->dist);  /* off-diagonal median */
  if (medianD <= 0.0 || !isfinite(medianD)) 
    data->pointscale = 1.0;   /* safe backup */
  else if (data->hyperbolic == TRUE)
    data->pointscale = POINTSPAN_HYP / (medianD * sqrt(data->negcurvature));
  else
    data->pointscale = POINTSPAN_EUC / medianD; 
}

/* Like nj_var_sample, but use rejection sampling to sharpen the
   approximate posterior */
List *nj_var_sample_rejection(int nsamples, multi_MVN *mmvn,
                              CovarData *data, TreeModel *mod,
                              FILE *logf) {
  List *retval = lst_new_ptr(nsamples), *init_samples = lst_new_ptr(nsamples * 10);
  Vector *points = vec_new(mmvn->d * mmvn->n);
  double logM = -INFTY, lnl, mvnl, u, avelnl = 0;
  double ll[nsamples * 10], mvnll[nsamples * 10];
  int i, ntot;

  /* collect an initial sample to approximate the upper bound on the
     log ratio of densities */
  for (i = 0; i < nsamples * 10; i++) {
    do {
      nj_sample_points(mmvn, data, points, NULL, NULL);
      nj_points_to_distances(points, data);
      mod->tree = nj_inf(data->dist, data->names, NULL, data);

      if (data->crispr_mod != NULL)
        ll[i] = cpr_compute_log_likelihood(data->crispr_mod, NULL);
      else
        ll[i] = nj_compute_log_likelihood(mod, data, NULL);

      if (!isfinite(ll[i])) tr_free(mod->tree);
    } while (!isfinite(ll[i]));  /* need for the crispr case */

    lst_push_ptr(init_samples, mod->tree);  
    mvnll[i] = mmvn_log_dens(mmvn, points);
    if (ll[i] - mvnll[i] > logM) 
      logM = ll[i] - mvnll[i];
  }

  /* subsample from those by rejection sampling */
  logM += log(1.05); /* provide a little buffer */
  for (i = 0; i < lst_size(init_samples) &&
         lst_size(retval) < nsamples; i++) {
    u = unif_rand();
    if (u < exp(ll[i] - mvnll[i] - logM)) {
      lst_push_ptr(retval, lst_get_ptr(init_samples, i));
      avelnl += ll[i];
    }
    else
      tr_free(lst_get_ptr(init_samples, i));
  }
  ntot = i; /* total number examined */
  for (; i < lst_size(init_samples); i++) /* free any remaining trees */
    tr_free(lst_get_ptr(init_samples, i));
  
  /* if necessary, continue until target is met */
  while (lst_size(retval) < nsamples) {
    ntot++;
    do {
      nj_sample_points(mmvn, data, points, NULL, NULL);
      nj_points_to_distances(points, data);
      mod->tree = nj_inf(data->dist, data->names, NULL, data);

      if (data->crispr_mod != NULL)
        lnl = cpr_compute_log_likelihood(data->crispr_mod, NULL);
      else
        lnl = nj_compute_log_likelihood(mod, data, NULL);

      if (!isfinite(lnl)) tr_free(mod->tree);
    } while (!isfinite(lnl)); /* need for crispr case */
    
    mvnl = mmvn_log_dens(mmvn, points);
    u = unif_rand();
    if (u < exp(lnl - mvnl - logM)) {
      lst_push_ptr(retval, mod->tree);
      avelnl += lnl;
    }
    else
      tr_free(mod->tree);
  } 

  if (logf != NULL)
    fprintf(logf, "# Rejection sampling from %d trees (acceptance rate %.3f); avelnl: %f\n",
            ntot, nsamples*1.0/ntot, avelnl/nsamples);
  
  lst_free(init_samples);
  vec_free(points);
  return(retval);
}

/* alternative versions of gradient calculation. Can be cross-checked
   for debugging */

/* compute the gradient of the log likelihood with respect to the
   individual points by a very simple, fully numerical method.
   Returns log likelihood of model as by product. */
double nj_dL_dx_dumb(Vector *x, Vector *dL_dx, TreeModel *mod, 
                     CovarData *data) {
  double ll, ll_base, xorig, deriv;
  int i, k;
  int n = data->nseqs; /* number of taxa */
  int d = data->dim; /* dimensionality; have to accommodate diagonal case */
  TreeNode *tree, *orig_tree;   /* has to be rebuilt repeatedly; restore at end */

  assert(data->msa != NULL && data->crispr_mod == NULL);
  /* this version not set up for crispr data */
         
  /* set up tree model and get baseline log likelihood */
  nj_points_to_distances(x, data);    
  tree = nj_inf(data->dist, data->msa->names, NULL, data);
  orig_tree = tr_create_copy(tree);   /* restore at the end */
  nj_reset_tree_model(mod, tree);
  ll_base = nj_compute_log_likelihood(mod, data, NULL);

  /* Now perturb each point and propagate perturbation through distance
     calculation, neighbor-joining reconstruction, and likelihood
     calculation on tree */  
  for (i = 0; i < n; i++) {  
    for (k = 0; k < d; k++) {
      int idx = i*d + k;

      xorig = vec_get(x, idx);
      vec_set(x, idx, xorig + DERIV_EPS);

      nj_points_to_distances(x, data); 
      tree = nj_inf(data->dist, data->msa->names, NULL, data);
      nj_reset_tree_model(mod, tree);      
      ll = nj_compute_log_likelihood(mod, data, NULL);
      
      deriv = (ll - ll_base) / DERIV_EPS; 

      vec_set(dL_dx, idx, deriv);
      vec_set(x, idx, xorig); /* restore orig */
    }
  }
  nj_reset_tree_model(mod, orig_tree);
  return ll_base;
}

/* compute the gradient of the log likelihood with respect to the
   individual branch lengths.  This version uses numerical methods
   (mostly useful for testing analytical version) */
double nj_dL_dt_num(Vector *dL_dt, TreeModel *mod, CovarData *data) {
  int nodeidx;
  double ll, ll_base;
  List *traversal;

  if (data->crispr_mod != NULL)
    ll_base = cpr_compute_log_likelihood(data->crispr_mod, NULL);
  else
    ll_base = nj_compute_log_likelihood(mod, data, NULL);

  if (!isfinite(ll_base)) /* can happen with crispr; force calling
                             code to deal with it */
    return ll_base;
  
  /* perturb each branch and recompute likelihood */
  traversal = mod->tree->nodes;
  assert(dL_dt->size == lst_size(traversal) - 1); 
  vec_zero(dL_dt);
  for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
    TreeNode *node = lst_get_ptr(traversal, nodeidx);
    double orig_t = node->dparent;

    if (node == mod->tree || node == mod->tree->rchild)
      continue; /* only consider one branch beneath the root because
                   implicitly unrooted */
    
    node->dparent += DERIV_EPS;

    if (data->crispr_mod != NULL)
      ll = cpr_compute_log_likelihood(data->crispr_mod, NULL);
    else
      ll = nj_compute_log_likelihood(mod, data, NULL);

    if (!isfinite(ll)) /* can happen with crispr; force calling
                          code to deal with it */
      return ll;
    
    vec_set(dL_dt, nodeidx, (ll - ll_base) / DERIV_EPS);
    node->dparent = orig_t;
  }
  return ll_base;  
}

/* compute the Jacobian matrix for 2n-3 branch lengths wrt n-choose-2
   pairwise distances.  This version uses numerical methods and is
   intended for validation of the analytical version */
void nj_dt_dD_num(Matrix *dt_dD, Matrix *D, TreeModel *mod, CovarData *data) {
  TreeNode *tree, *orign, *node;
  int i, j, n = data->msa->nseqs, nodeidx;
  List *trav_tree, *trav_orig;
  
  /* perturb each pairwise distance and measure effect on each branch */
  trav_orig = mod->tree->nodes;
  mat_zero(dt_dD);
  for (i = 0; i < n; i++) {
    for (j = i+1; j < n; j++) {
      double orig_d = mat_get(D, i, j);
      mat_set(D, i, j, orig_d + DERIV_EPS);
      tree = nj_inf(D, data->msa->names, NULL, data);

      /* compare the trees, branch by branch */
      /* we will assume the same topology although that will
         occasionally not be true; good enough for sanity checking */
      trav_tree = tree->nodes;
      for (nodeidx = 0; nodeidx < lst_size(trav_orig); nodeidx++) {
        node = lst_get_ptr(trav_tree, nodeidx);

        if (node == tree || node == tree->rchild) /* unrooted tree */
          continue;
        
        orign = lst_get_ptr(trav_orig, nodeidx);
        
        if (node->id != orign->id) continue;
        
        mat_set(dt_dD, nodeidx, nj_i_j_to_dist(i, j, n),
                (node->dparent - orign->dparent) / DERIV_EPS);
      }
      
      mat_set(D, i, j, orig_d);
      tr_free(tree);
    }
  }
}

/* map the indices of two taxa, i, j to a unique index
   for their pairwise distance.  Unique index will fall densely
   between 0 and n-choose-2 - 1 */
int nj_i_j_to_dist(int i, int j, int n) {
  double ii, jj;
  if (i < j) { ii = i; jj = j; }
  else { ii = j; jj = i; }
  return ((2*n - ii - 1)*ii / 2 + (jj - ii - 1));
}

/* reverse the mapping above: map an index for a pairwise distance to
   the indices for two sequences (s.t. j > i) */
void nj_dist_to_i_j(int pwidx, int *i, int *j, int n) {
  int rowstart;
  *i = 0;
  while (pwidx >= (2*n - (*i) - 1) * (*i) / 2 + (n - (*i) - 1))
    (*i)++;
  rowstart = (2*n - (*i) - 1) * (*i) / 2;
  *j = *i + 1 + (pwidx - rowstart);
}

/* compute the gradient of the log likelihood with respect to the
   individual points by the chain rule and using analytical methods
   for each component.  Fastest but most complicated and error-prone
   version. */
double nj_dL_dx_smartest(Vector *x, Vector *dL_dx, TreeModel *mod,
                         CovarData *data) {
  int n = data->nseqs, nbranches = 2*n-2,  /* have to work with the rooted tree here */
    ndist = n * (n-1) / 2, ndim = data->nseqs * data->dim;
  Vector *dL_dt = vec_new(nbranches);
  Matrix *dt_dD = mat_new(nbranches, ndist);
  Vector *dL_dD = vec_new(ndist);
  TreeNode *tree;
  double ll_base;
  int i, j, d;
  
  /* set up baseline objects */
  nj_points_to_distances(x, data);    
  tree = nj_inf(data->dist, data->names, dt_dD, data);
  nj_reset_tree_model(mod, tree);

  /* TEMPORARY: compare to numerical gradient */
  /* fprintf(stdout, "dt_dD (analytical):\n"); */
  /* mat_print(dt_dD, stdout); */
  /* nj_dt_dD_num(dt_dD, data->dist, mod, data); */
  /* fprintf(stdout, "dt_dD (numerical):\n"); */
  /* mat_print(dt_dD, stdout); */
  /* exit(0); */

  /* calculate log likelihood and analytical gradient */
  if (data->crispr_mod != NULL)
    ll_base = cpr_compute_log_likelihood(data->crispr_mod, dL_dt);
  else
    ll_base = nj_compute_log_likelihood(mod, data, dL_dt);

  if (!isfinite(ll_base)) /* can happen with crispr; force calling
                             code to deal with it */
    return ll_base;
  
  /* TEMPORARY: compare dL_dt to numerical version */
  /* fprintf(stdout, "dL_dt (analytical):\n"); */
  /* vec_print(dL_dt, stdout); */
  /* nj_dL_dt_num(dL_dt, mod, data); */
  /* fprintf(stdout, "dL_dt (numerical):\n"); */
  /* vec_print(dL_dt, stdout); */
  /* exit(0); */
  
  /* apply chain rule to get dL/dD gradient (a vector of dim ndist) */
  mat_vec_mult_transp(dL_dD, dt_dD, dL_dt);
  /* (note taking transpose of both vector and matrix and expressing
     result as column vector) */

  /* finally multiply by dD/dx to obtain final gradient.  This part is
     different for the euclidean and hyperbolic geometries */
  vec_zero(dL_dx);

  if (data->hyperbolic) {
    
    /* first precompute x0[i] = sqrt(1 + ||x_i||^2) */
    double *x0 = (double*)smalloc(n * sizeof(double));
    double alpha = 1.0 / sqrt(data->negcurvature);   /* curvature radius */
    for (i = 0; i < n; i++) {
      double ss = 1.0;
      int base = i * data->dim;
      for (d = 0; d < data->dim; d++) {
        double xid = vec_get(x, base + d);
        ss += xid * xid;
      }
      x0[i] = sqrt(ss);
    }

    /* accumulate pairwise contributions */
    for (i = 0; i < n; i++) {
      double denom_inv, pref;
      int base_i = i * data->dim;
      for (j = i + 1; j < n; j++) {
        int base_j = j * data->dim;

        /* weight = dL/dD_ij */
        double weight = vec_get(dL_dD, nj_i_j_to_dist(i, j, n));

        /* down-weight saturated pairs with large distance */
        double Dij = mat_get(data->dist, i, j);
        if (Dij > 10) weight *= (10 / Dij);           /* soft clip: in (0,1] */
 
        /* u = x0_i*x0_j - <x_i, x_j>  (equals -Lorentz inner product) */
        double dot_spatial = 0.0;
        for (d = 0; d < data->dim; d++) 
          dot_spatial += vec_get(x, base_i + d) * vec_get(x, base_j + d);

        /* prefactor; clamp sqrt(u^2 - 1) for stability */
        double u = x0[i] * x0[j] - dot_spatial;
        
        denom_inv = d_acosh_du_stable(u);         /* = d/du acosh(u) */        
        pref = (alpha / data->pointscale) * denom_inv;
        
        /* dD/dx_i and dD/dx_j contributions */
        for (d = 0; d < data->dim; d++) {
          int idx_i = base_i + d;
          int idx_j = base_j + d;
          double xid = vec_get(x, idx_i);
          double xjd = vec_get(x, idx_j);

          double gi = pref * (-xjd + (x0[j] / x0[i]) * xid);  /* dD_ij/dx_i^d */
          double gj = pref * (-xid + (x0[i] / x0[j]) * xjd);  /* dD_ij/dx_j^d */

          /* accumulate weighted by w_ij = dL/dD_ij */
          vec_set(dL_dx, idx_i, vec_get(dL_dx, idx_i) + weight * gi);
          vec_set(dL_dx, idx_j, vec_get(dL_dx, idx_j) + weight * gj);
        }
      }
    }

    /* add a small radius prior to prevent points from "ballooning" away from zero */
    const double lambda_base = 1e-5; 
    const double lambda_eff  = lambda_base / (data->pointscale*data->pointscale);
      
    for (i = 0; i < ndim; i++)
      vec_set(dL_dx, i, vec_get(dL_dx, i) - 2.0 * lambda_eff * vec_get(x, i));

    sfree(x0);
  }
  else { /* euclidean version is simpler */
    for (i = 0; i < n; i++) {
      int base_i = i * data->dim;
      for (j = i + 1; j < n; j++) {
        int base_j = j * data->dim;
        double dist_ij = mat_get(data->dist, i, j);
        double weight = vec_get(dL_dD, nj_i_j_to_dist(i, j, n));

        if (dist_ij < 1e-15) dist_ij = 1e-15;
        
        for (d = 0; d < data->dim; d++) {
          int idx_i = base_i + d;
          int idx_j = base_j + d;
            
          double coord_diff = vec_get(x, idx_i) - vec_get(x, idx_j);
          double grad_contrib = weight * coord_diff / (dist_ij * data->pointscale * data->pointscale);
          /* need two factors of pointscale, one for the coord_diff, one for the distance */
            
          vec_set(dL_dx, idx_i, vec_get(dL_dx, idx_i) + grad_contrib);
          vec_set(dL_dx, idx_j, vec_get(dL_dx, idx_j) - grad_contrib);
        }
      }
    }

    /* in case of radial flow, need one more step in the chain rule.
       Note that this is supported only in the Euclidean case for now;
       need to move this call outside of 'else' if hyperbolic support
       is added */
    if (data->rf != NULL) {
      /* here what we called dL_dx is actually the gradient wrt the
         radially transformed points, y.  We need to back-propagate
         through the radial flow to obtain the real dL_dx */
      Vector *dL_dy = vec_create_copy(dL_dx);
      rf_backprop(data->rf, x, dL_dx, dL_dy);
      /* note that the gradients wrt the parameters a, b, and ctr are
         computed as side-effects and stored inside rf */
      vec_free(dL_dy);
    }

    /* similar for planar flow */
    if (data->pf != NULL) {
      Vector *dL_dy = vec_create_copy(dL_dx);
      pf_backprop(data->pf, x, dL_dx, dL_dy);
      vec_free(dL_dy);
    }
  }
  
  vec_free(dL_dt);
  mat_free(dt_dD);
  vec_free(dL_dD);

  return ll_base;  
}

/* function used inside NJ algorithm to enable backpropagation of
   derivatives through the algorithm.  Jk is a matrix such that
   element [ij][ab] represents the partial derivative of the distance
   between i and j on iteration k to the distance between a and b at
   the start of the algorithm. The next value, Jnext, is defined
   recursively from Jk.  The last set of identified neighbors are
   denoted f and g, and u denotes the new node created to replace f
   and g. */
void nj_backprop(double *Jk, double *Jnext, int n, int f, int g, int u,
                 Vector *active) {
  int i, a, b;
  int total_nodes = 2*n - 2; /* total possible in final tree */
  int npairs_large = (total_nodes * (total_nodes - 1)) / 2;
  int npairs_small = (n * (n - 1)) / 2;
  
  /* most of Jk will be unchanged so start by copying the whole thing
     efficiently */
  memcpy(Jnext, Jk, npairs_large * npairs_small * sizeof(double));
  
  /* now update distances involving new node u */
  for (i = 0; i < total_nodes; i++) {
    if (vec_get(active, i) == FALSE || i == f || i == g || i == u) continue;  

    int idx_ui = nj_i_j_to_dist(u, i, total_nodes);
    int idx_fi = nj_i_j_to_dist(f, i, total_nodes);
    int idx_gi = nj_i_j_to_dist(g, i, total_nodes);
    int idx_fg = nj_i_j_to_dist(f, g, total_nodes);

    for (a = 0; a < n; a++) {
      for (b = a + 1; b < n; b++) {
        int idx_ab = nj_i_j_to_dist(a, b, n);

        /* recursive update rule for new distance from u to i */
        Jnext[idx_ui*npairs_small+idx_ab] = 0.5 * (Jk[idx_fi*npairs_small+idx_ab] +
                                                   Jk[idx_gi*npairs_small+idx_ab]
                                                   - Jk[idx_fg*npairs_small+idx_ab]);
      }
    }
  }
}

/* helper for nj_backprop_sparse; see below */
static inline void nj_backprop_fast_linear_comb(const SparseVector *rf,   // row idx_fi
                                                const SparseVector *rg,   // row idx_gi
                                                const SparseVector *rfg,  // row idx_fg
                                                SparseVector *rout        // row idx_ui (will be overwritten)
                                                ) {
  /* assumes sorted already */
  spvec_zero(rout); /* clear destination */

  /* use direct access to underlying arrays */
  const SparseVectorElement *af = (SparseVectorElement*)rf->elementlist->array;
  const SparseVectorElement *ag = (SparseVectorElement*)rg->elementlist->array;
  const SparseVectorElement *ah = (SparseVectorElement*)rfg->elementlist->array;
  int nf = lst_size(rf->elementlist),
      ng = lst_size(rg->elementlist),
      nh = lst_size(rfg->elementlist);
  int i = 0, j = 0, k = 0;

  /* merge union of indices with 0.5*(f + g - fg) */
  while (i < nf || j < ng || k < nh) {
    int ci = (i<nf) ? af[i].idx : INT_MAX;
    int cj = (j<ng) ? ag[j].idx : INT_MAX;
    int ck = (k<nh) ? ah[k].idx : INT_MAX;
    int c  = (ci<cj ? (ci<ck ? ci:ck) : (cj<ck?cj:ck));

    double vf = (ci==c) ? af[i].val : 0.0;
    double vg = (cj==c) ? ag[j].val : 0.0;
    double vh = (ck==c) ? ah[k].val : 0.0;
    double v  = 0.5 * (vf + vg - vh);

    if (v != 0.0) spvec_set_sorted(rout, c, v);

    if (ci==c) ++i;
    if (cj==c) ++j;
    if (ck==c) ++k;
  }
}

/* version of function above that uses a sparse matrix implementation
   to avoid explosion in size of Jk and Jnext */
void nj_backprop_sparse(SparseMatrix *Jk, SparseMatrix *Jnext, int n, int f, int g, int u,
                        Vector *active) {
  int i;
  int total_nodes = 2*n - 2; /* total possible in final tree */

  /* double buffering Jk and Jnext to avoid expensive copy.  Keep
     track of which destination rows are rebuilt */
  unsigned char *touched = (unsigned char*)calloc(Jk->nrows, 1);
  if (!touched) die("nj_backprop_sparse: out of memory\n");
  
  /* now update distances involving new node u */
  for (i = 0; i < total_nodes; i++) {
    if (vec_get(active, i) == FALSE || i == f || i == g || i == u) continue;  

    int idx_ui = nj_i_j_to_dist(u, i, total_nodes);
    int idx_fi = nj_i_j_to_dist(f, i, total_nodes);
    int idx_gi = nj_i_j_to_dist(g, i, total_nodes);
    int idx_fg = nj_i_j_to_dist(f, g, total_nodes);

    /* ensure destination row is deep copy before writing */
    spmat_replace_row_empty(Jnext, idx_ui, lst_size(Jk->rows[idx_fi]->elementlist) +
                            lst_size(Jk->rows[idx_gi]->elementlist) +
                            lst_size(Jk->rows[idx_fg]->elementlist));
    
    /* use helper function to combine values from three rows in one pass */
    nj_backprop_fast_linear_comb(Jk->rows[idx_fi], Jk->rows[idx_gi],
                                 Jk->rows[idx_fg], Jnext->rows[idx_ui]);

    touched[idx_ui] = 1;
  }

  /* Alias all untouched rows from Jk into Jnext (avoid data copy) */
  for (int r = 0; r < Jk->nrows; r++) {
    if (touched[r] == 1) continue;
    if (Jnext->rows[r] == Jk->rows[r]) continue;  /* already aliased (rare) */
    spvec_release(Jnext->rows[r]);
    Jnext->rows[r] = Jk->rows[r];
    spvec_retain(Jnext->rows[r]);
  }

  free(touched);
}

/* initialize Jk at the beginning of the NJ alg */
void nj_backprop_init(double *Jk, int n) {
  int i, j, total_nodes = 2*n - 2;
  int npairs_large = (total_nodes * (total_nodes - 1)) / 2;
  int npairs_small = (n * (n - 1)) / 2;

  memset(Jk, 0, npairs_large * npairs_small * sizeof(double));
  
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      int idx_ij_row = nj_i_j_to_dist(i, j, total_nodes);
      int idx_ij_col = nj_i_j_to_dist(i, j, n);
      Jk[idx_ij_row*npairs_small+idx_ij_col] = 1;  /* deriv of d_{ij} wrt itself */
    }
  }
}

/* version that uses sparse matrix */
void nj_backprop_init_sparse(SparseMatrix *Jk, int n) {
  int i, j, total_nodes = 2*n - 2;

  spmat_zero(Jk); 
  
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      int idx_ij_row = nj_i_j_to_dist(i, j, total_nodes);
      int idx_ij_col = nj_i_j_to_dist(i, j, n);
      spmat_set_sorted(Jk, idx_ij_row, idx_ij_col, 1.0);  /* deriv of d_{ij} wrt itself */
    }
  }
}

/* for use in backpropation with NJ.  Sets appropriate elements of
   Jacobian dt_dD after two neighbors f and g are joined */
void nj_backprop_set_dt_dD(double *Jk, Matrix *dt_dD, int n, int f, int g,
                           int branch_idx_f, int branch_idx_g, Vector *active) {
  int a, b, m;
  int total_nodes = 2*n - 2;
  int idx_fg = nj_i_j_to_dist(f, g, total_nodes);
  int npairs_small = (n * (n - 1)) / 2;
  int nk = vec_sum(active);

  /* the final call, with nk = 2, is a special case */
  if (nk == 2) {
    /* directly set the final branch derivative */
    for (a = 0; a < n; a++) {
      for (b = a + 1; b < n; b++) {
        int idx_ab = nj_i_j_to_dist(a, b, n);
        
        /* branch derivative is equal to the value in Jk times 1/2
           because of the way we split the last branch in the unrooted
           tree */
        mat_set(dt_dD, branch_idx_f, idx_ab, 0.5 * Jk[idx_fg*npairs_small+idx_ab]);
      }
    }
    return;
  }
    
  /* branch derivative for f -> u */
  for (a = 0; a < n; a++) {
    for (b = a + 1; b < n; b++) {
      int idx_ab = nj_i_j_to_dist(a, b, n);
      double sum_diff = 0;

      for (m = 0; m < total_nodes; m++) {
        if (vec_get(active, m) == FALSE || m == f || m == g)
          continue;

        int idx_fm = nj_i_j_to_dist(f, m, total_nodes);
        int idx_gm = nj_i_j_to_dist(g, m, total_nodes);
        sum_diff += Jk[idx_fm*npairs_small+idx_ab] - Jk[idx_gm*npairs_small+idx_ab];
      }

      mat_set(dt_dD, branch_idx_f, idx_ab, 0.5 * Jk[idx_fg*npairs_small+idx_ab] +
              (0.5 / (nk - 2)) * sum_diff);
      assert(isfinite(mat_get(dt_dD, branch_idx_f, idx_ab)));
    }
  }

  /* branch derivative for g -> u */
  for (a = 0; a < n; a++) {
    for (b = a + 1; b < n; b++) {
      int idx_ab = nj_i_j_to_dist(a, b, n);
      mat_set(dt_dD, branch_idx_g, idx_ab, Jk[idx_fg*npairs_small+idx_ab] -
        mat_get(dt_dD, branch_idx_f, idx_ab));
    }
  }
}

/* version that uses sparse matrix */
void nj_backprop_set_dt_dD_sparse(SparseMatrix *Jk, Matrix *dt_dD, int n, int f, int g,
                                  int branch_idx_f, int branch_idx_g, Vector *active) {
  int total_nodes = 2*n - 2;
  int idx_fg = nj_i_j_to_dist(f, g, total_nodes);
  int nk = vec_sum(active);
  int n_ab = n*(n-1)/2;
  double *sum_diff = malloc(sizeof(double) * n_ab);

  /* the final call, with nk = 2, is a special case */
  if (nk == 2) {
    /* just copy 0.5 * row idx_fg into branch f */
    const SparseVector *rfg = Jk->rows[idx_fg];

    /* first zero the whole dt_dD row */
    for (int ab = 0; ab < n_ab; ab++) mat_set(dt_dD, branch_idx_f, ab, 0.0);
    const SparseVectorElement *a = (SparseVectorElement*)rfg->elementlist->array;
    int nz = lst_size(rfg->elementlist);
    /* now fill in non-zero entries */
    for (int t = 0; t < nz; t++)
      mat_set(dt_dD, branch_idx_f, a[t].idx, 0.5 * a[t].val);

    free(sum_diff);
    return;
  }

  for (int ab = 0; ab < n_ab; ab++) sum_diff[ab] = 0.0;
  
  /* for each active m (excluding f,g), accumulate sparse diffs */
  for (int m = 0; m < total_nodes; m++) {
    if (vec_get(active, m) == FALSE || m == f || m == g)
      continue;

    int idx_fm = nj_i_j_to_dist(f, m, total_nodes);
    int idx_gm = nj_i_j_to_dist(g, m, total_nodes);

    const SparseVectorElement *rf = (SparseVectorElement*)Jk->rows[idx_fm]->elementlist->array;
    const SparseVectorElement *rg = (SparseVectorElement*)Jk->rows[idx_gm]->elementlist->array;
    int nf = lst_size(Jk->rows[idx_fm]->elementlist);
    int ng = lst_size(Jk->rows[idx_gm]->elementlist);
    int i = 0, j = 0;

    /* branch derivative for f->u; first merge the two rows then
       scatter (+1 for f, -1 for g) */
    while (i < nf || j < ng) {
      int cf = (i<nf) ? rf[i].idx : INT_MAX;
      int cg = (j<ng) ? rg[j].idx : INT_MAX;
      if (cf == cg) {
        sum_diff[cf] += (rf[i].val - rg[j].val);
        i++; j++;
      }
      else if (cf < cg) {
        sum_diff[cf] += rf[i].val;
        i++;
      }
      else {
        sum_diff[cg] -= rg[j].val;
        j++;
      }
    }
  }

  /* set dt_dD row for branch f: 0.5*Jk[idx_fg,:] + (0.5/(nk-2))*sum_diff */
  const SparseVectorElement *ah = (SparseVectorElement*)Jk->rows[idx_fg]->elementlist->array;
  int nh = lst_size(Jk->rows[idx_fg]->elementlist);

  /* start from (0.5/(nk-2))*sum_diff (dense) */
  double scale = 0.5 / (nk - 2);
  for (int ab = 0; ab < n_ab; ab++)
    mat_set(dt_dD, branch_idx_f, ab, scale * sum_diff[ab]);

  /* now add the sparse 0.5 * Jk[idx_fg,:] */
  for (int t = 0; t < nh; t++) {
    int ab = ah[t].idx;
    mat_set(dt_dD, branch_idx_f, ab,
            mat_get(dt_dD, branch_idx_f, ab) + 0.5 * ah[t].val);
  }

  /* branch derivative for g->u; g = Jk[idx_fg,:] - branch f 
     do it densely; add sparse afterwards */
  for (int ab = 0; ab < n_ab; ab++) {
    double jf = mat_get(dt_dD, branch_idx_f, ab);
    mat_set(dt_dD, branch_idx_g, ab, -jf); /* temporarily 0, will add Jk[idx_fg,:] */ 
  }

  for (int t = 0; t < nh; t++) {
    int ab = ah[t].idx;
    mat_set(dt_dD, branch_idx_g, ab,
            mat_get(dt_dD, branch_idx_g, ab) + ah[t].val);
  }
  
  free(sum_diff);
}

/* wrapper for various distance-based tree inference algorithms */
TreeNode *nj_inf(Matrix *D, char **names, Matrix *dt_dD,
                 CovarData *covar_data) {
  if (covar_data->ultrametric) {
    TreeNode *t = upgma_fast_infer(D, names, dt_dD);

    if (covar_data->no_zero_br == TRUE)
      nj_repair_zero_br(t);
    return t;
  }
  else
    return nj_fast_infer(D, names, dt_dD);
}

/* helper functions for nuisance parameters in variational
   inference. For now these include only the HKY ti/tv parameter for
   DNA models and the silencing rate and leading branch length for
   CRISPR models */
int nj_get_num_nuisance_params(TreeModel *mod, CovarData *data) {
  int retval = 0;

  if (data->crispr_mod != NULL)
    retval += 2;
  else if (mod->subst_mod == HKY85)
    retval += 1;

  if (data->rf != NULL)
    retval += data->rf->ctr->size + 2;

  if (data->pf != NULL)
    retval += data->pf->ndim * 2 + 1;
  
  return retval;
}

char *nj_get_nuisance_param_name(TreeModel *mod, CovarData *data, int idx) {
  assert(idx >= 0);
  if (data->crispr_mod != NULL) {
    if (idx == 0)
      return "nu";
    else if (idx == 1)
      return ("lead_t");
    else idx -= 2;  /* incrementally subtract each set of indices */
  }
  else if (mod->subst_mod == HKY85) {
    if (idx == 0) 
      return "kappa";
    else idx -= 1;
  }
  
  if (data->rf != NULL) {
    if (idx < data->rf->ctr->size) {
      char *tmp = smalloc(15 * sizeof(char));
      snprintf(tmp, 15, "rf_ctr[%d]", idx);
      return tmp;
    }
    else {
      idx -= data->rf->ctr->size;
      if (idx == 0)
        return "rf_a";
      else if (idx == 1)
        return "rf_b";
      else idx -= 2;
    }
  }

  if (data->pf != NULL) {
    if (idx < data->pf->ndim) {
      char *tmp = smalloc(15 * sizeof(char));
      snprintf(tmp, 15, "pf_u[%d]", idx);
      return tmp;
    }
    else if (idx < 2 * data->pf->ndim) {
      idx -= data->pf->ndim;
      char *tmp = smalloc(15 * sizeof(char));
      snprintf(tmp, 15, "pf_w[%d]", idx);
      return tmp;
    }
    else {
      idx -= 2 * data->pf->ndim;
      if (idx == 0)
        return "pf_b";
      else
        die("ERROR in nj_get_nuisance_param_name: index out of bounds.\n");
    }
  }
    
  return NULL;
}

/* update nuis_grad based on current gradients */
void nj_update_nuis_grad(TreeModel *mod, CovarData *data, Vector *nuis_grad) {
  int idx = 0;
  if (data->crispr_mod != NULL) {
    vec_set(nuis_grad, idx++, data->crispr_mod->deriv_sil);
    vec_set(nuis_grad, idx++, data->crispr_mod->deriv_leading_t);
  }
  else if (mod->subst_mod == HKY85) {
    vec_set(nuis_grad, idx++, data->deriv_hky_kappa);
  }

  if (data->rf != NULL) {
    for (int i = 0; i < data->rf->ctr_grad->size; i++)
      vec_set(nuis_grad, idx++, vec_get(data->rf->ctr_grad, i));
    vec_set(nuis_grad, idx++, data->rf->a_grad);
    vec_set(nuis_grad, idx++, data->rf->b_grad);
  }

  if (data->pf != NULL) {
    for (int i = 0; i < data->pf->ndim; i++)
      vec_set(nuis_grad, idx++, vec_get(data->pf->u_grad, i));
    for (int i = 0; i < data->pf->ndim; i++)
      vec_set(nuis_grad, idx++, vec_get(data->pf->w_grad, i));
    vec_set(nuis_grad, idx++, data->pf->b_grad);
  }
  
  assert(idx == nuis_grad->size);
}

/* save current values of nuisance params */
void nj_save_nuis_params(Vector *stored_vals, TreeModel *mod, CovarData *data) {
  int idx = 0;
  if (data->crispr_mod != NULL) {
    vec_set(stored_vals, idx++, data->crispr_mod->sil_rate);
    vec_set(stored_vals, idx++, data->crispr_mod->leading_t);
  }
  else if (mod->subst_mod == HKY85) 
    vec_set(stored_vals, idx++, data->hky_kappa);

  if (data->rf != NULL) {
    for (int i = 0; i < data->rf->ctr->size; i++)
      vec_set(stored_vals, idx++, vec_get(data->rf->ctr, i));
    vec_set(stored_vals, idx++, data->rf->a);
    vec_set(stored_vals, idx++, data->rf->b);
  }

  if (data->pf != NULL) {
    for (int i = 0; i < data->pf->ndim; i++)
      vec_set(stored_vals, idx++, vec_get(data->pf->u, i));
    for (int i = 0; i < data->pf->ndim; i++)
      vec_set(stored_vals, idx++, vec_get(data->pf->w, i));
    vec_set(stored_vals, idx++, data->pf->b);
  }
    
  assert(idx == stored_vals->size);
}

/* update all nuisance parameters based on vector of stored values */
void nj_update_nuis_params(Vector *stored_vals, TreeModel *mod, CovarData *data) {
  int idx = 0;
  if (data->crispr_mod != NULL) {
    data->crispr_mod->sil_rate = vec_get(stored_vals, idx++);
    data->crispr_mod->leading_t = vec_get(stored_vals, idx++);
  }
  else if (mod->subst_mod == HKY85) {
    data->hky_kappa = vec_get(stored_vals, idx++);
    tm_set_HKY_matrix(mod, data->hky_kappa, -1);
    tm_scale_rate_matrix(mod);
    mm_diagonalize(mod->rate_matrix);
  }

  if (data->rf != NULL) {
    /* TEMPORARY: disable update of center */
    for (int i = 0; i < data->rf->ctr->size; i++)
      idx++;
      /* vec_set(data->rf->ctr, i, vec_get(stored_vals, idx++)); */
    data->rf->a = vec_get(stored_vals, idx++);
    data->rf->b = vec_get(stored_vals, idx++);
    rf_update(data->rf);
  }

  if (data->pf != NULL) {
    for (int i = 0; i < data->pf->ndim; i++)
      vec_set(data->pf->u, i, vec_get(stored_vals, idx++)); 
    for (int i = 0; i < data->pf->ndim; i++)
      vec_set(data->pf->w, i, vec_get(stored_vals, idx++)); 
    data->pf->b = vec_get(stored_vals, idx++);
  }
    
  assert(idx == stored_vals->size);
}

/* add to single nuisance parameter */
void nj_nuis_param_pluseq(TreeModel *mod, CovarData *data, int idx, double inc) {
  if (data->crispr_mod != NULL) {
    if (idx == 0) {
      data->crispr_mod->sil_rate += inc;
      if (data->crispr_mod->sil_rate < 0)
        data->crispr_mod->sil_rate = 0;
    }
    else if (idx == 1) {
      data->crispr_mod->leading_t += inc;
      if (data->crispr_mod->leading_t < 0)
        data->crispr_mod->leading_t = 0;
    }
    else
      idx -= 2; /* subtract for below */
  }
  else if (mod->subst_mod == HKY85) {
    if (idx == 0) {
      data->hky_kappa += inc;
      if (data->hky_kappa < 0)
        data->hky_kappa = 0;
      tm_set_HKY_matrix(mod, data->hky_kappa, -1);
      tm_scale_rate_matrix(mod);
      mm_diagonalize(mod->rate_matrix);
    }
    else idx -= 1;      
  }

  if (data->rf != NULL) {
    if (idx < data->rf->ctr->size)
      ;  /* TEMPORARY */
      /* vec_set(data->rf->ctr, idx, vec_get(data->rf->ctr, idx) + inc); */
    else {
      idx -= data->rf->ctr->size;
      if (idx == 0)
        data->rf->a += inc;
      else if (idx == 1)
        data->rf->b += inc;
      else
        idx -= 2;
    }
  }

  if (data->pf != NULL) {
    if (idx < data->pf->ndim)
      vec_set(data->pf->u, idx, vec_get(data->pf->u, idx) + inc);
    else if (idx < 2 * data->pf->ndim) {
      idx -= data->pf->ndim;
      vec_set(data->pf->w, idx, vec_get(data->pf->w, idx) + inc);
    }
    else {
      idx -= 2 * data->pf->ndim;
      if (idx == 0)
        data->pf->b += inc;
      else
        die("ERROR in nj_nuis_param_pluseq: index out of bounds.\n");
    }
  }
}

/* return value of single nuisance parameter */
double nj_nuis_param_get(TreeModel *mod, CovarData *data, int idx) {
  if (data->crispr_mod != NULL) {
    if (idx == 0)
      return data->crispr_mod->sil_rate;
    else if (idx == 1)
      return data->crispr_mod->leading_t;
    else
      idx -= 2; /* subtract for below */
  }
  else if (mod->subst_mod == HKY85) {
    if (idx == 0)
      return data->hky_kappa;
    else idx -= 1;      
  }

  if (data->rf != NULL) {
    if (idx < data->rf->ctr->size) 
      return vec_get(data->rf->ctr, idx);
    else {
      idx -= data->rf->ctr->size;
      if (idx == 0)
        return data->rf->a;
      else if (idx == 1)
        return data->rf->b;
      else
        idx -= 2;
    }
  }

  if (data->pf != NULL) {
    if (idx < data->pf->ndim) 
      return vec_get(data->pf->u, idx);
    else if (idx < 2 * data->pf->ndim) {
      idx -= data->pf->ndim;
      return vec_get(data->pf->w, idx);
    }
    else {
      idx -= 2 * data->pf->ndim;
      if (idx == 0)
        return data->pf->b;
      else
        die("ERROR in nuis_param_get: index out of bounds.\n");
    }
  }
  
  return -1;
}

void nj_repair_zero_br(TreeNode *t) {
  for (int nodeidx = 0; nodeidx < lst_size(t->nodes); nodeidx++) {
    TreeNode *n = lst_get_ptr(t->nodes, nodeidx);
    if (n->parent != NULL && n->dparent <= 0)
      n->dparent = 1e-3;
  }
}

/* sample points from variational distribution.  This is a wrapper
   that encapsulates the radial flow layer (if selected) and
   antithetic sampling.  If points_std is non-NULL, it will be used to
   store the baseline standard normal variate for use in downstream
   calculations in variational inference. Antithetic sampling is only
   used in this case */
void nj_sample_points(multi_MVN *mmvn, CovarData *data, Vector *points, Vector *points_std,
                      double *logdet) {
  static int i = 0;
  static Vector *cachedpoints = NULL, *cachedstd = NULL;
  Vector *tmppoints = (data->rf != NULL || data->pf != NULL) ?
    vec_new(points->size) : NULL;
      
  if (points_std == NULL) 
    mmvn_sample(mmvn, points); /* simple in this case */
  else {
    /* otherwise we have to make use of caching for antithetic sampling */
    if (cachedpoints != NULL && cachedpoints->size != points->size) {
      vec_free(cachedpoints);
      vec_free(cachedstd);
      cachedpoints = NULL; /* force realloc */
    }
    if (cachedpoints == NULL) {
      cachedpoints = vec_new(points->size);
      cachedstd = vec_new(points->size);
      i = 0; /* force new sample */
    }
    
    if (i % 2 == 0) { /* new sample, update caches */
      mmvn_sample_anti_keep(mmvn, points, cachedpoints, points_std);
      vec_copy(cachedstd, points_std);
    }
    else { /* just use cache to define sample */
      vec_copy(points, cachedpoints);
      vec_copy(points_std, cachedstd);
      vec_scale(points_std, -1.0);
    }
    i++;
  }
  
  /* finally apply the normalizing flows if activated */
  if (data->rf != NULL) {
    double ldet = rf_forward(data->rf, tmppoints, points);
    vec_copy(points, tmppoints);
    if (logdet != NULL) (*logdet) = ldet; /* for use in calling code */
  }

  if (data->pf != NULL) {
    double ldet = pf_forward(data->pf, tmppoints, points);
    vec_copy(points, tmppoints);
    if (logdet != NULL) (*logdet) += ldet; /* for use in calling code */
  }

  if (tmppoints != NULL) vec_free(tmppoints);    
}

