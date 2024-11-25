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
#include "phast/stacks.h"
#include "phast/trees.h"
#include "phast/misc.h"
#include "phast/stringsplus.h"
#include "phast/nj.h"
/* #include <gsl.h>*/
#include <gsl/gsl_randist.h>

/* Reset Q matrix based on distance matrix.  Assume upper triangular
   square Q and D.  Only touches active rows and columns of Q and D up
   to maxidx. As a side-effect set u and v to the indices of the
   closest neighbors.  Also update sums to sum of distances from each
   node */
/* FIXME: pass in max n to consider, avoid wasted operations */
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

  mat_set(D, u, w, 0.5 * mat_get(D, u, v) +
	  1.0/(2.0*(n-2)) * (vec_get(sums, u) - vec_get(sums, v)));

  mat_set(D, v, w, mat_get(D, u, v) - mat_get(D, u, w));
    
	     
  vec_set(active, u, FALSE);
  vec_set(active, v, FALSE);

  for (k = 0; k < w; k++) {
    if (vec_get(active, k) == TRUE) {
      double du, dv;

      /* needed because of upper triangular constraint */
      du = (u < k ? mat_get(D, u, k) : mat_get(D, k, u));
      dv = (v < k ? mat_get(D, v, k) : mat_get(D, k, v));
      
      mat_set(D, k, w, 0.5 * (du + dv - mat_get(D, u, v)));
    }
  }
  
  vec_set(active, w, TRUE);
}


/* Main function to infer the tree from a starting distance matrix */
TreeNode* nj_infer_tree(Matrix *initD, char **names) {
    int n = initD->nrows;
    int N = 2*n - 2;   /* number of nodes in unrooted tree */
    int i, j, u = -1, v = -1, w;
    Matrix *D, *Q;
    Vector *sums, *active;
    List *nodes;  /* could just be an array */
    TreeNode *node_u, *node_v, *node_w, *root;
    
    if (initD->nrows != initD->ncols || n < 3)
      die("ERROR nj_infer_tree: bad distance matrix\n");

    /* create a larger distance matrix of dimension N x N to
       accommodate internal nodes; also set up list of active nodes
       and starting tree nodes */
    D = mat_new(N, N); mat_zero(D);
    active = vec_new(N); vec_set_all(active, FALSE);
    sums = vec_new(N); vec_zero(sums);
    nodes = lst_new_ptr(N);
    
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
    
    /* finish set up of tree */
    root->nnodes = N+1;
    
    lst_free(nodes);
    vec_free(active);
    vec_free(sums);
    mat_free(D);
    mat_free(Q);

    
    return root;

    /* FIXME: perhaps use a mapping from source indices to Q indices
       so that Q can stay small?  keep track of size of Q on each
       iteration */
}


/* compute pairwise distance between two DNA seqs using the
   Jukes-Cantor model */
double nj_compute_JC_dist(MSA *msa, int i, int j) {
  int k, diff = 0, n = 0;
  for (k = 0; k < msa->length; k++) {
    if (msa->seqs[i][k] == GAP_CHAR || msa->seqs[j][k] == GAP_CHAR ||
	msa->is_missing[(int)msa->seqs[i][k]] ||
	msa->is_missing[(int)msa->seqs[j][k]])
      continue;
    n++;
    if (msa->seqs[i][k] != msa->seqs[j][k])
      diff++;
  }
  return -0.75 * log(1 - 4.0/3 * diff/n);   /* Jukes-Cantor correction */
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


/* sample a vector from a multivariate normal distribution with the
   given mean vector and covariance matrix. Uses the GSL under the hood */
void nj_sample_mvn(Vector *mu, Matrix *sigma, Vector *retval) {
  gsl_vector *gsl_mu = gsl_vector_alloc(mu->size);
  gsl_matrix *gsl_sigma = gsl_matrix_alloc(sigma->nrow, sigma->ncol);
  static gsl_rng *r = NULL;
  int i, j;

  if (r == NULL) 
    r = gsl_rng_alloc(gsl_rng_taus);
  
  /* convert to gsl objects */
  for (i = 0; i < mu->size; i++)
    gsl_vec_set(gsl_mu, i, vec_get(mu, i));
  for (i = 0; i < sigma->nrow; i++)
    for (j = 0; j < sigma->ncol; j++)
      gsl_mat_set(gsl_sigma, i, j, mat_get(sigma, i, j));
 
  gsl_linalg_cholesky_decomp1(gsl_sigma);
  gsl_ran_multivariate_gaussian(r, gsl_mu, gsl_sigma, gsl_retval);

  /* copy result back to phast vector */
  for (i = 0; i < mu->size; i++)
    vec_set(retval, i, gsl_vector_get(gsl_retval, i));

  gsl_vector_free(gsl_mu);
  gsl_matrix_free(gsl_sigma);
}

/* convert an nd-dimensional vector to an nxn upper triangular
   distance matrix.  Assumes each taxon is represented as a point in
   d-dimensional space and computes Euclidean distances between these
   points */ 
void nj_points_to_distances(Vector *points, Vector *D) {
  int i, j, k, vidx1, vidx2, n, d;
  double sum;

  n = D->nrow;
  d = points->size / n;
  
  if (points->size != n*d || D->nrow != D->ncol) {
    die("ERROR nj_points_to_distances: bad dimensions\n");
  }

  mat_zero(D);
  for (i = 0; i < n; i++) {
    vidx1 = i*d;
    for (j = i+1; j < n; j++) {
      vidx2 = j*d;
      sum = 0;
      for (k = 0; k < d; k++) {
	sum += pow(vec_get(points, vidx1 + k) - vec_get(points, vidx2 + k), 2);
      }
      mat_set(D, vidx1, vidx2, sqrt(sum));
    }
  }
}

/* sample a tree from a multivariate normal averaging distribution
   with the given mean vector and covariance matrix.  Each taxon is
   represented as a point in a d-dimensional, distances are assumed to
   be Euclidean, and tree is computed by neighbor-joining */
TreeNode* nj_mvn_sample_tree(Vector *mu, Matrix *sigma, int n, char **names) {
  int d = mu->size / n;
  TreeNode *tree;
  Matrix *D;
  Vector *points;

  if (mu->size != n*d || sigma->nrow != mu->size || sigma->ncol != mu->size)
    die("ERROR nj_mvn-sample_tree: bad dimensions\n");

  points = vec_new(mu->size);
  D = mat_new(n, n);
  
  nj_sample_mvn(mu, sigma, points);
  nj_points_to_distances(points, D);
  tree = nj_infer_tree(D, names);
  
  vec_free(points);
  mat_free(D);

  return(tree); 
}

/* compute the gradient of the log likelihood for a tree model with
   respect to the free parameters for the MVN averaging distribution.
   This version uses numerical methods */

/*
take a tree model as input (?)    [try to reuse if possible]
also take alignment
  
loop through params
perturb each param by small epsilon
call sample tree
replace tree in model
recalculate likelihood
calculate derivative and put in return vector

  be careful about other parts of tree model; maybe instantiate new one each time?
*/


/* optimize variational model by gradient ascent.  Takes initial tree model and alignment and distance matrix, dimensionality of Euclidean space to work in */

/*
  create starting values of mu from distances.  can do from gram matrix but requires more code.  for now just initialize randomly.  better test anyway


  iterate until (approximate) convergence
  inner loop: iterate over minibatch of MVN samples
  for each sample:
  sample tree, build tree model, compute gradient wrt free parameters [above]
  add sparsity penalty based on current draw
  average all of those values for the minibatch and update parameters at given learning rate
  report parameter values to a log file
*/
