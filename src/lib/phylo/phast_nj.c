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
#include "phast/tree_likelihoods.h"
#include "phast/eigen.h"
#include "phast/sufficient_stats.h"
#include "phast/markov_matrix.h"
#include "phast/mvn.h"

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
  if (mat_get(D, u, w) < 0)
    mat_set(D, u, w, 0);
  if (mat_get(D, v, w) < 0)
    mat_set(D, v, w, 0);
  
  vec_set(active, u, FALSE);
  vec_set(active, v, FALSE);

  for (k = 0; k < w; k++) {
    if (vec_get(active, k) == TRUE) {
      double du, dv;

      /* needed because of upper triangular constraint */
      du = (u < k ? mat_get(D, u, k) : mat_get(D, k, u));
      dv = (v < k ? mat_get(D, v, k) : mat_get(D, k, v));
      
      mat_set(D, k, w, 0.5 * (du + dv - mat_get(D, u, v)));

      if (mat_get(D, k, w) < 0)
	mat_set(D, k, w, 0);
    }
  }
  
  vec_set(active, w, TRUE);
}


/* Main function to infer the tree from a starting distance matrix.
   Does not alter the provided distance matrix */
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
    tr_set_nnodes(root);

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
  
  assert(tree->nodes != NULL);  /* assume list of nodes exists */
  
  for (i = 0; i < tree->nnodes; i++) {
    n1 = lst_get_ptr(tree->nodes, i);
    if (n1->lchild == NULL && n1->rchild == NULL)
      lst_push_ptr(leaves, n1);
  }

  if (lst_size(leaves) != n)
    die("ERROR in nj_tree_to_distances: number of names must match number of leaves in tree.\n");
  
  D = mat_new(n, n);
  mat_zero(D);

  seq_idx = nj_build_seq_idx(leaves, names);
  
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
  for (id = 0; id <= root->nnodes; id++)
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
   d-dimensional space and computes Euclidean distances between these
   points */ 
void nj_points_to_distances(Vector *points, Matrix *D) {
  int i, j, k, vidx1, vidx2, n, d;
  double sum;

  n = D->nrows;
  d = points->size / n;
  
  if (points->size != n*d || D->nrows != D->ncols) {
    die("ERROR nj_points_to_distances: bad dimensions\n");
  }

  mat_zero(D);
  for (i = 0; i < n; i++) {
    vidx1 = i*d;
    for (j = i+1; j < n; j++) {
      vidx2 = j*d;
      sum = 0;
      for (k = 0; k < d; k++) {
        sum += pow(vec_get(points, vidx1 + k) -
                   vec_get(points, vidx2 + k), 2);
      }
      mat_set(D, i, j, sqrt(sum));
    }
  }
}

/* convert an nd-dimensional vector to an nxn upper triangular
   distance matrix.  Assumes each taxon is represented as a point in
   d-dimensional space and computes hyperbolic distances between these
   points */ 
void nj_points_to_distances_hyperbolic(Vector *points, Matrix *D, double negcurvature) {
  int i, j, k, vidx1, vidx2, n, d;
  double lor_inner, ss1, ss2, x0_1, x0_2;

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
        lor_inner += vec_get(points, vidx1 + k) * vec_get(points, vidx2 + k);
        ss1 += pow(vec_get(points, vidx1 + k), 2);
        ss2 += pow(vec_get(points, vidx2 + k), 2);
      }
      x0_1 = sqrt(ss1); /* the 0th dimension for each point is determined by the
                           others, to stay on the hyperboloid */
      x0_2 = sqrt(ss2);

      lor_inner -= x0_1 * x0_2;  /* last term of Lorentz inner product */
      if (fabs(1.0 + lor_inner) < 1.0e-8)
        mat_set(D, i, j, 0);
      else
        mat_set(D, i, j, 1/sqrt(negcurvature) * acosh(-lor_inner));
      /* distance between two points on the sheet, scaled by the curvature */
    }
  }
}

/* compute the gradient of the log likelihood for a tree model with
   respect to the free parameters of the MVN averaging distribution,
   starting from a given MVN sample (points). Returns log likelihood
   of current model, which is computed as a by-product.  This version
   uses numerical methods */
double nj_compute_model_grad(TreeModel *mod, MVN *mvn, MSA *msa,
                             unsigned int hyperbolic, double negcurvature,
                             Vector *points, Vector *grad, Matrix *D,
                             Vector *sigmapar, enum covar_type covar_param,
                             CovarData *data) {
  int n = msa->nseqs;
  int d = mvn->dim / n;
  int i, k;
  double porig, ll_base, ll, deriv, stdrv, sd, lambda_grad, sqrtl;
  TreeNode *tree, *orig_tree;   /* has to be rebuilt repeatedly; restore at end */
  
  if (mvn->dim != n*d)
    die("ERROR in nj_compute_model_grad: bad parameters\n");

  if (covar_param == DIST) {
    if (grad->size != mvn->dim + 1)
      die("ERROR in nj_compute_model_grad: bad gradient dimension.\n");
    if (data->cholL == NULL)
      die("ERROR in nj_compute_model_grad: Cholesky decomposition required in DIST case.\n");
  }
  else if (grad->size != 2*mvn->dim)  /* DIAG case */
    die("ERROR in nj_compute_model_grad: bad gradient dimension.\n");
  
  /* set up tree model and get baseline log likelihood */
  if (hyperbolic)
    nj_points_to_distances_hyperbolic(points, D, negcurvature);
  else
    nj_points_to_distances(points, D);
    
  tree = nj_infer_tree(D, msa->names);
  orig_tree = tr_create_copy(tree);   /* restore at the end */
  nj_reset_tree_model(mod, tree);
  ll_base = nj_compute_log_likelihood(mod, msa, NULL);
  
  /* Perturb each point and propagate perturbation through distance
     calculation, neighbor-joining reconstruction, and likelihood
     calculation on tree */

  vec_zero(grad);
  lambda_grad = 0;
  if (covar_param == DIST)
    sqrtl = sqrt(data->lambda);

  for (i = 0; i < n; i++) {  
    for (k = 0; k < d; k++) {
      int pidx = i*d + k;

      porig = vec_get(points, pidx);
      vec_set(points, pidx, porig + DERIV_EPS);

      if (hyperbolic)
        nj_points_to_distances_hyperbolic(points, D, negcurvature); 
      else
       nj_points_to_distances(points, D);

      tree = nj_infer_tree(D, msa->names);
      nj_reset_tree_model(mod, tree);      
      ll = nj_compute_log_likelihood(mod, msa, NULL);
      deriv = (ll - ll_base) / DERIV_EPS; 

      /* the partial derivative wrt the mean parameter is equal to the
         derivative with respect to the point, because the point is just a
         translation of a 0-mean MVN variable via the reparameterization
         trick */
      vec_set(grad, pidx, deriv);

      /* the partial derivative wrt the variance parameter, however,
         is more complicated, because of the reparameterization trick */

      /* first rederive the original standard normal rv */
      sd = sqrt(vec_get(sigmapar, pidx));
      stdrv = (porig - vec_get(mvn->mu, pidx)) / sd; 
      
      if (covar_param == DIAG)
        /* in the DIAG case, the partial derivative wrt the
           corresponding variance parameter can be computed directly
           based on a single point and coordinate */        
        vec_set(grad, (i+n)*d + k, deriv * 0.5 * stdrv / sd);
      
      else 
        /* in the DIST case, we add to a running total and update at the end */
        lambda_grad += deriv * 0.5 * stdrv / sqrtl * vec_get(data->chol_colsum, i);
      
      vec_set(points, pidx, porig); /* restore orig */
    }
  }
  if (covar_param == DIST) /* in this case, need to update the final
                              gradient component corresponding to the
                              lambda parameter */
    vec_set(grad, n*d, lambda_grad);

  nj_reset_tree_model(mod, orig_tree);
  return ll_base;
}  


/* optimize variational model by stochastic gradient ascent using the
   Adam algorithm.  Takes initial tree model and alignment and
   distance matrix, dimensionality of Euclidean space to work in.
   Note: alters distance matrix */
void nj_variational_inf(TreeModel *mod, MSA *msa, Matrix *D, MVN *mvn,
                        int dim, unsigned int hyperbolic, double negcurvature,
                        int nminibatch, double learnrate, int nbatches_conv,
                        int min_nbatches, Vector *sigmapar, enum covar_type covar_param,
                        CovarData *data, FILE *logf) {

  Vector *points, *grad, *avegrad, *m, *m_prev, *v, *v_prev, *best_mu, *best_sigmapar;
  int n = msa->nseqs, i, j, t, stop = FALSE, bestt = -1, graddim;
  double ll, avell, kld, avekld, bestelb = -INFTY, bestll = -INFTY, bestkld = -INFTY,
    running_tot = 0, last_running_tot = -INFTY;
  
  if (mvn->dim != n*dim)
    die("ERROR in nj_variational_inf: bad dimensions\n");
  
  points = vec_new(mvn->dim);
  graddim = (covar_param == DIST ? mvn->dim + 1 : 2*mvn->dim);
  grad = vec_new(graddim);  
  avegrad = vec_new(graddim);
  m = vec_new(graddim);
  m_prev = vec_new(graddim);
  v = vec_new(graddim);
  v_prev = vec_new(graddim);

  best_mu = vec_new(mvn->dim);
  vec_copy(best_mu, mvn->mu);
  best_sigmapar = vec_new(sigmapar->size);
  vec_copy(best_sigmapar, sigmapar);

  /* set up log file */
  if (logf != NULL) {
    fprintf(logf, "# nj_var logfile\n");
    fprintf(logf, "state\tll\tkld\telb\t");
    for (j = 0; j < mvn->dim; j++)
      fprintf(logf, "mu.%d\t", j);
    for (i = 0; i < D->nrows; i++)
      for (j = i+1; j < D->ncols; j++)
        fprintf(logf, "D.%d.%d\t", i, j);
    fprintf(logf, "\n");
  }
  
  /* initialize moments for Adam algorithm */
  vec_zero(m);  vec_zero(m_prev);
  vec_zero(v);  vec_zero(v_prev);
  t = 0;
  
  do {
    vec_zero(avegrad);
    avell = 0;
    avekld = 0;

    for (i = 0; i < nminibatch; i++) {
      /* sample points from MVN averaging distribution */
      mvn_sample(mvn, points);
      
      /* compute likelihood and gradient */
      ll = nj_compute_model_grad(mod, mvn, msa, hyperbolic, negcurvature, 
                                 points, grad, D, sigmapar, covar_param, data);

      /* compute the KLD (equation 7, Doersch arXiv 2016) */
      kld = 0;
      for (j = 0; j < mvn->dim; j++) {
        kld += 0.5 * (vec_get(sigmapar, j) + vec_get(mvn->mu, j) * vec_get(mvn->mu, j)); 
         /* 1/2 trace of sigma and inner product of mu with itself */

        if (vec_get(sigmapar, j) > 0)
          kld -= 0.5 * log(vec_get(sigmapar, j)); /* contribution to log determinant of sigma */
      }
      
      kld -= 0.5 * dim;  /* 1/2 of dimension */

      avell += ll;
      avekld += kld;
      
      /* incorporate gradient of KLD.  Must be subtracted rather than added */
      for (j = 0; j < grad->size; j++) {
        double gj = 0.0;

        if (j < n*dim)  /* partial deriv wrt mu_j is just mu_j */
          gj = -1.0*vec_get(mvn->mu, j);
        else {            /* partial deriv wrt sigma_j is more
                           complicated because of the trace and log
                           determinant */
          gj = 0.5 * (-1.0 + 1.0/vec_get(sigmapar, j-mvn->mu->size));   /* first term trace, second log det */
        }
        vec_set(grad, j, vec_get(grad, j) + gj);
      }

      /* add gradient to running total */
      vec_plus_eq(avegrad, grad); 
    }

    /* divide by nminibatch to get expected gradient */
    vec_scale(avegrad, 1.0/nminibatch);
    avell /= nminibatch;
    avekld /= nminibatch;

    /* store parameters if best yet */
    if (avell - avekld > bestelb) {
      bestelb = avell - avekld;
      bestll = avell;  /* not necessarily best ll but ll corresponding to bestelb */
      bestkld = avekld;  /* same comment */
      bestt = t;
      vec_copy(best_mu, mvn->mu);
      vec_copy(best_sigmapar, sigmapar);
    }
    
    /* Adam updates; see Kingma & Ba, arxiv 2014 */
    t++;
    for (j = 0; j < avegrad->size; j++) {   
      double mhatj, vhatj;      
      vec_set(m, j, ADAM_BETA1 * vec_get(m_prev, j) + (1.0 - ADAM_BETA1) * vec_get(avegrad, j));
      vec_set(v, j, ADAM_BETA2 * vec_get(v_prev, j) +
              (1.0 - ADAM_BETA2) * vec_get(avegrad, j) * vec_get(avegrad, j));
      mhatj = vec_get(m, j) / (1.0 - pow(ADAM_BETA1, t));
      vhatj = vec_get(v, j) / (1.0 - pow(ADAM_BETA2, t));

      /* update mu or sigma, depending on parameter index */
      if (j < n*dim)
        vec_set(mvn->mu, j, vec_get(mvn->mu, j) + learnrate * mhatj / (sqrt(vhatj) + ADAM_EPS)); 
      else {
        vec_set(sigmapar, j-mvn->dim, vec_get(sigmapar, j-mvn->dim) +
                learnrate * mhatj / (sqrt(vhatj) + ADAM_EPS)); 

        /* don't allow sigma to go negative */
        if (vec_get(sigmapar, j-mvn->dim) < MIN_VAR)
          vec_set(sigmapar, j-mvn->dim, MIN_VAR);
      }
    }
    nj_update_covariance(mvn, sigmapar, covar_param, data);
    
    vec_copy(m_prev, m);
    vec_copy(v_prev, v);
    
    /* report to log file */
    if (logf != NULL) {
      fprintf(logf, "%d\t%f\t%f\t%f\t", t, avell, avekld, avell - avekld);
      for (j = 0; j < mvn->mu->size; j++)
        fprintf(logf, "%f\t", vec_get(mvn->mu, j));
      if (hyperbolic)
        nj_points_to_distances_hyperbolic(mvn->mu, D, negcurvature);
      else
        nj_points_to_distances(mvn->mu, D);
      for (i = 0; i < D->nrows; i++)
        for (j = i+1; j < D->ncols; j++)
          fprintf(logf, "%f\t", mat_get(D, i, j));
      fprintf(logf, "\n");
    }
    
    /* check total elb every nbatches_conv to decide whether to stop */
    running_tot += avell - avekld;
    if (t % nbatches_conv == 0) {
      if (logf != NULL)
        fprintf(logf, "# Average for last %d: %f\n", nbatches_conv, running_tot/nbatches_conv);
      if (t >= min_nbatches && running_tot <= last_running_tot)
        stop = TRUE;
      else {
        last_running_tot = running_tot;
        running_tot = 0;
      }
    }
    
  } while(stop == FALSE);
    
  vec_copy(mvn->mu, best_mu);
  vec_copy(sigmapar, best_sigmapar);
  nj_update_covariance(mvn, sigmapar, covar_param, data);
  
  if (logf != NULL) {
    fprintf(logf, "# Reverting to parameters from iteration %d; ELB: %.2f, LNL: %.2f, KLD: %.2f, ",
            bestt+1, bestelb, bestll, bestkld);
     fprintf(logf, "mu:\t");
     vec_print(mvn->mu, logf);
  }
  
  vec_free(grad);
  vec_free(avegrad);
  vec_free(points);
  vec_free(m);
  vec_free(m_prev);
  vec_free(v);
  vec_free(v_prev);
  vec_free(best_mu);
  vec_free(best_sigmapar);
}

/* sample a list of trees from the approximate posterior distribution
   and return as a new list.  If logdens is non-null, return
   corresponding vector of log densities for the samples */
List *nj_var_sample(int nsamples, int dim, MVN *mvn, char** names,
                    unsigned int hyperbolic, double negcurvature, Vector *logdens) {
  List *retval = lst_new_ptr(nsamples);
  int i, n = mvn->dim / dim;
  Matrix *D = mat_new(n, n);
  TreeNode *tree;
  Vector *points = vec_new(mvn->dim);
  
  for (i = 0; i < nsamples; i++) {
     mvn_sample(mvn, points);

     if (logdens != NULL) 
       vec_set(logdens, i, mvn_log_dens(mvn, points));
     /* FIXME: need Jacobian in hyperbolic case */
     
     if (hyperbolic)
       nj_points_to_distances_hyperbolic(points, D, negcurvature);
     else
       nj_points_to_distances(points, D);

     tree = nj_infer_tree(D, names);
     lst_push_ptr(retval, tree);
  }
  
  mat_free(D);
  vec_free(points);
  return(retval);
}

/* return a single tree representing the approximate posterior mean */
TreeNode *nj_mean(Vector *mu, int dim, char **names, unsigned int hyperbolic,
                  double negcurvature) {
  int n = mu->size / dim;
  Matrix *D = mat_new(n, n);
  TreeNode *tree;
  
  if (n * dim != mu->size)
    die("ERROR in nj_mean: bad dimensions\n");

  if (hyperbolic)
    nj_points_to_distances_hyperbolic(mu, D, negcurvature);
  else
    nj_points_to_distances(mu, D);
  
  tree = nj_infer_tree(D, names);
  
  mat_free(D);
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

   
/* comparison function for sorting by eigenvalues (below) */
/** Structure representing complex number */
typedef struct {
  int idx;
  Vector *evals;		       
} Evidx;

int nj_eigen_compare_desc(const void* ptr1, const void* ptr2) {
  Evidx *idx1 = *((Evidx**)ptr1);
  Evidx *idx2 = *((Evidx**)ptr2);
  double eval1 = vec_get(idx1->evals, idx1->idx);
  double eval2 = vec_get(idx2->evals, idx2->idx);
  if (eval1 == eval2) return 0;
  else if (eval1 < eval2) return 1;
  return -1;
}

/* generate an approximate mu and sigma from a distance matrix, for
   use in initializing the variational inference algorithm.  */
void nj_estimate_mvn_from_distances(Matrix *D, int dim, MVN *mvn,
                                    Vector *sigmapar, enum covar_type covar_param,
                                    CovarData *data) {
  int n = D->nrows;
  Matrix *Dsq, *G, *revec_real;
  Vector *eval_real;
  int i, j, d, N;
  List *eiglst;
  double rowsum_orig = 0, rowsum_new = 0, x = 0, x2 = 0;
    
  if (D->nrows != D->ncols || mvn->dim != n * dim)
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

  /* build gram matrix */
  G = mat_new(n, n);
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      mat_set(G, i, j, mat_get(Dsq, i, 1) + mat_get(Dsq, 1, j) -
	      mat_get(Dsq, i, j));

  /* find eigendecomposition of G */
  eval_real = vec_new(n);
  revec_real = mat_new(n, n);
  mat_diagonalize_sym(G, eval_real, revec_real);
  
  /* sort eigenvalues from largest to smallest */
  /* (eigenvalues are now guaranteed to be returned from smallest to
     largest so the sorting step could be avoided but will leave it
     for now */
  eiglst = lst_new_ptr(n);
  for (i = 0; i < n; i++) {
    Evidx *obj = malloc(sizeof(Evidx));
    obj->idx = i;
    obj->evals = eval_real;
    lst_push_ptr(eiglst, obj);
  }
  lst_qsort(eiglst, nj_eigen_compare_desc);
  
  /* create a vector of points based on the first 'dim' eigenvalues */
  for (d = 0; d < dim; d++) {
    Evidx *obj = lst_get_ptr(eiglst, d);
    double evalsqrt = sqrt(vec_get(eval_real, obj->idx));

    /* product of evalsqrt and corresponding column of revec will define
       the dth component of each point */ 
    for (i = 0; i < n; i++)
      vec_set(mvn->mu, i*dim + d, evalsqrt * mat_get(revec_real, i, obj->idx));
  }  

  /* rescale the matrix to match the original distance matrix */
  /* use sum of first row to normalize */
  nj_points_to_distances(mvn->mu, Dsq);   /* reuse Dsq here, no longer needed */
  for (j = 1; j < n; j++) {
    rowsum_orig += mat_get(D, 0, j);
    rowsum_new += mat_get(Dsq, 0, j);
  }
  vec_scale(mvn->mu, rowsum_orig/rowsum_new);
  
  for (i = 0; i < n; i++)
    free((Evidx*)lst_get_ptr(eiglst, i));
  lst_free(eiglst);

  /* initialize sigma to the identity scaled by 1/n of the variance
     across pairwise distances */
  for (i = 0; i < n; i++) {
    for (j = i+1; j < n; j++) {
      x += mat_get(D, i, j);
      x2 += mat_get(D, i, j) * mat_get(D, i, j);
    }
  }
  N = n * (n-1)/2;

  vec_set_all(sigmapar, 1.0/N * (x2/N - x*x/(N*N)));
  if (covar_param == DIST) /* set up the Laplacian pseudoinverse */
    nj_laplacian_pinv(data);
  nj_update_covariance(mvn, sigmapar, covar_param, data);
  
  mat_free(Dsq);
  mat_free(G);
  vec_free(eval_real);
  mat_free(revec_real);
}

/* generate an approximate mu and sigma from a distance matrix, for
   use in initializing the variational inference algorithm. In this
   version, use the 'hydra' algorithm to solve the problem
   approximately in hyperbolic space (Keller-Ressel & Nargang,
   arXiv:1903.08977, 2019) */
void nj_estimate_mvn_from_distances_hyperbolic(Matrix *D, int dim, MVN *mvn,
                                               double negcurvature,
                                               Vector *sigmapar, enum covar_type covar_param,
                                               CovarData *data) {
  int n = D->nrows;
  Matrix *A, *revec_real;
  Vector *eval_real;
  int i, j, d, N;
  List *eiglst;
  double x = 0, x2 = 0;
    
  if (D->nrows != D->ncols || mvn->dim != n * dim)
    die("ERROR in nj_estimate_points_from_distances_hyperbolic: bad dimensions\n");
  
  /* build matrix A of transformed distances; note that D is upper
     triangular but A must be symmetric */
  A = mat_new(n, n);
  for (i = 0; i < n; i++) {
    mat_set(A, i, i, 1);
    for (j = i + 1; j < n; j++) {
      double a = cosh(sqrt(negcurvature) * mat_get(D, i, j));
      mat_set(A, i, j, a);
      mat_set(A, j, i, a);
    }
  }
  
  /* find eigendecomposition of A */
  eval_real = vec_new(n);
  revec_real = mat_new(n, n);
  mat_diagonalize_sym(A, eval_real, revec_real);
  
  /* sort eigenvalues from largest to smallest */ 
  eiglst = lst_new_ptr(n);
  for (i = 0; i < n; i++) {
    Evidx *obj = malloc(sizeof(Evidx));
    obj->idx = i;
    obj->evals = eval_real;
    lst_push_ptr(eiglst, obj);
  }
  lst_qsort(eiglst, nj_eigen_compare_desc);

  /* create a vector of points based on eigenvalues (n-dim+1) through
     n (1st dimension will be implicit) */
  for (d = 0; d < dim; d++) {
    int eigd = n - dim + d;
    Evidx *obj = lst_get_ptr(eiglst, eigd); 
    double ev = -vec_get(eval_real, obj->idx);
    if (ev < 0) ev = 0;

    /* product of evalsqrt and corresponding column of revec will define
       the dth component of each point */ 
    for (i = 0; i < n; i++)
      vec_set(mvn->mu, i*dim + d, sqrt(ev) * mat_get(revec_real, i, obj->idx));
  }  

  for (i = 0; i < n; i++)
    free((Evidx*)lst_get_ptr(eiglst, i));
  lst_free(eiglst);

  /* initialize sigma to the identity scaled by 1/n of the variance
     across pairwise distances */
  for (i = 0; i < n; i++) {
    for (j = i+1; j < n; j++) {
      x += mat_get(D, i, j);
      x2 += mat_get(D, i, j) * mat_get(D, i, j);
    }
  }
  N = n * (n-1)/2;

  vec_set_all(sigmapar, 1.0/N * (x2/N - x*x/(N*N)));
  if (covar_param == DIST) /* set up the Laplacian pseudoinverse */
    nj_laplacian_pinv(data);
  nj_update_covariance(mvn, sigmapar, covar_param, data);
  
  mat_free(A);
  vec_free(eval_real);
  mat_free(revec_real);
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
   makes use a scaling factor to avoid underflow with large trees.  If
   branchgrad is non-null, it will be populated with the gradient of
   the log likelihood with respect to the individual branches of the
   tree, in post-order.  Note that this function uses natural log
   rather than log2.  */
double nj_compute_log_likelihood(TreeModel *mod, MSA *msa, Vector *branchgrad) {

  int i, j, k, nodeidx, rcat = 0, tupleidx;
  int nstates = mod->rate_matrix->size;
  TreeNode *n;
  double total_prob = 0;
  List *traversal;
  double **pL = NULL, **pLbar = NULL;
  double log_scale;
  double scaling_threshold = DBL_MIN;
  double ll = 0;
  double tmp[nstates];
  MarkovMatrix *temp_subst_mat = NULL;
  
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

  if (branchgrad != NULL) {
    if (branchgrad->size != mod->tree->nnodes)
      die("ERROR in nj_compute_log_likelihood: size of branchgrad must equal number of nodes in tree\n");
    vec_zero(branchgrad);
    pLbar = smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++)
      pLbar[j] = smalloc((mod->tree->nnodes+1) * sizeof(double));
    temp_subst_mat = mm_new(nstates, msa->alphabet, CONTINUOUS);
  }
  
  tm_set_subst_matrices(mod);  /* just call this in all cases; we'll be tweaking the model a lot */

  /* get sequence index if not already there */
  if (mod->msa_seq_idx == NULL)
    tm_build_seq_idx(mod, msa);

  traversal = tr_postorder(mod->tree);
  for (tupleidx = 0; tupleidx < msa->ss->ntuples; tupleidx++) {
    log_scale = 0;
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
	for (i = 0; i < nstates; i++) {
	  double totl = 0, totr = 0;
	  for (j = 0; j < nstates; j++)
	    totl += pL[j][n->lchild->id] *
	      mm_get(lsubst_mat, i, j);
	  
	  for (k = 0; k < nstates; k++)
	    totr += pL[k][n->rchild->id] *
	      mm_get(rsubst_mat, i, k);
          
	  if (totl * totr < scaling_threshold) {
	    pL[i][n->id] = (totl / scaling_threshold) * totr;
	    log_scale -= log(scaling_threshold);
	  }
	  else {
	    pL[i][n->id] = totl * totr;
	  }
	}
      }
    }
  
    /* termination */
    total_prob = 0;
    for (i = 0; i < nstates; i++)
      total_prob += vec_get(mod->backgd_freqs, i) *
	pL[i][mod->tree->id] * mod->freqK[rcat];
    
    ll += (log(total_prob) - log_scale) * msa->ss->counts[tupleidx];

    assert(isfinite(ll));


    /* to compute gradients efficiently, need to make a second pass
       across the tree */
    if (branchgrad != NULL) {
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

	  /* FIXME: need to do scale factor here also; separate variable?  how to propagate correctly? */
	  
	  for (j = 0; j < nstates; j++) { /* parent state */
	    tmp[j] = 0;
	    for (k = 0; k < nstates; k++) { /* sibling state */
	      tmp[j] += pLbar[j][n->parent->id] *
		pL[k][sibling->id] * mm_get(sib_subst_mat, j, k);
	    }
	  }

	  for (i = 0; i < nstates; i++) { /* child state */
	    pLbar[i][n->id] = 0;
	    for (j = 0; j < nstates; j++) { /* parent state */
	      pLbar[i][n->id] +=
		tmp[j] * mm_get(par_subst_mat, j, i);
	    }
	  }
	}
      }

      /* now compute branchwise derivatives in postorder in a final pass */
      traversal = tr_postorder(mod->tree);
      for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
	TreeNode *par;
	double base_prob = 0, new_prob = 0;
	
	n = lst_get_ptr(traversal, nodeidx);
	par = n->parent;
	
	if (par == NULL)  /* should be the last one visited but let's be safe */
	  continue;
       
	/* recompute total probability based on current node, to
	   avoid numerical errors */
	for (i = 0; i < nstates; i++) /* parent state */
	  for (j = 0; j < nstates; j++) /* child state */
	    base_prob += pL[j][n->id] * pLbar[i][par->id] *
	      mm_get(mod->P[n->id][rcat], i, j);
       
	/* recompute subst_mat with perturbed branch length */
	mm_cpy(temp_subst_mat, mod->P[n->id][rcat]);
	mm_exp(mod->P[n->id][rcat], mod->rate_matrix, n->dparent + DERIV_EPS);
	
	for (i = 0; i < nstates; i++) /* parent state */
	  for (j = 0; j < nstates; j++) /* child state */
	    new_prob += pL[j][n->id] * pLbar[i][par->id] *
	      mm_get(mod->P[n->id][rcat], i, j);

	/* restore original subst_mat */
	mm_cpy(mod->P[n->id][rcat], temp_subst_mat);
	
	vec_set(branchgrad, nodeidx, (log(new_prob) - log(base_prob)) / DERIV_EPS);
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
    mm_free(temp_subst_mat);
  }

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
                           TreeModel *mod, MSA *msa, FILE *logf) {
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
    ll = nj_compute_log_likelihood(mod, msa, NULL);
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
    fprintf(logf, "# Importance sampling from %d eligible trees out of %d; maxlnl = %f; ave sampled lnl = %f\n",
            count, lst_size(trees), maxll, sampll/nsamples);

  vec_free(weights);
  vec_free(lls);
  
  return(retval);
}

/* define new vector of covariance parameters depending on parameterization type */
Vector *nj_new_sigma_params(int ntips, int dim, enum covar_type covar_param) {
  int npars = 1;
  Vector *retval;
  
  if (covar_param == DIAG)
    npars = ntips * dim;

  retval = vec_new(npars);

  return retval;
}

/* update covariance matrix based on the parameters and (optionally)
   an auxiliary data object. In DIST case, leaves underlying distance
   matrix and derived matrices unchanged */
void nj_update_covariance(MVN *mvn, Vector *sigma_params, 
                          enum covar_type covar_param, CovarData *data) {
  int i;

  mat_zero(mvn->sigma);
  if (covar_param == DIAG) {
    assert(sigma_params->size == mvn->dim);
    for (i = 0; i < mvn->dim; i++)
      mat_set(mvn->sigma, i, i, vec_get(sigma_params, i));
  }
  else { /* DIST case */
    int dim, n, d, j;

    assert(sigma_params->size == 1 && data != NULL);
    dim = mvn->dim / data->Lapl_pinv->nrows; /* dimensionality of embedding */
    n = data->Lapl_pinv->nrows;  /* number of taxa */    
    data->lambda = vec_get(sigma_params, 0);

    /* Lapl_inv is n x n but sigma is n*dim x n*dim; we need to
       project upward.  Covariance terms for unequal coordinate axes
       (dimensions) are all zero assuming standard Brownian motion */
    for (i = 0; i < n; i++) {   /* row of Lapl_inv [taxa] */
      for (j = 0; j < n; j++) { /* col of Lapl_inv [taxa] */
        for (d = 0; d < dim; d++) {  /* dimension for row and col
                                        (only care about case where
                                        they are the same) */
          int sig_row = i*dim + d, sig_col = j*dim + d; /* corresponding row and 
                                                           col of sigma */
          mat_set(mvn->sigma, sig_row, sig_col, mat_get(data->Lapl_pinv, i, j));
          /* cases where dimensions different implicitly left as zero */
        }
      }
    }
    mat_scale(mvn->sigma, data->lambda);
  }
  mvn_update_type(mvn);
}

/* create a new CovarData object of the desired dimension for use in
   the DIST parameterization of covariance */
CovarData *nj_new_covar_data(Matrix *dist) {
  CovarData *retval = smalloc(sizeof(CovarData));   /* check */
  retval->lambda = LAMBDA_INIT;
  retval->dist = dist;
  retval->Lapl_pinv = mat_new(dist->nrows, dist->ncols);
  retval->cholL = mat_new(dist->nrows, dist->ncols);
  retval->chol_colsum = vec_new(dist->ncols);
  return (retval);
}

void nj_dump_covar_data(CovarData *data, FILE *F) {
  fprintf(F, "CovarData\nlambda: %f\n", data->lambda);
  fprintf(F, "distance matrix:\n");
  mat_print(data->dist, F);
  fprintf(F, "Laplacian pseudoinverse:\n");
  mat_print(data->Lapl_pinv, F);
  fprintf(F, "Cholesky lower-triangular:\n");
  mat_print(data->cholL, F);
  fprintf(F, "Cholesky column sums:\n");
  vec_print(data->chol_colsum, F);
}

/* define Laplacian pseudoinverse from distance matrix, for use with
   DIST parameterization of covariance.  Also compute Cholesky
   decomposition for use in gradient calculations. Store everything in
   CovarData */
void nj_laplacian_pinv(CovarData *data) {
  int i, j, dim = data->dist->nrows;
  double grandmean = 0, mindiag = INFTY, val, x;
  Vector *row_mean = vec_new(dim);
  
  /* define Laplacian pseudoinverse as double centered version of
     distance matrix */
  
  /* first compute column/row means and grand mean */
  vec_zero(row_mean);
  for (i = 0; i < dim; i++) {
    double rowsum = 0;
    for (j = 0; j < dim; j++) {
      val = (j >= i ? mat_get(data->dist, i, j) : mat_get(data->dist, j, i));
      /* stored in upper triangular form */
      rowsum += val;
      grandmean += val;
    }
    vec_set(row_mean, i, rowsum/dim);
  }
  grandmean /= (dim * dim);

  /* now double center */
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      x = (j >= i ? mat_get(data->dist, i, j) : mat_get(data->dist, j, i));
      val = -0.5 * (x - vec_get(row_mean, i) - vec_get(row_mean, j) + grandmean);

      if (i == j)
        val += 1.0e-6; /* CHECK: avoids singular matrix */
      
      mat_set(data->Lapl_pinv, i, j, val);
      if (i == j && val < mindiag)
        mindiag = val;
    }
  }
      
  /* this matrix is only defined up to a constant shift.  For it
     to define a valid covariance matrix (up to a scale constant) all
     values on the main diagonal must be nonnegative.  We can shift it
     so that the minimum such element is zero, thereby fixing one
     point and defining the others relative to it */
  mat_add_const(data->Lapl_pinv, -mindiag);    
  
  /* finally recompute the Cholesky decomposition */
  /*  retval = mat_cholesky(data->cholL, data->Lapl_pinv);
  if (retval != 0)
    die("ERROR in nj_laplacian_pinv. Cannot compute Cholesky decomposition of Laplacian pseudoinverse.\n");
  */
  
  /* also recompute the column sums of the Cholesky matrix for use in
     gradient calculations */
  vec_zero(data->chol_colsum);
  for (j = 0; j < dim; j++)  /* column */
    for (i = 0; i < dim; i++)  /* row */
      vec_set(data->chol_colsum, j, mat_get(data->cholL, i, j));

  vec_free(row_mean);
}
