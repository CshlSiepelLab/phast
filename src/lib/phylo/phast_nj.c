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
/*
  void nj_sample_mvn(Vector *mu, Matrix *sigma, Vector *retval) {
  gsl_vector *gsl_mu = gsl_vector_alloc(mu->size);
  gsl_matrix *gsl_sigma = gsl_matrix_alloc(sigma->nrow, sigma->ncol);
  static gsl_rng *r = NULL;
  int i, j;

  if (r == NULL) 
    r = gsl_rng_alloc(gsl_rng_taus);
  
  for (i = 0; i < mu->size; i++)
    gsl_vec_set(gsl_mu, i, vec_get(mu, i));
  for (i = 0; i < sigma->nrow; i++)
    for (j = 0; j < sigma->ncol; j++)
      gsl_mat_set(gsl_sigma, i, j, mat_get(sigma, i, j));
 
  gsl_linalg_cholesky_decomp1(gsl_sigma);
  gsl_ran_multivariate_gaussian(r, gsl_mu, gsl_sigma, gsl_retval);

  for (i = 0; i < mu->size; i++)
    vec_set(retval, i, gsl_vector_get(gsl_retval, i));

  gsl_vector_free(gsl_mu);
  gsl_matrix_free(gsl_sigma);
}
*/

/* sample a vector from a standard multivariate normal distribution,
   with zero mean and identity covariance.  */
void nj_sample_std_mvn(Vector *retval) {
  int i;
  static int seeded = 0;
  double u1, u2, z1, z2;
  
  if (!seeded) {
    srandom((unsigned int)time(NULL));
    seeded = 1;
  }

  /* draw indep samples from standard normal using Box-Muller transform */
  for (i = 0; i < retval->size; i += 2) {
    u1 = unif_rand();
    u2 = unif_rand();
    z1 = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
    z2 = sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
    vec_set(retval, i, z1);
    if (i+1 < retval->size)
      vec_set(retval, i+1, z2);
  }  
}

/* sample a vector from a multivariate normal distribution with mean mu and covariance sigma. */
/* FIXME: for now this assumes diagonal sigma */
void nj_sample_mvn(Vector *mu, Matrix *sigma, Vector *retval) {
  int i;
  
  if (mu->size != sigma->nrows || mu->size != sigma->ncols || retval->size != mu->size)
    die("ERROR in nj_sample_mvn: bad dimension\n");
  
  nj_sample_std_mvn(retval);
  for (i = 0; i < retval->size; i++) {
    double v = vec_get(mu, i) + sqrt(mat_get(sigma, i, i)) * vec_get(retval, i);
    vec_set(retval, i, v);
  }
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
      mat_set(D, vidx1, vidx2, sqrt(sum));
    }
  }
}

/* compute the gradient of the log likelihood for a tree model with
   respect to the free parameters of the MVN averaging distribution,
   starting from a given MVN sample (points). Returns log likelihood
   of current model, which is computed as a by-product.  This version
   uses numerical methods */
double nj_compute_model_grad(TreeModel *mod, Vector *mu, Matrix *sigma, MSA *msa,
                           Vector *points, Vector *grad, Matrix *D) {
  int n = msa->nseqs;
  int d = mu->size / n;
  int i, k;
  double porig, ll_base, ll, deriv, vorig, stdrv;
  TreeNode *tree;
  
  if (mu->size != n*d || sigma->nrows != n || sigma->ncols != n ||
      grad->size != mu->size)
    die("ERROR in nj_compute_model_grad: bad parameters\n");

  /* set up tree model and get baseline log likelihood */
  /* FIXME: can we pass some of this in and avoid one call? */
  nj_points_to_distances(points, D);
  tree = nj_infer_tree(D, msa->names);
  tr_free(mod->tree);   /* FIXME: where first allocated? */
  mod->tree = tree;
  tm_set_subst_matrices(mod);
  ll_base = tl_compute_log_likelihood(mod, msa, NULL, NULL, -1, NULL);

  /* Perturb each point and propagate perturbation through distance
     calculation, neighbor-joining reconstruction, and likelihood
     calculation on tree */
 
  for (i = 0; i < n; i++) {
    for (k = 0; k < d; k++) {
      porig = vec_get(points, i*d + k);
      vec_set(points, i*d + k, porig + DERIV_EPS);
      nj_points_to_distances(points, D);
      tree = nj_infer_tree(D, msa->names);

      tr_free(mod->tree);   /* FIXME */
      mod->tree = tree;
      tm_set_subst_matrices(mod);  
      ll = tl_compute_log_likelihood(mod, msa, NULL, NULL, -1, NULL);
      
      deriv = (ll - ll_base) / DERIV_EPS; /* CHECK */

      /* the partial derivative wrt the mean parameter is equal to the
	 derivative with respect to the point, because the mean is just a
	 translation of a 0-mean MVN variable via the reparameterization
	 trick */
      vec_set(grad, i*d + k, deriv);

      /* the partial derivative wrt the variance parameter, however,
	 has an additional factor */

      vorig = mat_get(sigma, i*d + k, i*d + k);
      stdrv = (vorig - vec_get(mu, i*d + k)) / sqrt(mat_get(sigma, i*d + k, i*d + k));  
      vec_set(grad, i*d + k + n, deriv * 0.5 * sqrt(vorig) * stdrv);
      /* CHECK: need a factor here equal to the original
	 standardized random variate sampled; have I recomputed correctly? */
    }
  }
  return ll_base;
}  


/* optimize variational model by stochastic gradient ascent.  Takes
   initial tree model and alignment and distance matrix,
   dimensionality of Euclidean space to work in.  Note: alters
   distance matrix */

void nj_variational_inf(TreeModel *mod, MSA *msa, Matrix *D, int dim, int nminibatch, double learnrate) {

  Vector *mu, *points, *grad, *avegrad;
  Matrix *sigma;
  int n = msa->nseqs, i, stop = FALSE;
  double ll;
  
  mu = vec_new(n*dim);
  sigma = mat_new(n*dim, n*dim);
  points = vec_new(n*dim);
  grad = vec_new(n*dim*2);
  avegrad = vec_new(n*dim*2);
  
  /* create starting values of mu from distances.  can do from gram
     matrix but requires more code.  for now just initialize
     randomly. */
  nj_sample_std_mvn(mu);
  mat_set_identity(sigma);

  do {
    for (i = 0; i < nminibatch; i++) {
      /* sample points from MVN averaging distribution */
      nj_sample_mvn(mu, sigma, points);
      
      /* compute likelihood and gradient */
      ll = nj_compute_model_grad(mod, mu, sigma, msa, points, grad, D);

      /* add gradient to running total */
      vec_plus_eq(avegrad, grad);
      
      /* add terms from entropy */
    }

    /* divide by nminibatch to get expected gradient */
    vec_scale(avegrad, 1.0/nminibatch);
    
    /* FIXME report parameter values to a log file */

    /* update mu and sigma based on gradient and current learning rate */
    for (i = 0; i < n*dim; i++) {
      vec_set(mu, i, vec_get(mu, i) + learnrate * vec_get(avegrad, i));
      mat_set(sigma, i, i, mat_get(sigma, i, i) + learnrate * vec_get(avegrad, i + n*dim));
    }
    
    /* FIXME: what are stopping criteria? */    

  } while(!stop);
    

  vec_free(grad);
  vec_free(avegrad);
  vec_free(points);
  vec_free(mu);
  mat_free(sigma);
}


