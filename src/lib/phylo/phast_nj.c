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

/* sample a vector from a multivariate normal distribution with mean
   mu and covariance sigma. */
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
      mat_set(D, i, j, sqrt(sum));
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
  double porig, ll_base, ll, deriv, stdrv, sd;
  TreeNode *tree, *orig_tree;   /* has to be rebuilt repeatedly; restore at end */
  
  if (mu->size != n*d || sigma->nrows != n*d || sigma->ncols != n*d ||
      grad->size != 2*mu->size)
    die("ERROR in nj_compute_model_grad: bad parameters\n");

  /* set up tree model and get baseline log likelihood */
  nj_points_to_distances(points, D);
  tree = nj_infer_tree(D, msa->names);
  orig_tree = tr_create_copy(tree);   /* restore at the end */
  nj_reset_tree_model(mod, tree);
  ll_base = tl_compute_log_likelihood(mod, msa, NULL, NULL, -1, NULL);
  
  /* Perturb each point and propagate perturbation through distance
     calculation, neighbor-joining reconstruction, and likelihood
     calculation on tree */
 
  for (i = 0; i < n; i++) {
    for (k = 0; k < d; k++) {
      int pidx = i*d + k;
      porig = vec_get(points, pidx);
      vec_set(points, pidx, porig + DERIV_EPS);
      nj_points_to_distances(points, D);
      tree = nj_infer_tree(D, msa->names);
      nj_reset_tree_model(mod, tree);      
      ll = tl_compute_log_likelihood(mod, msa, NULL, NULL, -1, NULL);
      deriv = (ll - ll_base) / DERIV_EPS; 

      /* the partial derivative wrt the mean parameter is equal to the
         derivative with respect to the point, because the point is just a
         translation of a 0-mean MVN variable via the reparameterization
         trick */
      vec_set(grad, pidx, deriv);

      /* the partial derivative wrt the variance parameter, however,
         has an additional factor */
      sd = sqrt(mat_get(sigma, pidx, pidx));
      stdrv = (porig - vec_get(mu, pidx)) / sd;  /* orig standard normal rv */
      vec_set(grad, (i+n)*d + k, deriv * 0.5 * stdrv / sd);
      
      vec_set(points, pidx, porig); /* restore orig */
    }
  }
  nj_reset_tree_model(mod, orig_tree);
  return ll_base;
}  


/* optimize variational model by stochastic gradient ascent using the
   Adam algorithm.  Takes initial tree model and alignment and
   distance matrix, dimensionality of Euclidean space to work in.
   Note: alters distance matrix */
void nj_variational_inf(TreeModel *mod, MSA *msa, Matrix *D, Vector *mu, Matrix *sigma,
                        int dim, int nminibatch, double learnrate, FILE *logf) {

  Vector *points, *grad, *avegrad, *m, *m_prev, *v, *v_prev;
  int n = msa->nseqs, i, j, t, stop = FALSE;
  double ll, avell, det, running_tot = 0, last_running_tot = -INFTY;
  
  /* TEMPORARY: make this a parameter */
  int nbatches_conv = 10;
  
  if (mu->size != n*dim || sigma->nrows != n*dim || sigma->ncols != n*dim)
    die("ERROR in nj_variational_inf: bad dimensions\n");
  
  points = vec_new(mu->size);
  grad = vec_new(2*mu->size);
  avegrad = vec_new(2*mu->size);
  m = vec_new(2*mu->size);
  m_prev = vec_new(2*mu->size);
  v = vec_new(2*mu->size);
  v_prev = vec_new(2*mu->size);
  
  /* initialize moments for Adam algorithm */
  vec_zero(m);  vec_zero(m_prev);
  vec_zero(v);  vec_zero(v_prev);
  t = 0;
  
  do {
    vec_zero(avegrad);
    avell = 0;

    for (i = 0; i < nminibatch; i++) {
      /* sample points from MVN averaging distribution */
      nj_sample_mvn(mu, sigma, points);
      
      /* compute likelihood and gradient */
      ll = nj_compute_model_grad(mod, mu, sigma, msa, points, grad, D);

      /* add terms for KLD (equation 7, Doersch arXiv 2016) */
      det = 1;
      for (j = 0; j < sigma->nrows; j++) {
        ll -= 0.5 * (mat_get(sigma, j, j) + vec_get(mu, j) * vec_get(mu, j)); 
         /* 1/2 trace of sigma and inner product of mu with itself*/
        
        det *= mat_get(sigma, j, j); /* also compute determinant of sigma */
      }
      
      ll += 0.5 * dim;  /* 1/2 of dimension */
      ll += 0.5 * log(det); /* 1/2 of log determinant of sigma */
      
      avell += ll;

      /* add gradient of KLD */
      for (j = 0; j < grad->size; j++) {
        double gj;

        if (j < n*dim)  /* partial deriv wrt mu_j is just -mu_j */
          gj = -1.0*vec_get(mu, j);
        else            /* partial deriv wrt sigma_j is more
                           complicated because of the log
                           determinant */
          gj = 0.5 * (-1.0 + 1.0/mat_get(sigma, j-mu->size, j-mu->size));  

        vec_set(grad, j, vec_get(grad, j) + gj);
      }

      /* add gradient to running total */
      vec_plus_eq(avegrad, grad); 
    }

    /* divide by nminibatch to get expected gradient */
    vec_scale(avegrad, 1.0/nminibatch);
    avell /= nminibatch;
    
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
        vec_set(mu, j, vec_get(mu, j) + learnrate * mhatj / (sqrt(vhatj) + ADAM_EPS)); 
        /*        vec_set(mu, j, vec_get(mu, j) + learnrate * vec_get(avegrad, j)); *//* FIXME: temporary */ 
      else {
        mat_set(sigma, j-mu->size, j-mu->size, mat_get(sigma, j-mu->size, j-mu->size) +
                learnrate * mhatj / (sqrt(vhatj) + ADAM_EPS)); 
        /* mat_set(sigma, j-mu->size, j-mu->size, mat_get(sigma, j-mu->size, j-mu->size) +
           learnrate * vec_get(avegrad, j)); */

        /* don't allow sigma to go negative */
        if (mat_get(sigma, j-mu->size, j-mu->size) < MIN_VAR)
          mat_set(sigma, j-mu->size, j-mu->size, MIN_VAR);
      }
    }
    vec_copy(m_prev, m);
    vec_copy(v_prev, v);
    
    /* report gradient and parameters to a log file */
    if (logf != NULL) {
      fprintf(logf, "***\nIteration %d: ELB = %f\n", t, avell);
      fprintf(logf, "mu:\t");
      vec_print(mu, logf);
      fprintf(logf, "sigma:\n");
      mat_print(sigma, logf);
      fprintf(logf, "gradient:\t");
      vec_print(avegrad, logf);
      fprintf(logf, "Adam m:\t");
      vec_print(m, logf);
      fprintf(logf, "Adam v:\t");
      vec_print(v, logf);
      tr_print(logf, nj_mean(mu, dim, msa->names), TRUE);
      fprintf(logf, "distance matrix:\n");
      nj_points_to_distances(mu, D);
      mat_print(D, logf);
    }

    /* check total likelihood every nbatches_conv to decide whether to stop */
    running_tot += avell;
    if (t % nbatches_conv == 0) {
      if (running_tot <= last_running_tot)
	stop = TRUE; /* FIXME: better to revert to previous parameters? */
      else {
	last_running_tot = running_tot;
	running_tot = 0;
      }
    }
    
  } while(stop == FALSE);
    

  vec_free(grad);
  vec_free(avegrad);
  vec_free(points);
  vec_free(m);
  vec_free(m_prev);
  vec_free(v);
  vec_free(v_prev);
}

/* sample a list of trees from the approximate posterior distribution
   and return as a new list */
List *nj_var_sample(int nsamples, int dim, Vector *mu, Matrix *sigma, char** names) {
  List *retval = lst_new_ptr(nsamples);
  int i, n = mu->size / dim;
  Matrix *D = mat_new(n, n);
  TreeNode *tree;
  Vector *points = vec_new(n*dim);
  
  if (n * dim != mu->size || mu->size != sigma->nrows || mu->size != sigma->ncols)
    die("ERROR in nj_var_sample: bad dimensions\n");
 
  for (i = 0; i < nsamples; i++) {
     nj_sample_mvn(mu, sigma, points);
     nj_points_to_distances(points, D);
     tree = nj_infer_tree(D, names);
     lst_push_ptr(retval, tree);
  }
  
  mat_free(D);
  vec_free(points);
  return(retval);
}

/* return a single tree representing the approximate posterior mean */
TreeNode *nj_mean(Vector *mu, int dim, char **names) {
  int n = mu->size / dim;
  Matrix *D = mat_new(n, n);
  TreeNode *tree;
  
  if (n * dim != mu->size)
    die("ERROR in nj_mean: bad dimensions\n");
  
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
     call to tm_reset_tree; we just need to force the substitution
     matrices to be recomputed */
  tm_set_subst_matrices(mod);
}

/* FIXME: create a struct that contains the index and a pointer to the eval vector */
   
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
  return (eval2 - eval1);
}

/* generate an approximate mu and sigma from a distance matrix, for
   use in initializing the variational inference algorithm.  */
void nj_estimate_mvn_from_distances(Matrix *D, int dim, Vector *mu, Matrix *sigma) {
  int n = D->nrows;
  Matrix *Dsq, *G, *revec_real;
  Zvector *eval;
  Zmatrix *revec, *levec;
  Vector *eval_real;
  int i, j, d, N;
  List *eiglst;
  double rowsum_orig = 0, rowsum_new = 0, x = 0, x2 = 0;
    
  if (D->nrows != D->ncols || mu->size != n * dim || sigma->nrows != mu->size ||
      sigma->ncols != mu->size)
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
  eval = zvec_new(n);
  revec = zmat_new(n, n);
  levec = zmat_new(n, n);
  mat_diagonalize(G, eval, revec, levec);
  
  /* convert eigenvalues and right eigenvectors to real numbers; will
     fail if they have imaginary component but they should not because
     G is symmetric by construction */
  eval_real = vec_new(n);
  revec_real = mat_new(n, n);
  zvec_as_real(eval_real, eval, TRUE);
  zmat_as_real(revec_real, revec, TRUE);
  
  /* sort eigenvalues from largest to smallest */   /* FIXME: only keep nonzero? */
  eiglst = lst_new_ptr(n);
  for (i = 0; i < n; i++) {
    Evidx *obj = malloc(sizeof(Evidx));
    obj->idx = i;
    obj->evals = eval_real;
    lst_push_ptr(eiglst, obj);
  }
  lst_qsort(eiglst, nj_eigen_compare_desc);

  /* FIXME: what to do with zero eigenvalues?  random numbers? */
  
  /* create a vector of points based on the first 'dim' eigenvalues */
  for (d = 0; d < dim; d++) {
    Evidx *obj = lst_get_ptr(eiglst, d);
    double evalsqrt = sqrt(vec_get(eval_real, obj->idx));

    /* product of evalsqrt and corresponding column of revec will define
       the dth component of each point */ 
    for (i = 0; i < n; i++)
      vec_set(mu, i*dim + d, evalsqrt * mat_get(revec_real, i, obj->idx));
  }  

  /* rescale the matrix to match the original distance matrix */
  /* use sum of first row to normalize */
  nj_points_to_distances(mu, Dsq);   /* reuse Dsq here, no longer needed */
  for (j = 1; j < n; j++) {
    rowsum_orig += mat_get(D, 0, j);
    rowsum_new += mat_get(Dsq, 0, j);
  }
  vec_scale(mu, rowsum_orig/rowsum_new);
  
  for (i = 0; i < n; i++)
    free((Evidx*)lst_get_ptr(eiglst, i));
  lst_free(eiglst);

  /* initialize sigma to the identity scaled by 1/n of the variance
     across pairwise distances */
  mat_set_identity(sigma);
  for (i = 0; i < n; i++) {
    for (j = i+1; j < n; j++) {
      x += mat_get(D, i, j);
      x2 += mat_get(D, i, j) * mat_get(D, i, j);
    }
  }
  N = n * (n-1)/2;
  mat_scale(sigma, 1.0/N * (x2/N - x*x/(N*N)));
  
  mat_free(Dsq);
  mat_free(G);
  zvec_free(eval);
  zmat_free(revec);
  zmat_free(levec);
  vec_free(eval_real);
  mat_free(revec_real);
}
