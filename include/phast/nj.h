/* PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file nj.h
    Simple neighbor-joining tree inference  
    @ingroup phylo
*/

#ifndef NJ_H
#define NJ_H

#include <stdio.h>
#include <phast/matrix.h>
#include <phast/msa.h>
#include <phast/trees.h>
#include <phast/tree_model.h>
#include <phast/mvn.h>
#include <phast/multi_mvn.h>
#include <phast/sparse_matrix.h>
#include <phast/crispr.h>
#include <phast/misc.h>

/* for numerical derivatives */
#define DERIV_EPS 1e-5

/* tuning parameters for Adam algorithm.  These will be kept at the
   default values.  The learning rate (called alpha) will be passed in
   as a parameter */
#define ADAM_BETA1 0.9
//#define ADAM_BETA2 0.999
#define ADAM_BETA2 0.9
#define ADAM_EPS 1e-8

/* rescale embedding space so that maximum distances are equal to
   these values in the Euclidean and hyperbolic cases,
   respectively. Helps address the problem that branch lengths tend to
   be small so means and variances can get close to zero.  But scaling
   needs to be different for the two geometries */
#define POINTSPAN_EUC 50
#define POINTSPAN_HYP 3

/* use this as a floor for variance parameters.  Avoids drift to ever
   smaller values */
#define VARFLOOR 1.0e-3

/* initialization of lambda, which is scale factor for covariance
   matrix in DIST and CONST parameterizations */
#define LAMBDA_INIT 1.0e-5

/* types of parameterization for covariance matrix: constant (and
   diagonal), diagonal with free variances, proportional to Laplacian
   pseudoinverse based on pairwise distances, or low-rank
   approximation to full matrix */
enum covar_type {CONST, DIAG, DIST, LOWR};
  
/* auxiliary data for parameterization of covariance matrix in DIST
   case */
typedef struct {
  enum covar_type type; /* type of parameterization */
  int nseqs; /* number of taxa in tree */
  int dim; /* dimension of point embedding */
  enum mvn_type mvn_type;
  Vector *params; /* vector of free parameters */
  double lambda;  /* scale parameter for covariance matrix 
                     (DIST or CONST cases) */
  double pointscale; /* scale factor for geometry */
  unsigned int natural_grad; /* whether to rescale for natural
                                gradients during optimization */
  double kld_upweight; /* optional upweighting factor for KLD in ELBO */
  Matrix *dist;   /* distance matrix on which covariance is based */
  int lowrank;  /* dimension of low-rank approximation if LOWR or -1
                   otherwise */
  Matrix *Lapl_pinv;  /* Laplacian pseudoinverse (DIST) */
  Vector *Lapl_pinv_evals; /* eigendecomposition of Lapl_pinv (DIST) */
  Matrix *Lapl_pinv_evecs;
  Matrix *R; /* used for LOWR; has dimension lowrank x nseqs */
  double sparsity; /* multiplier for sparsity penalty */
  double penalty; /* the current value of the penalty */
  unsigned int hyperbolic; /* whether or not hyperbolic geometry is used */
  double negcurvature; /* for hyperbolic case */
  MSA *msa; /* multiple alignment under analysis if available */
  CrisprMutTable *crispr_muts; /* alternative CRISPR mutation table */
  char **names;
} CovarData;

/* for use with min-heap in fast nj algorithm */
typedef struct NJHeapData {
  double val;
  int i, j; 
  int rev_i, rev_j; // for lazy validation
} NJHeapNode;

void nj_resetQ(Matrix *Q, Matrix *D, Vector *active, Vector *sums, int *u,
	       int *v, int maxidx);

void nj_updateD(Matrix *D, int u, int v, int w, Vector *active, Vector *sums);

TreeNode* nj_infer_tree(Matrix *initD, char **names, Matrix *dt_dD);

TreeNode* nj_fast_infer(Matrix *initD, char **names, Matrix *dt_dD);

NJHeapNode* nj_heap_computeQ(int i, int j, int n, Matrix *D,
                             Vector *sums, int *rev);

double nj_compute_JC_dist(MSA *msa, int i, int j);

Matrix *nj_compute_JC_matr(MSA *msa);

Matrix *nj_tree_to_distances(TreeNode *tree, char **names, int n);

double nj_distance_on_tree(TreeNode *root, TreeNode *n1, TreeNode *n2);

void nj_points_to_distances(Vector *points, CovarData *data);

void nj_points_to_distances_euclidean(Vector *points, CovarData *data);

void nj_points_to_distances_hyperbolic(Vector *points, CovarData *data);

double nj_compute_model_grad(TreeModel *mod, multi_MVN *mmvn, 
                             Vector *points, Vector *points_std,
                             Vector *grad, CovarData *data);

double nj_compute_model_grad_check(TreeModel *mod, multi_MVN *mmvn, 
                                   Vector *points, Vector *points_std,
                                   Vector *grad, CovarData *data);

void nj_variational_inf(TreeModel *mod, multi_MVN *mmvn,
                        int nminibatch, double learnrate, int nbatches_conv,
                        int min_nbatches, CovarData *data, FILE *logf);

List *nj_var_sample(int nsamples, multi_MVN *mmvn, CovarData *data,
                    char** names, Vector *logdens);

TreeNode *nj_mean(Vector *mu, char **names, CovarData *data);

void nj_reset_tree_model(TreeModel *mod, TreeNode *newtree);

void nj_estimate_mmvn_from_distances(CovarData *data, multi_MVN *mmvn);

void nj_estimate_mmvn_from_distances_euclidean(CovarData *data, multi_MVN *mmvn);

void nj_estimate_mmvn_from_distances_hyperbolic(CovarData *data, multi_MVN *mmvn);

void nj_test_D(Matrix *D);

double nj_compute_log_likelihood(TreeModel *mod, MSA *msa, Vector *branchgrad);

int *nj_build_seq_idx(List *leaves, char **names) ;

int nj_get_seq_idx(char **names, char *name, int n);

List *nj_importance_sample(int nsamples, List *trees, Vector *logdens,
                           TreeModel *mod, MSA *msa, FILE *logf);

void nj_update_covariance(multi_MVN *mmvn, CovarData *data);

CovarData *nj_new_covar_data(enum covar_type covar_param, Matrix *dist,
                             int dim, MSA *msa, CrisprMutTable *crispr_muts,
                             char **names, unsigned int natural_grad,
                             double kld_upweight, int rank,
                             double sparsity, unsigned int hyperbolic,
                             double negcurvature);

void nj_dump_covar_data(CovarData *data, FILE *F);

void nj_laplacian_pinv(CovarData *data);

void nj_mmvn_to_distances(multi_MVN *mmvn, CovarData *data);

void nj_set_kld_grad_LOWR(Vector *kldgrad, multi_MVN *mvn);

void nj_rescale_grad(Vector *grad, Vector *rsgrad, multi_MVN *mmvn,
                     CovarData *data);

void nj_set_sparsity_penalty_LOWR(Vector *grad, multi_MVN *mmvn,
                                  CovarData *data);

void nj_set_LASSO_penalty_LOWR(Vector *grad, multi_MVN *mmvn,
                               CovarData *data);

void nj_set_pointscale(CovarData *data);

List *nj_var_sample_rejection(int nsamples, multi_MVN *mmvn,
                              CovarData *data, TreeModel *mod,
                              FILE *logf);

double nj_dL_dx_dumb(Vector *x, Vector *dL_dx, TreeModel *mod, 
                     CovarData *data);

double nj_dL_dt_num(Vector *dL_dt, TreeModel *mod, MSA *msa);

void nj_dt_dD_num(Matrix *dt_dD, Matrix *D, TreeModel *mod, MSA *msa);

int nj_i_j_to_dist(int i, int j, int n);

void nj_dist_to_i_j(int pwidx, int *i, int *j, int n);

double nj_dL_dx_smartest(Vector *x, Vector *dL_dx, TreeModel *mod, 
                         CovarData *data);


void nj_backprop(double *Jk, double *Jnext, int n, int f, int g, int u,
                 Vector *active);

void nj_backprop_sparse(SparseMatrix *Jk, SparseMatrix *Jnext, int n, int f, int g, int u,
                        Vector *active);

void nj_backprop_init(double *Jk, int n);

void nj_backprop_init_sparse(SparseMatrix *Jk, int n);

void nj_backprop_set_dt_dD(double *Jk, Matrix *dt_dD, int n, int f, int g,
                           int branch_idx_f, int branch_idx_g, Vector *active);

void nj_backprop_set_dt_dD_sparse(SparseMatrix *Jk, Matrix *dt_dD, int n, int f, int g,
                                  int branch_idx_f, int branch_idx_g, Vector *active);

#endif
