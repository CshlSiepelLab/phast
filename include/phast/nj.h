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
#include <phast/misc.h>

#define DERIV_EPS 1e-5

/* tuning parameters for Adam algorithm.  These will be kept at the
   default values.  The learning rate (called alpha) will be passed in
   as a parameter */
#define ADAM_BETA1 0.9
//#define ADAM_BETA2 0.999
#define ADAM_BETA2 0.9
#define ADAM_EPS 1e-8

/* don't allow variance terms to get smaller than this value */
#define MIN_VAR 1e-6

/* types of parameterization for covariance matrix: constant (and
   diagonal), diagonal with free variances, proportional to Laplacian
   pseudoinverse based on pairwise distances, or low-rank
   approximation to full matrix */
enum covar_type {CONST, DIAG, DIST, LOWR};
  
/* auxiliary data for parameterization of covariance matrix in DIST
   case */
#define LAMBDA_INIT 0.001
typedef struct {
  enum covar_type type; /* type of parameterization */
  int nseqs; /* number of taxa in tree */
  int dim; /* dimension of point embedding */
  enum mvn_type mvn_type;
  Vector *params; /* vector of free parameters */
  double lambda;  /* scale parameter for covariance matrix 
                     (DIST or CONST cases) */
  Matrix *dist;   /* distance matrix on which covariance is based */
  int lowrank;  /* dimension of low-rank approximation if LOWR or -1
                   otherwise */
  Matrix *Lapl_pinv;  /* Laplacian pseudoinverse (DIST) */
  Vector *Lapl_pinv_evals; /* eigendecomposition of Lapl_pinv (DIST) */
  Matrix *Lapl_pinv_evecs;
  Vector *Lapl_pinv_sqrt_evals; /* precompute for efficiency (DIST) */
  Matrix *R; /* used for LOWR; has dimension lowrank x nseqs */
  double sparsity; /* multiplier for sparsity penalty */
  double penalty; /* the current value of the penalty */
} CovarData;

void nj_resetQ(Matrix *Q, Matrix *D, Vector *active, Vector *sums, int *u,
	       int *v, int maxidx);

void nj_updateD(Matrix *D, int u, int v, int w, Vector *active, Vector *sums);

TreeNode* nj_infer_tree(Matrix *initD, char **names);

double nj_compute_JC_dist(MSA *msa, int i, int j);

Matrix *nj_compute_JC_matr(MSA *msa);

Matrix *nj_tree_to_distances(TreeNode *tree, char **names, int n);

double nj_distance_on_tree(TreeNode *root, TreeNode *n1, TreeNode *n2);

void nj_points_to_distances(Vector *points, Matrix *D,
                            double negcurvature, unsigned int use_hyperbolic);

void nj_points_to_distances_euclidean(Vector *points, Matrix *D);

void nj_points_to_distances_hyperbolic(Vector *points, Matrix *D,
                                       double negcurvature);

double nj_compute_model_grad(TreeModel *mod, multi_MVN *mmvn, MSA *msa,
                             unsigned int hyperbolic, double negcurvature,
                             Vector *points, Vector *points_std,
                             Vector *grad, Matrix *D,
                             CovarData *data);

double nj_compute_model_grad_check(TreeModel *mod, multi_MVN *mmvn, MSA *msa,
                                   unsigned int hyperbolic, double negcurvature,
                                   Vector *points, Vector *points_std,
                                   Vector *grad, Matrix *D,
                                   CovarData *data);

void nj_variational_inf(TreeModel *mod, MSA *msa, Matrix *D, multi_MVN *mmvn,
                        int dim, unsigned int hyperbolic, double negcurvature,
                        int nminibatch, double learnrate, int nbatches_conv,
                        int min_nbatches, CovarData *data, FILE *logf);

List *nj_var_sample(int nsamples, int dim, multi_MVN *mmvn,
                    char** names, unsigned int hyperbolic,
                    double negcurvature, Vector *logdens);

TreeNode *nj_mean(Vector *mu, int dim, char **names,
                  unsigned int hyperbolic, double negcurvature);

void nj_reset_tree_model(TreeModel *mod, TreeNode *newtree);

void nj_estimate_mmvn_from_distances(Matrix *D, int dim, multi_MVN *mmvn,
                                     double negcurvature, CovarData *data,
                                     unsigned int use_hyperbolic);

void nj_estimate_mmvn_from_distances_euclidean(Matrix *D, int dim, multi_MVN *mmvn,
                                               CovarData *data);

void nj_estimate_mmvn_from_distances_hyperbolic(Matrix *D, int dim, multi_MVN *mmvn,
                                               double negcurvature,
                                               CovarData *data);

void nj_test_D(Matrix *D);

double nj_compute_log_likelihood(TreeModel *mod, MSA *msa, Vector *branchgrad);

int *nj_build_seq_idx(List *leaves, char **names) ;

int nj_get_seq_idx(char **names, char *name, int n);

List *nj_importance_sample(int nsamples, List *trees, Vector *logdens,
                           TreeModel *mod, MSA *msa, FILE *logf);

void nj_update_covariance(multi_MVN *mmvn, CovarData *data);

CovarData *nj_new_covar_data(enum covar_type covar_param, Matrix *dist,
                             int dim, int rank, double sparsity);

void nj_dump_covar_data(CovarData *data, FILE *F);

void nj_laplacian_pinv(CovarData *data);

void nj_mmvn_to_distances(multi_MVN *mmvn, Matrix *D, unsigned int hyperbolic,
                          double negcurvature);

void nj_set_kld_grad_LOWR(Vector *kldgrad, multi_MVN *mvn);

void nj_rescale_grad(Vector *grad, Vector *rsgrad, multi_MVN *mmvn,
                     CovarData *data);

void nj_set_sparsity_penalty_LOWR(Vector *grad, multi_MVN *mmvn,
                                  CovarData *data);

#endif
