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

#define DERIV_EPS 1e-5

/* tuning parameters for Adam algorithm.  These will be kept at the
   default values.  The learning rate (called alpha) will be passed in
   as a parameter */
#define ADAM_BETA1 0.9
#define ADAM_BETA2 0.999
#define ADAM_EPS 1e-8

/* don't allow variance terms to get smaller than this value */
#define MIN_VAR 1e-6

/* types of parameterization for covariance matrix: diagonal only or
   version proportional to Laplacian pseudoinverse based on pairwise
   distances */
enum covar_type {DIAG, DIST};
  
/* auxiliary data for parameterization of covariance matrix in DIST
   case */
#define LAMBDA_INIT 0.1
typedef struct {
  double lambda;  /* scale parameter for covariance matrix */
  Matrix *dist;   /* distance matrix on which covariance is based */
  Matrix *Lapl_pinv;  /* Laplacian pseudoinverse */
  Vector *Lapl_pinv_evals; /* eigendecomposition of Lapl_pinv */
  Matrix *Lapl_pinv_evecs;
  Vector *Lapl_pinv_sqrt_evals; /* precompute for efficiency */
} CovarData;
  
/* bundle of MVNs sharing the same covariance matrix; for efficiency in DIST case */
typedef struct {
  int ntips; /* number of taxa */
  int d; /* number of dimensions, equal to number of component MVNs */
  Vector **mu; /* array of mean vectors */
  MVN *mvn; /* MVN object with shared covariance and associated preprocessing */
  enum covar_type type;
} multi_MVN;

void nj_resetQ(Matrix *Q, Matrix *D, Vector *active, Vector *sums, int *u,
	       int *v, int maxidx);

void nj_updateD(Matrix *D, int u, int v, int w, Vector *active, Vector *sums);

TreeNode* nj_infer_tree(Matrix *initD, char **names);

double nj_compute_JC_dist(MSA *msa, int i, int j);

Matrix *nj_compute_JC_matr(MSA *msa);

Matrix *nj_tree_to_distances(TreeNode *tree, char **names, int n);

double nj_distance_on_tree(TreeNode *root, TreeNode *n1, TreeNode *n2);

void nj_points_to_distances(Vector *points, Matrix *D);

void nj_points_to_distances_hyperbolic(Vector *points, Matrix *D,
                                       double negcurvature);

double nj_compute_model_grad(TreeModel *mod, multi_MVN *mmvn, MSA *msa,
                             unsigned int hyperbolic, double negcurvature,
                             Vector *points, Vector *grad, Matrix *D,
                             Vector *sigmapar, enum covar_type covar_param,
                             CovarData *data);

void nj_variational_inf(TreeModel *mod, MSA *msa, Matrix *D, multi_MVN *mmvn,
                        int dim, unsigned int hyperbolic, double negcurvature,
                        int nminibatch, double learnrate, int nbatches_conv,
                        int min_nbatches, Vector *sigmapar, enum covar_type covar_param,
                        CovarData *data, FILE *logf);

List *nj_var_sample(int nsamples, int dim, multi_MVN *mmvn,
                    char** names, unsigned int hyperbolic,
                    double negcurvature, Vector *logdens);

TreeNode *nj_mean(multi_MVN *mmvn, int dim, char **names,
                  unsigned int hyperbolic, double negcurvature);

void nj_reset_tree_model(TreeModel *mod, TreeNode *newtree);

void nj_estimate_mmvn_from_distances(Matrix *D, int dim, multi_MVN *mmvn,
                                    Vector *sigmapar,
                                    enum covar_type covar_param,
                                    CovarData *data);

void nj_estimate_mmvn_from_distances_hyperbolic(Matrix *D, int dim, multi_MVN *mmvn,
                                               double negcurvature,
                                               Vector *sigmapar,
                                               enum covar_type covar_param,
                                               CovarData *data);

void nj_test_D(Matrix *D);

double nj_compute_log_likelihood(TreeModel *mod, MSA *msa, Vector *branchgrad);

int *nj_build_seq_idx(List *leaves, char **names) ;

int nj_get_seq_idx(char **names, char *name, int n);

List *nj_importance_sample(int nsamples, List *trees, Vector *logdens,
                           TreeModel *mod, MSA *msa, FILE *logf);

Vector *nj_new_sigma_params(int ntips, int dim, enum covar_type covar_param);

void nj_update_covariance(multi_MVN *mmvn, Vector *sigma_params, 
                          enum covar_type covar_param, CovarData *data);

CovarData *nj_new_covar_data(Matrix *dist);

void nj_dump_covar_data(CovarData *data, FILE *F);

void nj_laplacian_pinv(CovarData *data);

multi_MVN *nj_multi_mvn_new(MVN *mvn, int d, enum covar_type type);

void nj_multi_mvn_sample(multi_MVN *mmvn, Vector *retval);

void nj_multi_mvn_combine_means(multi_MVN *mmvn, Vector *mu);

double nj_multi_mvn_log_dens(multi_MVN *mmvn, Vector *x);

double nj_multi_mvn_log_det(multi_MVN *mmvn);

void nj_multi_mvn_project_down(multi_MVN *mmvn, Vector *x_full,
                               Vector *x_d, int d);

void nj_multi_mvn_project_up(multi_MVN *mmvn, Vector *x_d,
                             Vector *x_full, int d);

void nj_multi_mvn_rederive_std(multi_MVN *mmvn, Vector *points,
                               Vector *points_std);

double nj_multi_mvn_trace(multi_MVN *mmvn);

double nj_multi_mvn_mu2(multi_MVN *mmvn);

void nj_multi_mvn_print(multi_MVN *mmvn, FILE *F);

#endif
