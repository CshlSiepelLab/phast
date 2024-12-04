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

#define DERIV_EPS 1e-5

/* tuning parameters for Adam algorithm.  These will be kept at the
   default values.  The learning rate (called alpha) will be passed in
   as a parameter */
#define ADAM_BETA1 0.9
#define ADAM_BETA2 0.999
#define ADAM_EPS 1e-8

/* don't allow variance terms to get smaller than this value */
#define MIN_VAR 1e-6


void nj_resetQ(Matrix *Q, Matrix *D, Vector *active, Vector *sums, int *u,
	       int *v, int maxidx);

void nj_updateD(Matrix *D, int u, int v, int w, Vector *active, Vector *sums);

TreeNode* nj_infer_tree(Matrix *initD, char **names);

double nj_compute_JC_dist(MSA *msa, int i, int j);

Matrix *nj_compute_JC_matr(MSA *msa);

/* void nj_sample_mvn(Vector *mu, Matrix *sigma, Vector *retval);*/

void nj_sample_std_mvn(Vector *retval);

void nj_sample_mvn(Vector *mu, Matrix *sigma, Vector *retval);

void nj_points_to_distances(Vector *points, Matrix *D);

/* TreeNode* nj_mvn_sample_tree(Vector *mu, Matrix *sigma, int n, char **names); */

double nj_compute_model_grad(TreeModel *mod, Vector *mu, Matrix *sigma, MSA *msa,
			     Vector *points, Vector *grad, Matrix *D);

void nj_variational_inf(TreeModel *mod, MSA *msa, Matrix *D, Vector *mu, Matrix *sigma,
                        int dim, int nminibatch, double learnrate, FILE *logf);

List *nj_var_sample(int nsamples, int dim, Vector *mu, Matrix *sigma,
                    char** names);

TreeNode *nj_mean(Vector *mu, int dim, char **names);

void nj_reset_tree_model(TreeModel *mod, TreeNode *newtree);

void nj_estimate_mvn_from_distances(Matrix *D, int dim, Vector *mu, Matrix *sigma);

void nj_test_D(Matrix *D);

#endif
