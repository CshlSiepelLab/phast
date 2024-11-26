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

#define DERIV_EPS 1e-6

void nj_resetQ(Matrix *Q, Matrix *D, Vector *active, Vector *sums, int *u,
	       int *v, int maxidx);

void nj_updateD(Matrix *D, int u, int v, int w, Vector *active, Vector *sums);

TreeNode* nj_infer_tree(Matrix *initD, char **names);

double nj_compute_JC_dist(MSA *msa, int i, int j);

Matrix *nj_compute_JC_matr(MSA *msa);

/* void nj_sample_mvn(Vector *mu, Matrix *sigma, Vector *retval);*/

void nj_sample_std_mvn(Vector *retval);

void nj_points_to_distances(Vector *points, Matrix *D);

/* TreeNode* nj_mvn_sample_tree(Vector *mu, Matrix *sigma, int n, char **names); */

void nj_compute_model_grad(TreeModel *mod, Vector *mu, Matrix *sigma, MSA *msa,
                           Vector *points, Vector *grad);

#endif
