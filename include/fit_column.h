/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: fit_column.h,v 1.15 2008-12-16 19:54:09 mt269 Exp $ */

#ifndef FIT_COL_H
#define FIT_COL_H

#include <misc.h>
#include <tree_model.h>
#include <msa.h>
#include <vector.h>
#include <matrix.h>
#include <complex_matrix.h>

typedef enum {ALL, SUBTREE} scale_type;
typedef enum {CON, ACC, NNEUT, CONACC} mode_type;

/* metadata for fitting scale factors to individual alignment columns */
typedef struct {
  TreeModel *mod;
  MSA *msa;
  int tupleidx;
  scale_type stype;             /* whether doing all-branches or
                                   supertree/subtree estimation  */
  mode_type mode;               /* type of parameter bounding  */
  int second_derivs;            /* whether or not second derivatives
                                   need to be computed */
  Vector *params;
  Vector *lb;
  Vector *ub;
  double init_scale;
  double init_scale_sub;
  Matrix ***PP;
  Matrix ***PPP;
  Matrix ***QQ;
  Matrix ***QQQ;
  Matrix ***RRR;
  Zvector *expdiag_z;
  Vector *expdiag_r;

  double ***fels_scratch;       /* scratch memory for Felsenstein's
                                   alg (likelihoods and derivatives) */
  int nfels_scratch;            /* number of scratch arrays (depends on mode) */
  Zmatrix *mat_scratch_z;         /* scratch memory for derivative computation */
  Zvector *vec_scratch1_z, *vec_scratch2_z;
  Vector *vec_scratch1_r, *vec_scratch2_r;
  double deriv2;                /* second deriv for 1d case */
} ColFitData;

/* data for grid of precomputed Fisher Information Matrices */
#define GRIDSIZE1 0.02
#define GRIDSIZE2 0.05
#define GRIDMAXLOG 3
typedef struct {
  double *scales;               /* scale factors in grid */
  int ngrid1;                   /* number of grid points between 0 and
                                   1 (linear) */
  int ngrid2;                   /* number of grid points above 1 (log
                                   linear) */
  int ngrid;                    /* total number of grid points */
  Matrix **fim;                 /* precomputed FIMs */
} FimGrid;

double col_compute_log_likelihood(TreeModel *mod, MSA *msa, int tupleidx,
                                  double **scratch);

double col_compute_likelihood(TreeModel *mod, MSA *msa, int tupleidx,
		              double **scratch);

double col_scale_derivs(ColFitData *d, double *first_deriv, 
                        double *second_deriv, double ***scratch);

double col_scale_derivs_subtree(ColFitData *d, Vector *gradient, 
                                Matrix *hessian, double ***scratch);

double col_likelihood_wrapper(Vector *params, void *data);

double col_likelihood_wrapper_1d(double x, void *data);

void col_grad_wrapper(Vector *grad, Vector *params, void *data, 
                      Vector *lb, Vector *ub);

void col_lrts(TreeModel *mod, MSA *msa, mode_type mode, double *tuple_pvals, 
              double *tuple_scales, double *tuple_llrs, FILE *logf);

void col_lrts_sub(TreeModel *mod, MSA *msa, mode_type mode, 
                  double *tuple_pvals, double *tuple_null_scales, 
                  double *tuple_scales, double *tuple_sub_scales, 
                  double *tuple_llrs, FILE *logf);

void col_score_tests(TreeModel *mod, MSA *msa, mode_type mode, 
                     double *tuple_pvals, double *tuple_derivs, 
                     double *tuple_teststats);

void col_score_tests_sub(TreeModel *mod, MSA *msa, mode_type mode,
                         double *tuple_pvals, double *tuple_null_scales, 
                         double *tuple_derivs, double *tuple_sub_derivs, 
                         double *tuple_teststats, FILE *logf);

ColFitData *col_init_fit_data(TreeModel *mod, MSA *msa, scale_type stype,
                              mode_type mode, int second_derivs);

void col_free_fit_data(ColFitData *d);

void col_gerp(TreeModel *mod, MSA *msa, mode_type mode, double *tuple_nneut, 
              double *tuple_nobs, double *tuple_nrejected, 
              double *tuple_nspecies, FILE *logf);

void col_find_missing_branches(TreeModel *mod, MSA *msa, int tupleidx, 
                               int *has_data, int *nspec);

void col_scale_derivs_num(ColFitData *d, double *first_deriv, 
                          double *second_deriv);

void col_scale_derivs_subtree_num(ColFitData *d, Vector *gradient, 
                                  Matrix *hessian);

Matrix *col_estimate_fim_sub(TreeModel *mod);

FimGrid *col_fim_grid_sub(TreeModel *mod);

void col_free_fim_grid(FimGrid *g);

double col_estimate_fim(TreeModel *mod);

Matrix *col_get_fim_sub(FimGrid *g, double scale);

int col_has_data(TreeModel *mod, MSA *msa, int tupleidx);
int col_has_data_sub(TreeModel *mod, MSA *msa, int tupleidx, List *inside,
		                     List *outside);
#endif
