/* $Id: fit_column.h,v 1.3 2008-02-18 22:03:36 acs Exp $
   Written by Adam Siepel, 2008 */

#ifndef FIT_COL_H
#define FIT_COL_H

#include <misc.h>
#include <tree_model.h>
#include <msa.h>
#include <vector.h>
#include <matrix.h>
#include <complex_matrix.h>

typedef enum {ALL, SUBTREE} scale_type;
typedef enum {CON, ACC, NNEUT} mode_type;

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
  Zvector *expdiag;

  double ***fels_scratch;       /* scratch memory for Felsenstein's
                                   alg (likelihoods and derivatives) */
  int nfels_scratch;            /* number of scratch arrays (depends on mode) */
  Zmatrix *mat_scratch;         /* scratch memory for derivative computation */
  Zvector *vec_scratch1, *vec_scratch2;
} ColFitData;

double col_compute_log_likelihood(TreeModel *mod, MSA *msa, int tupleidx,
                                  double **scratch);

double col_scale_derivs(ColFitData *d, double *first_deriv, 
                        double *second_deriv, double ***scratch);

double col_scale_derivs_subtree(ColFitData *d, Vector *gradient, 
                                Matrix *hessian, double ***scratch);

double col_likelihood_wrapper(Vector *params, void *data);

void col_grad_wrapper(Vector *grad, Vector *params, void *data, 
                      Vector *lb, Vector *ub);

void col_lrts(TreeModel *mod, MSA *msa, mode_type mode, double *tuple_pvals, 
              double *tuple_scales, double *tuple_llrs);

void col_lrts_sub(TreeModel *mod, MSA *msa, mode_type mode, 
                  double *tuple_pvals, double *tuple_null_scales, 
                  double *tuple_scales, double *tuple_sub_scales, 
                  double *tuple_llrs);

void col_score_tests(TreeModel *mod, MSA *msa, double *tuple_pvals, 
                     double *tuple_derivs, double *tuple_teststats);

void col_score_tests_sub(TreeModel *mod, MSA *msa, double *tuple_pvals, 
                         double *tuple_null_scales, double *tuple_derivs,
                         double *tuple_sub_derivs, double *tuple_teststats);

ColFitData *col_init_fit_data(TreeModel *mod, MSA *msa, scale_type stype,
                              mode_type mode, int second_derivs);

void col_free_fit_data(ColFitData *d);

#endif
