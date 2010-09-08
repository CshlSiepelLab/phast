/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#ifndef SUBST_DISTRIB
#define SUBST_DISTRIB

#include <vector.h>
#include <msa.h>
#include <tree_model.h>
#include <prob_vector.h>
#include <fit_column.h>

typedef struct {
  int njumps_max;
  double lambda;
  Matrix *R;
  TreeModel *mod;               
  Matrix **A;
  Matrix ***B;
  Matrix *M;
  Matrix ***branch_distrib;
  double epsilon;
} JumpProcess;
/* note: a jump process is defined wrt a whole tree model, not just a
   rate matrix */

typedef struct {
  double prior_mean, prior_var, post_mean, post_var, p_cons, p_anti_cons;
  int prior_min, prior_max, post_min, post_max;
} p_value_stats;

typedef struct {
  double prior_mean_left, prior_mean_right, prior_var_left, prior_var_right, 
    post_mean_left, post_mean_right, post_mean_tot, post_var_left, 
    post_var_right, post_var_tot, p_cons_right, p_cons_left, p_anti_cons_left, 
    p_anti_cons_right, cond_p_cons_left, cond_p_anti_cons_left,
    cond_p_cons_right, cond_p_anti_cons_right;
  int prior_min_left, prior_min_right, prior_max_left, prior_max_right, 
    post_min_left, post_min_right, post_min_tot, post_max_left, post_max_right,
    post_max_tot, cond_p_approx;
} p_value_joint_stats;

JumpProcess *sub_define_jump_process(TreeModel *mod, double epsilon, 
                                     double maxbranch);
void sub_free_jump_process(JumpProcess *jp);
Vector *sub_distrib_branch(JumpProcess *jp, double t);
Matrix **sub_distrib_branch_conditional(JumpProcess *jp, double t);
Vector *sub_prior_distrib_site(JumpProcess *jp);
Vector *sub_posterior_distrib_site(JumpProcess *jp, MSA *msa, int tuple_idx);
Vector *sub_prior_distrib_alignment(JumpProcess *jp, int nsites);
Vector *sub_posterior_distrib_alignment(JumpProcess *jp, MSA *msa);
void sub_pval_per_site(JumpProcess *jp, MSA *msa, mode_type side,
                       int fit_model, double *prior_mean, double *prior_var, 
                       double *pvals, double *post_mean, double *post_var,
                       FILE *logf);
void sub_pval_per_site_subtree(JumpProcess *jp, MSA *msa, mode_type mode, 
                               int fit_model, 
                               double *prior_mean_sub, double *prior_var_sub, 
                               double *prior_mean_sup, double *prior_var_sup, 
                               double *pvals, 
                               double *post_mean_sub, double *post_var_sub, 
                               double *post_mean_sup, double *post_var_sup, 
                               FILE *logf);
void sub_posterior_stats_alignment(JumpProcess *jp, MSA *msa, 
                                   double *mean, double *variance);
Matrix *sub_joint_distrib_site(JumpProcess *jp, MSA *msa, int tuple_idx);
Matrix *sub_prior_joint_distrib_alignment(JumpProcess *jp, int nsites);
Matrix *sub_posterior_joint_distrib_alignment(JumpProcess *jp, MSA *msa);
void sub_posterior_joint_stats_alignment(JumpProcess *jp, MSA *msa, 
                                         double *mean_tot, double *var_tot,
                                         double *mean_left, double *var_left,
                                         double *mean_right, double *var_right);
p_value_stats *sub_p_value_many(JumpProcess *jp, MSA *msa, List *feats, 
                                double ci);
p_value_joint_stats* sub_p_value_joint_many(JumpProcess *jp, MSA *msa, 
                                            List *feats, double ci,
                                            int max_convolve_size,
                                            FILE *timing_f);
void sub_reroot(TreeModel *mod, char *subtree_name);
#endif
