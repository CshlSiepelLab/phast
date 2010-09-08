/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: indel_mod.h,v 1.4 2008-11-12 02:07:59 acs Exp $ */

/* simple model of insertions and deletions, assumes given indel history */

#ifndef INDMOD_H
#define INDMOD_H

#include <indel_history.h>

#define NINDEL_STATES 3
//need to make sure ERROR isn't already defined from R libraries
#undef ERROR
typedef enum {MATCH, CHILDINS, CHILDDEL, SKIP, ERROR} col_type;

typedef struct {
  double alpha;
  double beta;
  double tau;
  double t;
  MarkovMatrix *probs;
  Matrix *log_probs;
  Vector *beg_probs;
  Vector *beg_log_probs;
} BranchIndelModel;

typedef struct {
  double alpha;
  double beta;
  double tau;
  double training_lnl;
  TreeNode *tree;
  BranchIndelModel **branch_mods;
} IndelModel;

typedef struct {
  Matrix *trans_counts;
  Vector *beg_counts;
} BranchIndelSuffStats;

typedef struct {
  TreeNode *tree;
  BranchIndelSuffStats **branch_counts;
} IndelSuffStats;

BranchIndelModel *im_new_branch(double alpha, double beta, double tau,
                                double t);
IndelModel *im_new_all(double alpha, double beta, double tau,  TreeNode *tree);
IndelModel *im_new(double *alpha, double *beta, double *tau, 
                   TreeNode *tree);
void im_free_branch(BranchIndelModel *bim);
void im_free(IndelModel *im);
void im_set_branch(BranchIndelModel *bim, double alpha, double beta, 
                   double tau, double t);
void im_set_all(IndelModel *im, double alpha, double beta, double tau, 
                TreeNode *tree);
void im_set(IndelModel *im, double *alpha, double *beta, double *tau, 
            TreeNode *tree);
double im_branch_column_logl(IndelHistory *ih, BranchIndelModel *bim, 
                             int child_id, double *col_logl);
double im_column_logl(IndelHistory *ih, IndelModel *im, double *col_logl);
double im_likelihood(IndelModel *im, IndelSuffStats *iss);
BranchIndelSuffStats *im_suff_stats_branch(IndelHistory *ih, int child_id);
IndelSuffStats *im_suff_stats(IndelHistory *ih);
void im_free_suff_stats(IndelSuffStats *iss);
double im_simulate_history(IndelModel *tim, int ncols);
void im_estimate(IndelModel *im, IndelHistory *ih, IndelSuffStats *ss, 
                 FILE *logf);
BranchIndelSuffStats *im_suff_stats_branch_cat(IndelHistory *ih, int child_id,
                                               int *categories, int do_cat);
IndelSuffStats *im_suff_stats_cat(IndelHistory *ih, int *categories, 
                                  int do_cat);

#endif
