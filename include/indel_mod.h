/* $Id: indel_mod.h,v 1.2 2005-09-06 07:05:01 acs Exp $
   Written by Adam Siepel, 2005
   Copyright 2005, Adam Siepel, University of California */

/* simple model of insertions and deletions, assumes given indel history */

#ifndef INDMOD_H
#define INDMOD_H

#include <indel_history.h>

#define NINDEL_STATES 3
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
IndelModel *im_new(double alpha, double beta, double tau,  TreeNode *tree);
void im_free_branch(BranchIndelModel *bim);
void im_free(IndelModel *im);
void im_set_branch(BranchIndelModel *bim, double alpha, double beta, 
                   double tau, double t);
void im_set(IndelModel *im, double alpha, double beta, double tau, 
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
