/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* simple model of insertions and deletions including a parameter for indels in
   motifs, assumes given indel history */

#ifndef DM_INDMOD_H
#define DM_INDMOD_H

#include <indel_history.h>
#include <msa.h>

#define NINDEL_STATES 3
typedef enum {MATCH, CHILDINS, CHILDDEL, SKIP, ERROR} col_type;

typedef struct {
  double alpha;
  double beta;
  double tau;
  double epsilon;
  double t;
  MarkovMatrix *probs;
  Matrix *log_probs;
  Vector *beg_probs;
  Vector *beg_log_probs;
} DMotifBranchIndelModel;

typedef struct {
  double alpha;
  double beta;
  double tau;
  double epsilon;
  double training_lnl;
  TreeNode *tree;
  DMotifBranchIndelModel **branch_mods;
} DMotifIndelModel;

typedef struct {
  Matrix *trans_counts;
  Vector *beg_counts;
} DMotifBranchIndelSuffStats;

typedef struct {
  TreeNode *tree;
  DMotifBranchIndelSuffStats **branch_counts;
} DMotifIndelSuffStats;

typedef enum {NEUT, CONS, DEATH, BIRTH} dmevent_t; /* Declared as an 
						      incomplete type in
						      dmotif_phmm.h */

DMotifBranchIndelModel *dmih_new_branch(double alpha, double beta, double tau,
					double epsilon, double t, int motif);
DMotifIndelModel *dmih_new_all(double alpha, double beta, double tau, 
			       double epsilon, TreeNode *tree);
DMotifIndelModel *dmih_new(double *alpha, double *beta, double *tau, 
			   double *epsilon, TreeNode *tree, int e, List *l,
			   int motif);
void dmih_free_branch(DMotifBranchIndelModel *bim);
void dmih_free(DMotifIndelModel *im);
void dmih_set_branch(DMotifBranchIndelModel *bim, double alpha, double beta, 
		     double tau, double epsilon, double t, int motif);
void dmih_set_all(DMotifIndelModel *im, double alpha, double beta, double tau, 
		  double epsilon, TreeNode *tree);
void dmih_set(DMotifIndelModel *im, double *alpha, double *beta, double *tau, 
	      double *epsilon, TreeNode *tree, int e, List *l, int motif);
double dmih_branch_column_logl(IndelHistory *ih, DMotifBranchIndelModel *bim, 
			       int child_id, double *col_logl);
double dmih_column_logl(IndelHistory *ih, DMotifIndelModel *im, 
			double *col_logl);
double dmih_likelihood(DMotifIndelModel *im, DMotifIndelSuffStats *iss);
DMotifBranchIndelSuffStats *dmih_suff_stats_branch(IndelHistory *ih, int child_id);
DMotifIndelSuffStats *dmih_suff_stats(IndelHistory *ih);
void dmih_free_suff_stats(DMotifIndelSuffStats *iss);
double dmih_simulate_history(DMotifIndelModel *tim, int ncols);
void dmih_estimate(DMotifIndelModel *im, IndelHistory *ih, 
		   DMotifIndelSuffStats *ss, FILE *logf);
DMotifBranchIndelSuffStats *dmih_suff_stats_branch_cat(IndelHistory *ih, 
						 int child_id, int *categories,
						 int do_cat);
DMotifIndelSuffStats *dmih_suff_stats_cat(IndelHistory *ih, int *categories, 
					  int do_cat);
MSA *dmih_as_alignment(IndelHistory *ih, MSA *msa, int debug);

#endif
