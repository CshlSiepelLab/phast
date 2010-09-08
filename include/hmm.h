/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: hmm.h,v 1.8 2009-02-03 22:11:28 agd27 Exp $ */

/* Library of functions relating to hidden Markov models.  Includes
 * simple reading and writing routines, as well as implementations of
 * the Viterbi algorithm, the forward algorithm, and the backward
 * algorithm.  Also includes a function to compute posterior
 * probabilities. */

#ifndef HMM_H
#define HMM_H

#include <matrix.h>
#include <markov_matrix.h>
#include <misc.h>
#include <lists.h>
#include <vector.h>
#include <prob_vector.h>

#define MAXSTATES 1000
#define BEGIN_STATE -99
#define END_STATE -98

#define BEGIN_TRANSITIONS_TAG "BEGIN_TRANSITIONS:"
#define END_TRANSITIONS_TAG "END_TRANSITIONS:"
#define TRANSITION_MATRIX_TAG "TRANSITION_MATRIX:"
#define EQ_FREQS_TAG "EQUILIBRIUM_FREQUENCIES:"

typedef enum {VITERBI, FORWARD, BACKWARD} hmm_mode;

/* NOTE: eventually need to be able to support an "adjacency list"
 * rather than an "adjacency matrix" representation of an HMM, for
 * better efficiency when there are many states and they are not fully
 * connected  */
typedef struct {
  int nstates;
  MarkovMatrix *transition_matrix;
  Matrix *transition_score_matrix; /* entries are logs of entries
				    * in transition matrix */
  Vector *begin_transitions, *end_transitions, 
    *begin_transition_scores, *end_transition_scores, *eq_freqs;
  List **predecessors, **successors;
  List *begin_successors, *end_predecessors;
} HMM;



HMM* hmm_new(MarkovMatrix *mm, Vector *eq_freqs,
             Vector *begin_transitions, 
             Vector *end_transitions);
HMM *hmm_new_nstates(int nstates, int begin, int end);
HMM *hmm_create_copy(HMM *src);
void hmm_free(HMM *hmm);
double hmm_get_transition_score(HMM *hmm, int from_state, int to_state);
HMM* hmm_new_from_file(FILE *F);
void hmm_print(FILE *F, HMM *hmm);
void hmm_viterbi(HMM *hmm, double **emission_scores, int seqlen, int *path);
double hmm_forward(HMM *hmm, double **emission_scores, int seqlen, 
                   double **forward_scores);
double hmm_backward(HMM *hmm, double **emission_scores, int seqlen,
                    double **backward_scores);
double hmm_posterior_probs(HMM *hmm, double **emission_scores, int seqlen,
                           double **posterior_probs);
void hmm_do_dp_forward(HMM *hmm, double **emission_scores, int seqlen, 
                       hmm_mode mode, double **full_scores, int **backptr);
void hmm_do_dp_backward(HMM *hmm, double **emission_scores, int seqlen, 
                        double **full_scores);
double hmm_max_or_sum(HMM *hmm, double **full_scores, double **emission_scores,
                      int **backptr, int i, int j, hmm_mode mode);

void hmm_dump_matrices(HMM *hmm, double **emission_scores, int seqlen,
                       double **full_scores, int **backptr);

void hmm_train_from_counts(HMM *hmm, Matrix *trans_counts, 
                           Matrix *trans_pseudocounts,
                           Vector *state_counts,
                           Vector *state_pseudocounts,
                           Vector *beg_counts, 
                           Vector *beg_pseudocounts);
void hmm_train_from_paths(HMM *hmm, int **path, int npaths,
                          Matrix *trans_pseudocounts, 
                          Vector *state_pseudocounts,int use_begin,
                          Vector *beg_pseudocounts);
void hmm_train_update_counts(Matrix *trans_counts, Vector *state_counts, 
                             Vector *beg_counts,
                             int *path, int len, int nstates);
HMM *hmm_create_trivial();
double hmm_path_likelihood(HMM *hmm, double **emission_scores, int seqlen, 
                           int *path);
double hmm_score_subset(HMM *hmm, double **emission_scores, List *states,
                        int begidx, int len);
double hmm_log_odds_subset(HMM *hmm, double **emission_scores, 
                           List *test_states, List *null_states,
                           int begidx, int len);
void hmm_cross_product(HMM *dest, HMM *src1, HMM *src2);
void hmm_reset(HMM *hmm);
HMM *hmm_reverse_compl(HMM *hmm, List *pivot_states, int *mapping);
void hmm_renormalize(HMM *hmm);
void hmm_stochastic_traceback(HMM *hmm, double **forward_scores, 
			      int seqlen, int *path);
void hmm_set_transition_score_matrix(HMM *hmm);

#endif
