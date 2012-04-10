/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: hmm.c,v 1.16 2009-03-09 16:33:04 agd27 Exp $ */

#include "hmm.h"
#include <math.h>
#include <misc.h>
#include "queues.h"
#include "stacks.h"
#include <vector.h>
#include <prob_vector.h>
#include <time.h>

/* Library of functions for manipulation of hidden Markov models.
   Includes simple reading and writing routines, as well as
   implementations of the Viterbi algorithm, the forward algorithm,
   and the backward algorithm.  Also includes a function to compute
   posterior probabilities. */

/* to do: handling of end states needs some work (not currently used
   much).  E.g., if there is an end state, the transition matrix will
   not be a true Markov matrix */


/* Creates a new HMM object based on a Markov matrix of transition
   probabilities, a vector of transitions from the begin state, and a
   vector of transitions to the end state.  All transitions are
   expected to be described as probabilities, not logs.  The
   vector of end-state transitions may be NULL, in which case it will
   be assumed that the end state is not of interest (e.g., the Viterbi
   and forward/backward algorithms will ignore the end state).  NULL
   may also be specified for the begin_transitions; in this case, they
   will be assumed to be uniform.   */
HMM* hmm_new(MarkovMatrix *mm, Vector *eq_freqs,
             Vector *begin_transitions, Vector *end_transitions) {
  HMM *hmm = (HMM*)smalloc(sizeof(HMM));
  int i;

/*   mm_validate(mm); */
  hmm->transition_matrix = mm;
  hmm->eq_freqs = eq_freqs;
  hmm->begin_transitions = begin_transitions;
  hmm->end_transitions = end_transitions;
  hmm->nstates = mm->size;
  hmm->transition_score_matrix = NULL;
  hmm->begin_transition_scores = hmm->end_transition_scores = NULL;
  hmm->predecessors = hmm->successors = NULL;
  hmm->begin_successors = hmm->end_predecessors = NULL;

  /* if begin_transitions are NULL, make them uniform */
  if (begin_transitions == NULL) {
    hmm->begin_transitions = vec_new(mm->size);
    for (i = 0; i < mm->size; i++) 
      vec_set(hmm->begin_transitions, i, 1.0/mm->size);
  }

  hmm_reset(hmm);

  return hmm;
}

/* Create a new HMM object with the specified number of states.
   Transition probabilities and equilibrium freqs will all be
   initialized to zero.  If begin == 1 (end == 1), the HMM will
   include a begin_transition (end_transition) vector, also
   initialized to zero; otherwise this vector will be NULL (see
   hmm_new) */
HMM *hmm_new_nstates(int nstates, int begin, int end) {
  Vector *eqfreqs = vec_new(nstates), *begv = NULL, *endv = NULL;
  vec_zero(eqfreqs);
  if (begin) {
    begv = vec_new(nstates);
    vec_zero(begv);
  }
  if (end) {
    endv = vec_new(nstates);
    vec_zero(endv);
  }
  return hmm_new(mm_new(nstates, NULL, DISCRETE), eqfreqs, begv, endv);
}

/* Create a copy of an HMM */
HMM *hmm_create_copy(HMM *src) {
  MarkovMatrix *transition_matrix = NULL;
  Vector *eq_freqs = NULL, *begin_transitions = NULL, 
    *end_transitions = NULL;

  if (src->transition_matrix != NULL) 
    transition_matrix = mm_create_copy(src->transition_matrix);

  if (src->eq_freqs != NULL) {
    eq_freqs = vec_new(src->nstates);
    vec_copy(eq_freqs, src->eq_freqs);
  }

  if (src->begin_transitions != NULL) {
    begin_transitions = vec_new(src->nstates);
    vec_copy(begin_transitions, src->begin_transitions);
  }

  if (src->end_transitions != NULL) {
    end_transitions = vec_new(src->nstates);
    vec_copy(end_transitions, src->end_transitions);
  }

  return hmm_new(transition_matrix, eq_freqs, begin_transitions, 
                 end_transitions);
}

/* Frees all memory associated with an HMM object */
void hmm_free(HMM *hmm) {
  int i;
  if (hmm->transition_matrix != NULL) 
    mm_free(hmm->transition_matrix);
  if (hmm->transition_score_matrix != NULL) 
    mat_free(hmm->transition_score_matrix);
  if (hmm->eq_freqs != NULL)
    vec_free(hmm->eq_freqs);
  if (hmm->begin_transitions != NULL)
    vec_free(hmm->begin_transitions);
  if (hmm->begin_transition_scores != NULL)
    vec_free(hmm->begin_transition_scores);
  if (hmm->end_transitions != NULL)
    vec_free(hmm->end_transitions);
  if (hmm->end_transition_scores != NULL)
    vec_free(hmm->end_transition_scores);
  for (i = 0; i < hmm->nstates; i++) {
    lst_free(hmm->predecessors[i]);
    lst_free(hmm->successors[i]);
  }
  lst_free(hmm->begin_successors);
  lst_free(hmm->end_predecessors);
  sfree(hmm->predecessors);
  sfree(hmm->successors);
  sfree(hmm);
}

/* Creates a new HMM object from the contents of a file.  File format
   is very simple; it includes specification of a transition matrix, a
   vector of transitions from the begin state, and (optionally) a
   vector of transitions to the end state.  The transition matrix
   should be specified using the conventions described in
   markov_matrix.c; the others should simply be lists of
   whitespace-separated numbers.  The description of each object must
   be preceded by a tag, as defined in hmm.h. */
HMM* hmm_new_from_file(FILE *F) {
  char tag[100];
  MarkovMatrix *mm=NULL;
  Vector *beg = NULL, *end = NULL, *eqfreqs = NULL;

  while (fscanf(F, "%s ", tag) != EOF) {
    if (!strcmp(tag, TRANSITION_MATRIX_TAG)) 
      mm = mm_new_from_file(F, DISCRETE);
    else if (!strcmp(tag, EQ_FREQS_TAG)) {
      if (mm == NULL)
        die("ERROR: transition matrix must come first in HMM file.\n");
      eqfreqs = vec_new(mm->size);
      vec_read(eqfreqs, F);
    }
    else if (!strcmp(tag, BEGIN_TRANSITIONS_TAG)) {
      if (mm == NULL)
        die("ERROR: transition matrix must come first in HMM file.\n");
      beg = vec_new(mm->size);
      vec_read(beg, F);
    }
    else if (!strcmp(tag, END_TRANSITIONS_TAG)) {
      if (mm == NULL)
        die("ERROR: transition matrix must come first in HMM file.\n");
      end = vec_new(mm->size);
      vec_read(end, F);
    }
  }

  return hmm_new(mm, eqfreqs, beg, end);
}

/* Prints textual representation of HMM object to file. */
void hmm_print(FILE *F, HMM *hmm) {
  fprintf(F, "%s\n", TRANSITION_MATRIX_TAG);
  mm_pretty_print(F, hmm->transition_matrix);
  if (hmm->eq_freqs != NULL) {
    fprintf(F, "%s\n", EQ_FREQS_TAG);
    vec_fprintf(hmm->eq_freqs, F, "%e ");
  }
  if (hmm->begin_transitions != NULL) {
    fprintf(F, "%s\n", BEGIN_TRANSITIONS_TAG);
    vec_fprintf(hmm->begin_transitions, F, "%e ");
  }
  if (hmm->end_transitions != NULL) {
    fprintf(F, "%s\n", END_TRANSITIONS_TAG);
    vec_fprintf(hmm->end_transitions, F, "%e ");
  }
}

/* Returns the log probability associated with the transition between
   any two states, including the special states BEGIN_STATE and
   END_STATE */
double hmm_get_transition_score(HMM *hmm, int from_state, int to_state) {
  int i, j;
  double prob;

  if (hmm->transition_matrix == NULL)
    die("ERROR hmm_get_transition_score: hmm->transition_matrix is NULL\n");
  if (from_state == BEGIN_STATE && to_state == END_STATE)
    die("ERROR: hmm_get_transition_score: from_state==BEGIN_STATE and to_state==END_STATE\n");
  
  if (from_state == BEGIN_STATE) {
    if (hmm->begin_transition_scores == NULL) {
				/* create "begin" score vector; use
				   logs */
      hmm->begin_transition_scores =
        vec_new(hmm->begin_transitions->size);
      for (i = 0; i < hmm->begin_transitions->size; i++) {
        prob = vec_get(hmm->begin_transitions, i);
        if (prob == 0)
          vec_set(hmm->begin_transition_scores, i, NEGINFTY);
        else
          vec_set(hmm->begin_transition_scores, i, log2(prob));
      }
    }
    return vec_get(hmm->begin_transition_scores, to_state);
  }
  else if (to_state == END_STATE) {
    if (hmm->end_transitions == NULL)
      return 0;
    if (hmm->end_transition_scores == NULL) {
      /* create "end" score vector; use
         logs */
      hmm->end_transition_scores =
        vec_new(hmm->end_transitions->size);
      for (i = 0; i < hmm->end_transitions->size; i++) {
        prob = vec_get(hmm->end_transitions, i);
        if (prob == 0)
          vec_set(hmm->end_transition_scores, i, NEGINFTY);
        else
          vec_set(hmm->end_transition_scores, i, log2(prob));
      }
    }
    return vec_get(hmm->end_transition_scores, from_state);
  }

  /* anything else must be a normal transition */
  if (hmm->transition_score_matrix == NULL) {
				/* create transition score matrix; use
				   logs */
    Matrix *m = mat_new(hmm->nstates, hmm->nstates);
    for (i = 0; i < hmm->nstates; i++) {
      for (j = 0; j < hmm->nstates; j++) {
        prob = mm_get(hmm->transition_matrix, i, j);
        if (prob == 0)
          mat_set(m, i, j, NEGINFTY);
        else
          mat_set(m, i, j, log2(prob));
      }
    }
    hmm->transition_score_matrix = m;
  }

  return (mat_get(hmm->transition_score_matrix, from_state, to_state));
}

/* Finds most probable path, according to the Viterbi algorithm.
   Emission scores must be passed in as a two dimensional matrix, with
   hmm->nstates rows and seqlen columns.  The array "path" must be
   allocated externally and be of length seqlen.  This array will be
   filled with integers indicating state numbers in the HMM. */
void hmm_viterbi(HMM *hmm, double **emission_scores, int seqlen, int *path) {

  double **full_scores;
  int **backptr;
  int i, j, len, bestidx;
  double besttran;

  /* set up necessary arrays */
  full_scores = (double**)smalloc(hmm->nstates * sizeof(double*));
  backptr = (int**)smalloc(hmm->nstates * sizeof(int*));
  len = seqlen;
  for (i = 0; i < hmm->nstates; i++) {
    full_scores[i] = (double*)smalloc(len * sizeof(double));
    backptr[i] = (int*)smalloc(len * sizeof(int));
  }

  /* fill array using DP */
  hmm_do_dp_forward(hmm, emission_scores, seqlen, VITERBI, full_scores, 
                    backptr);

  /* find path by backtracing */

  /* first find starting place */
  bestidx = 0; 
  besttran = hmm_get_transition_score(hmm, 0, END_STATE);
  for (i = 1; i < hmm->nstates; i++) {
    double thistran = hmm_get_transition_score(hmm, i, END_STATE);
    if (full_scores[i][len-1] + thistran > 
        full_scores[bestidx][len-1] + besttran) 
      bestidx = i;
                                /* note: when hmm->end_transitions ==
                                   NULL, besttran and thistran will
                                   always be zero (see function
                                   hmm_get_transition_score) */
  }

  /* now backtrace recursively */
  i = bestidx;
  j = seqlen - 1;
  while (i != -1) {
    path[j] = i;
    i = backptr[i][j];
    j--;
  }

  for (i = 0; i < hmm->nstates; i++) {
    sfree(full_scores[i]);
    sfree(backptr[i]);
  }
  sfree(full_scores);
  sfree(backptr);
}

/* Fills matrix of "forward" scores and returns total log probability
   of sequence.  As above, emission scores must be passed in as a two
   dimensional matrix with hmm->nstates rows and seqlen columns.  Here
   the array forward_scores must be allocated externally as well, to
   the same size.  It will be filled by this function. */
double hmm_forward(HMM *hmm, double **emission_scores, int seqlen, 
                   double **forward_scores) {
  double llh;
/*   int t0, t1; */

/*   t0 = (int)time(0); */
  hmm_do_dp_forward(hmm, emission_scores, seqlen, FORWARD, forward_scores, 
                    NULL);
  llh = hmm_max_or_sum(hmm, forward_scores, NULL, NULL, END_STATE, 
                        seqlen, FORWARD);
/*   t1 = (int)time(0); */
/*   fprintf(stderr, "Forward algorithm time elapsed: %d seconds\n", (t1 - t0)); */
  return llh;
}

/* Fills matrix of "backward" scores and returns total log probability
   of sequence.  As above, emission scores must be passed in as a two
   dimensional matrix with hmm->nstates rows and seqlen columns.  Here
   the array backward_scores must be allocated externally as well, to
   the same size.  It will be filled by this function. */
double hmm_backward(HMM *hmm, double **emission_scores, int seqlen,
                    double **backward_scores) {

  hmm_do_dp_backward(hmm, emission_scores, seqlen, backward_scores);

  return hmm_max_or_sum(hmm, backward_scores, emission_scores, NULL, 
                        BEGIN_STATE, -1, BACKWARD);
}

/* Fills matrix of posterior probabilities.  As above, emission scores
   must be passed in as a two dimensional matrix with hmm->nstates
   rows and seqlen columns.  Here the array posterior_probs_scores
   must be allocated externally as well, to the same size.  It will be
   filled by this function.  This function calls hmm_forward and
   hmm_backward, but it transparently handles the management of the
   arrays used by those routines.  NOTE: if the posterior probs for
   any state i are not desired, set posterior_probs[i] = NULL.  The
   return value is the log likelihood.  */
double hmm_posterior_probs(HMM *hmm, double **emission_scores, int seqlen,
                         double **posterior_probs) {
  int i, j, len;
  double logp_fw, logp_bw;
  double **forward_scores, **backward_scores;
  List *val_list;

  len = seqlen;

  /* allocate arrays for forward and backward algs */
  forward_scores = (double**)smalloc(hmm->nstates * sizeof(double*));
  backward_scores = (double**)smalloc(hmm->nstates * sizeof(double*));
  for (i = 0; i < hmm->nstates; i++) {
    forward_scores[i] = (double*)smalloc(len * sizeof(double));
    backward_scores[i] = (double*)smalloc(len * sizeof(double));
  }

  /* run forward and backward algs */
  logp_fw = hmm_forward(hmm, emission_scores, seqlen, forward_scores); 
  logp_bw = hmm_backward(hmm, emission_scores, seqlen, backward_scores);

  if (fabs(logp_fw - logp_bw) > 1.0)
    fprintf(stderr, "WARNING: forward and backward algorithms returned different total log\nprobabilities (%f and %f, respectively).\n", logp_fw, logp_bw);

  /* compute posterior probs */
  val_list = lst_new_dbl(hmm->nstates);
  for (j = 0; j < len; j++) {
    double this_logp;
    checkInterruptN(i, 1000);

    /* to avoid rounding errors, estimate total log prob
       separately for each column */
    lst_clear(val_list);
    for (i = 0; i < hmm->nstates; i++) 
      lst_push_dbl(val_list, (forward_scores[i][j] + backward_scores[i][j]));
    this_logp = log_sum(val_list);

    for (i = 0; i < hmm->nstates; i++) 
      if (posterior_probs[i] != NULL) /* indicates probs for this
                                         state are not desired */
        posterior_probs[i][j] = exp2(forward_scores[i][j] + 
                                     backward_scores[i][j] - this_logp);
  }

  for (i = 0; i < hmm->nstates; i++) {
    sfree(forward_scores[i]);
    sfree(backward_scores[i]);
  }
  sfree(forward_scores);
  sfree(backward_scores);
  lst_free(val_list);

  return logp_fw;
}

/* This is the core dynamic programming routine used by hmm_viterbi
   and hmm_forward.  It is not intended to be called directly. */
void hmm_do_dp_forward(HMM *hmm, double **emission_scores, int seqlen, 
                       hmm_mode mode, double **full_scores, int **backptr) {  

  int i, j;

  if (!(seqlen > 0 && hmm != NULL && hmm->nstates > 0 && 
	(mode == VITERBI || mode == FORWARD) && 
	full_scores != NULL && (mode != VITERBI || backptr != NULL)))
    die("ERROR hmm_do_dp_forward: bad params\n");

  /* initialization */
  for (i = 0; i < hmm->nstates; i++) {
    full_scores[i][0] = emission_scores[i][0] +
      hmm_get_transition_score(hmm, BEGIN_STATE, i);
    if (mode == VITERBI) backptr[i][0] = -1;
  }

  /* recursion */
  for (j = 1; j < seqlen; j++) {
    for (i = 0; i < hmm->nstates; i++) {
      full_scores[i][j] = emission_scores[i][j] + 
        hmm_max_or_sum(hmm, full_scores, emission_scores, backptr, 
                       i, j, mode);
    }
  }

#ifdef DEBUG
  hmm_dump_matrices(hmm, emission_scores, seqlen, full_scores, backptr);
#endif
}

/* This is the core dynamic programming routine used by hmm_backward.
   It is not intended to be called directly. */
void hmm_do_dp_backward(HMM *hmm, double **emission_scores,  int seqlen, 
                        double **full_scores) {  

  int i, j;

  if (!(seqlen > 0 && hmm != NULL && hmm->nstates > 0 && 
	full_scores != NULL))
    die("ERROR hmm_do_dp_backward: bad params\n");

  /* initialization */
  for (i = 0; i < hmm->nstates; i++)
    full_scores[i][seqlen-1] = hmm_get_transition_score(hmm, i, END_STATE);
                                /*  will be 0 when no end state */

  /* recursion */
  for (j = seqlen - 2; j >= 0; j--) {
    checkInterruptN(j, 1000);
    for (i = 0; i < hmm->nstates; i++) {
      full_scores[i][j] = 
        hmm_max_or_sum(hmm, full_scores, emission_scores, NULL, 
                       i, j, BACKWARD);
    }
  }
}

/* Finds max or sum of score/transition combination over all previous
   states (max for Viterbi, sum for forward/backward).  In Viterbi
   case, sets backpointer as a side-effect.  NOTE: 'i' is the present
   state and 'j' the present column.  The state 'i' may be END_STATE
   (only if mode == FORWARD) or BEGIN_STATE (only if mode ==
   BACKWARD).  */  
double hmm_max_or_sum(HMM *hmm, double **full_scores, double **emission_scores,
                      int **backptr, int i, int j, hmm_mode mode) { 
  int k;
  double retval = NEGINFTY;
  
  static List *l = NULL;

  if (l == NULL) {
    l = lst_new_dbl(hmm->nstates);
    set_static_var((void**)&l);
  }

  if (mode == VITERBI) {
    for (k = 0; k < lst_size(hmm->predecessors[i]); k++) {
      int pred;
      double candidate;
      pred = lst_get_int(hmm->predecessors[i], k);
      if (pred == BEGIN_STATE) continue;
      candidate = full_scores[pred][j-1] + 
        hmm_get_transition_score(hmm, pred, i);

      if (candidate > retval) {
        retval = candidate;
        backptr[i][j] = pred;
      }
    }
  }
  else if (mode == FORWARD) {
    List *pred_lst = (i == END_STATE ? hmm->end_predecessors :
      hmm->predecessors[i]);
    lst_clear(l);
    for (k = 0; k < lst_size(pred_lst); k++) {
      int pred;
      double candidate;
      pred = lst_get_int(pred_lst, k);
      if (pred == BEGIN_STATE) continue;
      candidate = full_scores[pred][j-1] +  
        hmm_get_transition_score(hmm, pred, i);
      lst_push_dbl(l, candidate);
    }
  }
  else {                        /* mode == BACKWARD */
    List *succ_lst = (i == BEGIN_STATE ? hmm->begin_successors :
      hmm->successors[i]);
    lst_clear(l);
    for (k = 0; k < lst_size(succ_lst); k++) {
      int succ;
      double candidate;
      succ = lst_get_int(succ_lst, k);
      if (succ == END_STATE) continue;
      candidate = emission_scores[succ][j+1] + full_scores[succ][j+1]
        + hmm_get_transition_score(hmm, i, succ);
      lst_push_dbl(l, candidate);
    }
  }    

  if (mode == FORWARD || mode == BACKWARD)
    retval = log_sum(l);

  return retval;
}

/* reset arcs of HMM according to a matrix of counts and, optionally,
   a specified matrix of pseudocounts.  If beg_counts is NULL, uniform
   beg transitions will be used.  Either pseudocounts parameter may be
   NULL.  TODO: allow end counts to be specified  */
void hmm_train_from_counts(HMM *hmm, Matrix *trans_counts, 
                           Matrix *trans_pseudocounts, 
                           Vector *state_counts,
                           Vector *state_pseudocounts,
                           Vector *beg_counts, 
                           Vector *beg_pseudocounts) {

  Matrix *countmat;
  Vector *statecount, *begcount;
  MarkovMatrix *mm;
  int i;
  double sum;

  /* add pseudocounts (if necessary) */
  if (trans_pseudocounts == NULL) 
    countmat = trans_counts;
  else {
    countmat = mat_new(hmm->nstates, hmm->nstates);
    mat_copy(countmat, trans_counts);
    mat_plus_eq(countmat, trans_pseudocounts);
  }
  if (state_pseudocounts == NULL)
    statecount = state_counts;
  else {
    statecount = vec_new(hmm->nstates);
    vec_copy(statecount, state_counts);
    vec_plus_eq(statecount, state_pseudocounts);
  }     
  if (!(beg_pseudocounts == NULL || beg_counts != NULL))
    die("ERROR hmm_train_from_counts: beg_pseudocounts==NULL=%i, beg_counts==NULL=%i\n",
	beg_pseudocounts==NULL, beg_counts==NULL);
  if (beg_pseudocounts == NULL)
    begcount = beg_counts;      /* could be NULL */
  else {
    begcount = vec_new(hmm->nstates);
    vec_copy(begcount, beg_counts);
    vec_plus_eq(begcount, beg_pseudocounts);
  }

  /* create a markov matrix */
  mm = mm_new_from_counts(countmat, hmm->transition_matrix->states); 
  
  /* free old data */
  if (hmm->transition_matrix) {
    mm_free(hmm->transition_matrix);
    hmm->transition_matrix = NULL; 
  }
  if (hmm->transition_score_matrix) {
    mat_free(hmm->transition_score_matrix);
    hmm->transition_score_matrix = NULL; 
  }
  if (hmm->eq_freqs) {
    vec_free(hmm->eq_freqs);
    hmm->eq_freqs = NULL; 
  }
  if (hmm->begin_transition_scores) {
    vec_free(hmm->begin_transition_scores);
    hmm->begin_transition_scores = NULL; 
  }
  if (hmm->end_transitions) {
    vec_free(hmm->end_transitions);
    hmm->end_transitions = NULL; 
  }
  if (hmm->end_transition_scores) {
    vec_free(hmm->end_transition_scores);
    hmm->end_transition_scores = NULL; 
  }

  /* set eq freqs */
  hmm->eq_freqs = vec_new(hmm->nstates);
  sum = 0;
  for (i = 0; i < hmm->nstates; i++)
    sum += vec_get(statecount, i);
  if (sum <= 0)
    die("ERROR hmm_train_from_counts sum=%f, should be >0\n", sum);
  vec_copy(hmm->eq_freqs, statecount);
  vec_scale(hmm->eq_freqs, 1/sum);

  /* set begin transitions */
  if (hmm->begin_transitions == NULL)
    hmm->begin_transitions = vec_new(hmm->nstates);
  if (begcount == NULL)         /* use uniform probs */
    for (i = 0; i < hmm->nstates; i++) 
      vec_set(hmm->begin_transitions, i, (double)1/hmm->nstates);
  else {                        /* use counts */
    sum = 0;
    for (i = 0; i < hmm->nstates; i++) sum += vec_get(begcount, i);
    for (i = 0; i < hmm->nstates; i++) 
      vec_set(hmm->begin_transitions, i, 
                     safediv(vec_get(begcount, i), sum));
  }

  hmm->transition_matrix = mm;

  if (trans_pseudocounts != NULL) mat_free(countmat);
  if (state_pseudocounts != NULL) vec_free(statecount);
  if (beg_pseudocounts != NULL) vec_free(begcount);
}

/* retrain HMM according to observed path and, optionally, specified
   matrix of pseudocounts.  Path indices are expected to range between 1
   and hmm->nstates.  Each path is expected to be terminated by -1. */
void hmm_train_from_paths(HMM *hmm, int **path, int npaths,
                          Matrix *trans_pseudocounts, 
                          Vector *state_pseudocounts,
                          int use_begin, Vector *beg_pseudocounts) {
  Matrix *trans_counts = mat_new(hmm->nstates, hmm->nstates);
  Vector *state_counts = vec_new(hmm->nstates);
  Vector *begcounts = use_begin ? vec_new(hmm->nstates) : NULL;
  int i, j;
  mat_zero(trans_counts);
  vec_zero(state_counts);
  if (begcounts != NULL) vec_zero(begcounts);
  for (i = 0; i < npaths; i++) {
    for (j = 0; path[i][j] != -1; j++) {
      if (!(path[i][j] >= 0 && path[i][j] < hmm->nstates))
	die("ERROR hmm_train_from_paths: path[%i][%i]=%i, should be in [0,%i)\n",
	    i, j, path[i][j], hmm->nstates);
      vec_set(state_counts, path[i][j], 
                     vec_get(state_counts, path[i][j]) + 1);
      if (j > 0)
        mat_set(trans_counts, path[i][j-1], path[i][j], 
                       mat_get(trans_counts, path[i][j-1], 
                                      path[i][j]) + 1);
    }
    if (use_begin && path[i][0] != -1) 
      vec_set(begcounts, path[i][0], 
                     vec_get(begcounts, path[i][0]) + 1);
  }
  hmm_train_from_counts(hmm, trans_counts, trans_pseudocounts, state_counts, 
                        state_pseudocounts, begcounts, beg_pseudocounts);
  mat_free(trans_counts);
  vec_free(state_counts);
  if (use_begin) vec_free(begcounts);
}

/* update counts according to new path */
void hmm_train_update_counts(Matrix *trans_counts, Vector *state_counts, 
                             Vector *beg_counts, int *path, int len, 
                             int nstates) {
  int j;
  for (j = 0; j < len; j++) {
    if (!(path[j] >= 0 && path[j] < nstates))
      die("ERROR hmm_train_update_counts: path[%i]=%i, should be in [0, %i)\n",
	  j, path[j], nstates);
    vec_set(state_counts, path[j], 
                   vec_get(state_counts, path[j]) + 1);
    if (j > 0)
      mat_set(trans_counts, path[j-1], path[j], 
                     mat_get(trans_counts, path[j-1], path[j]) + 1);
  }
  if (beg_counts != NULL && len > 0)
    vec_set(beg_counts, path[0], 
                   vec_get(beg_counts, path[0]) + 1);
  /* temporary */
  for (j = 0; j < state_counts->size; j++) 
    if (vec_get(state_counts, j) < 0)
      die("ERROR hmm_train_update_counts: state_counts[%j]=%f\n",
	  j, vec_get(state_counts, j));
}

/* for debugging */
void hmm_dump_matrices(HMM *hmm, double **emission_scores, int seqlen,
                       double **full_scores, int **backptr) {
  FILE *F = phast_fopen("hmm.debug", "w+");
  int i, j;
  char tmpstr[50];

  fprintf(F, "EMISSION SCORES:\n");
  fprintf(F, "    ");
  for (j = 0; j < seqlen; j++) 
    fprintf(F, "%10d ", j);
  fprintf(F, "\n");
  for (i = 0; i < hmm->nstates; i++) { 
    fprintf(F, "%2d: ", i);
    for (j = 0; j < seqlen; j++) {
      if (emission_scores[i][j] <= NEGINFTY)
        strcpy(tmpstr, "-INF");
      else 
        sprintf(tmpstr, "%.3f", emission_scores[i][j]);
      fprintf(F, "%10s ", tmpstr);
    }
    fprintf(F, "\n");
  }

  fprintf(F, "\nFULL SCORES:\n");
  fprintf(F, "    ");
  for (j = 0; j < seqlen; j++) 
    fprintf(F, "%10d ", j);
  fprintf(F, "\n");
  for (i = 0; i < hmm->nstates; i++) { 
    fprintf(F, "%2d: ", i);
    for (j = 0; j < seqlen; j++) {
      if (full_scores[i][j] <= NEGINFTY)
        strcpy(tmpstr, "-INF");
      else 
        sprintf(tmpstr, "%.3f", full_scores[i][j]);
      fprintf(F, "%10s ", tmpstr);
    }
    fprintf(F, "\n");
  }

  if (backptr != NULL) {
    fprintf(F, "\nBACK POINTERS:\n");
    fprintf(F, "    ");
    for (j = 0; j < seqlen; j++) 
      fprintf(F, "%10d ", j);
    fprintf(F, "\n");
    for (i = 0; i < hmm->nstates; i++) { 
      fprintf(F, "%2d: ", i);
      for (j = 0; j < seqlen; j++) {
        fprintf(F, "%10d ", backptr[i][j]);
      }
      fprintf(F, "\n");
    }
  }

  phast_fclose(F);
}

/* creates and returns a trivial HMM, with a single state that
   transitions to itself with probability 1; useful for causing
   general HMM-based models to collapse to simpler models  */
HMM *hmm_create_trivial() {
  MarkovMatrix *triv_mm = mm_new(1, "1", DISCRETE);
  Vector *triv_begin_trans = vec_new(1);
  Vector *triv_eqfreq = vec_new(1);
  mm_set(triv_mm, 0, 0, 1);
  vec_set(triv_begin_trans, 0, 1);
  vec_set(triv_eqfreq, 0, 1);
  return hmm_new(triv_mm, triv_eqfreq, triv_begin_trans, NULL);  
}

/* define an HMM as the cross product of two other HMMs.  The pair of
   state i from src1 and state j from src2 is represented in the new
   HMM by state i*src2->nstates + j */
void hmm_cross_product(HMM *dest, HMM *src1, HMM *src2) {
  int i, j, k, l;

  for (i = 0; i < src1->nstates; i++) {
    checkInterrupt();
    for (j = 0; j < src2->nstates; j++) 
      for (k = 0; k < src1->nstates; k++) 
        for (l = 0; l < src2->nstates; l++) 
          mm_set(dest->transition_matrix, i*src2->nstates + j, 
                 k*src2->nstates + l, mm_get(src1->transition_matrix, i, k) * 
                 mm_get(src2->transition_matrix, j, l));
  }

  if (src1->eq_freqs != NULL && src2->eq_freqs != NULL) {
    for (i = 0; i < src1->nstates; i++) 
      for (j = 0; j < src2->nstates; j++) 
        vec_set(dest->eq_freqs, i*src2->nstates + j, 
                       vec_get(src1->eq_freqs, i) * 
                       vec_get(src2->eq_freqs, j));
  }
  else dest->eq_freqs = NULL;

  for (i = 0; i < src1->nstates; i++) 
    for (j = 0; j < src2->nstates; j++) 
      vec_set(dest->begin_transitions, i*src2->nstates + j, 
                     vec_get(src1->begin_transitions, i) * 
                     vec_get(src2->begin_transitions, j));

  if (src1->end_transitions != NULL && src2->end_transitions != NULL) {
    for (i = 0; i < src1->nstates; i++) 
      for (j = 0; j < src2->nstates; j++) 
        vec_set(dest->end_transitions, i*src2->nstates + j, 
                       vec_get(src1->end_transitions, i) * 
                       vec_get(src2->end_transitions, j));
  }
  else dest->end_transitions = NULL;

  hmm_reset(dest);
}


/* compute the total log likelihood of a specified path */
double hmm_path_likelihood(HMM *hmm, double **emission_scores, int seqlen, 
                           int *path) {
  int i;
  double l = 0;
  if (seqlen <= 0) return 0;
  l = hmm_get_transition_score(hmm, BEGIN_STATE, path[0]) +
    emission_scores[path[0]][0];
  for (i = 1; i < seqlen; i++)
    l += hmm_get_transition_score(hmm, path[i-1], path[i]) +
      emission_scores[path[i]][i];
  l += hmm_get_transition_score(hmm, path[seqlen-1], END_STATE);
  return l;
}

/* Compute the total log likelihood of a subsequence of the input,
   using only the specified states.  Useful for scoring candidate
   predictions.  The parameter "states" must be a list of indices of
   states.  This routine is not as efficient as it could be (trying to
   reuse code).  */
double hmm_score_subset(HMM *hmm, double **emission_scores, List *states,
                        int begidx, int len) {
  double **forward_scores;
  double **dummy_emissions;
  int do_state[hmm->nstates];
  int i, j;
  double retval;
  Vector *orig_begin;
  MarkovMatrix *orig_trans;

  forward_scores = smalloc(hmm->nstates * sizeof(double*));
  dummy_emissions = smalloc(hmm->nstates * sizeof(double*));
  for (i = 0; i < hmm->nstates; i++) 
    forward_scores[i] = smalloc(len * sizeof(double));

  for (i = 0; i < hmm->nstates; i++) do_state[i] = 0;
  for (i = 0; i < lst_size(states); i++) do_state[lst_get_int(states, i)] = 1;

  /* set up a dummy emissions array */
  for (i = 0; i < hmm->nstates; i++) 
    dummy_emissions[i] = &(emission_scores[i][begidx]);

  /* need to tweak the begin transitions to be sure that the HMM can
     make it into the states in question.  We'll simply use a uniform
     distribution over the allowable states */
  orig_begin = hmm->begin_transitions;
  hmm->begin_transitions = vec_new(hmm->nstates);
  for (i = 0; i < hmm->nstates; i++) 
    vec_set(hmm->begin_transitions, i, 
                   do_state[i] ? 1.0/lst_size(states) : 0);

  /* also need to adjust the transition probs and renormalize them */
  orig_trans = hmm->transition_matrix;
  hmm->transition_matrix = mm_create_copy(hmm->transition_matrix);
  for (i = 0; i < hmm->nstates; i++) 
    for (j = 0; j < hmm->nstates; j++) 
      if (!do_state[i] || !do_state[j]) 
        mm_set(hmm->transition_matrix, i, j, 0);
                                /* hmm_renormalize will clean up rows
                                   of all zeros */
  hmm_renormalize(hmm);

  /* FIXME: now that I'm changing the transition probs anyway, maybe
     should just drop the extra states altogether (more efficient);
     wouldn't actually be that much more complicated */
  
  retval = hmm_forward(hmm, dummy_emissions, len, forward_scores);

  vec_free(hmm->begin_transitions);
  hmm->begin_transitions = orig_begin;
  mm_free(hmm->transition_matrix);
  hmm->transition_matrix = orig_trans;
  hmm_reset(hmm);

  for (i = 0; i < hmm->nstates; i++) 
    sfree(forward_scores[i]);
  sfree(forward_scores);
  sfree(dummy_emissions);

  return retval;
}

/* report a log odds score for a subsequence of the input, comparing
   the likelihood based on one subset of states (test_states) against
   the likelihood based on another subset (null_states).  Useful for
   scoring candidate predictions. */
double hmm_log_odds_subset(HMM *hmm, double **emission_scores, 
                                   List *test_states, List *null_states,
                                   int begidx, int len) {
  return (hmm_score_subset(hmm, emission_scores, test_states, begidx, len) -
          hmm_score_subset(hmm, emission_scores, null_states, begidx, len));
}


/* Reset various attributes that are derived from the underlying
   matrix of transitions.  Should be called after the matrix is
   changed for any reason.  Note: this routine assumes that
   hmm->nstates does not change. */
void hmm_reset(HMM *hmm) {
  int i, j;

  /* set up predecessor and successor lists */
  /* Note: lists for ordinary states include begin and end states;
     separate lists provide successors of begin and
     predecessors of end */
  if (hmm->predecessors == NULL) {
    hmm->predecessors = (List**)smalloc(hmm->nstates * sizeof(List*));
    for (i = 0; i < hmm->nstates; i++) 
      hmm->predecessors[i] = lst_new_int(hmm->nstates);
  }
  else {
    for (i = 0; i < hmm->nstates; i++)
      lst_clear(hmm->predecessors[i]);
  }
  if (hmm->successors == NULL) {
    hmm->successors = (List**)smalloc(hmm->nstates * sizeof(List*));
    for (i = 0; i < hmm->nstates; i++) 
      hmm->successors[i] = lst_new_int(hmm->nstates);
  }
  else {
    for (i = 0; i < hmm->nstates; i++)
      lst_clear(hmm->successors[i]);
  }
  for (i = 0; i < hmm->nstates; i++) {
    for (j = 0; j < hmm->nstates; j++) {
      if (mm_get(hmm->transition_matrix, i, j) > 0) {
        lst_push_int(hmm->predecessors[j], i);
        lst_push_int(hmm->successors[i], j);
      }
    }
  }
  if (hmm->begin_successors == NULL) 
    hmm->begin_successors = lst_new_int(hmm->nstates);
  else 
    lst_clear(hmm->begin_successors);
  if (hmm->end_predecessors == NULL) 
    hmm->end_predecessors = lst_new_int(hmm->nstates);
  else
    lst_clear(hmm->end_predecessors);
  for (i = 0; i < hmm->nstates; i++) {
    if (hmm->begin_transitions != NULL &&
        vec_get(hmm->begin_transitions, i) > 0) {
      lst_push_int(hmm->begin_successors, i);
      lst_push_int(hmm->predecessors[i], BEGIN_STATE);
    }
    if (hmm->end_transitions == NULL ||
        vec_get(hmm->end_transitions, i) > 0) {
      lst_push_int(hmm->end_predecessors, i);
      lst_push_int(hmm->successors[i], END_STATE);
    }
    /* Note that I've adopted the convention of considering all states
       predecessors of the end state, when no end state is required;
       this simplifies coding somewhat ... */
  }

  /* below is inefficient on repeated calls, but shouldn't be used
     heavily */
  if (hmm->transition_score_matrix != NULL) {
    mat_free(hmm->transition_score_matrix);
    hmm->transition_score_matrix = NULL;
  }
  if (hmm->begin_transition_scores != NULL) {
    vec_free(hmm->begin_transition_scores);
    hmm->begin_transition_scores = NULL;
  }
  if (hmm->end_transition_scores != NULL) {
    vec_free(hmm->end_transition_scores);
    hmm->end_transition_scores = NULL;
  }
}

/* Given an HMM, some of whose states represent strand-specific
   phenomena in DNA (e.g., coding regions, UTRs, introns), create a
   new HMM with a second version of all such states, corresponding to
   the reverse strand (assuming the original HMM describes the forward
   strand).  Transition probabilities between reverse-strand states
   will be a "reflection" of those for the forward-strand states,
   based on an assumption of reversibility.  No transitions will be
   allowed between forward- and reverse-strand states except via
   designated "pivot" states, which are not reflected (e.g.,
   background or intergenic states).  State 0 is automatically a pivot
   state; additional pivot states may be defined by the list
   "pivot_states", which, if non-NULL, is expected to contain integers
   corresponding to state indices.  All non-pivot states will be
   reflected.  States on the forward and reverse strands are assumed
   to have equal overall probability.  On return, the array "mapping"
   (must be preallocated to a size of hmm->nstates*2 -
   lst_size(pivot_states) - 1), defines the correspondence between
   forward and reverse states: mapping[i] == i if i is a forward or
   pivot state, and mapping[i] = -j if i is a reflection of j.  */ 
HMM *hmm_reverse_compl(HMM *hmm, List *pivot_states, int *mapping) {
  int npiv_states = 1 + (pivot_states != NULL ? lst_size(pivot_states) : 0);
  List *effective_pivot_states = lst_new_int(npiv_states);
  int nstates = hmm->nstates*2 - npiv_states;
  int is_pivot[nstates];
  HMM *retval = hmm_new_nstates(nstates, 
                                hmm->begin_transitions != NULL,
                                hmm->end_transitions != NULL);
  double prob;
  int i, j, lidx;

  lst_push_int(effective_pivot_states, 0);

  is_pivot[0] = 1;
  for (i = 1; i < nstates; i++) is_pivot[i] = 0;
  if (pivot_states != NULL) {
    for (i = 0; i < lst_size(pivot_states); i++) {
      is_pivot[lst_get_int(pivot_states, i)] = 1;
      lst_push_int(effective_pivot_states, lst_get_int(pivot_states, i));
    }
  }

  for (i = 0, j = hmm->nstates; i < hmm->nstates; i++) {
    mapping[i] = i;
    if (!is_pivot[i]) mapping[j++] = -i;
  }
  if (j != nstates)
    die("ERROR hmm_reverse_compl: j (%i) != nstates (%i)\n", j, nstates);

  /* copy transitions from forward states unchanged, and halve
     transition probs from pivot states to non-pivot (forward)
     states  */
  for (i = 0; i < hmm->nstates; i++) {
    for (j = 0; j < hmm->nstates; j++) {
      prob = mm_get(hmm->transition_matrix, i, j);
      if (is_pivot[i] && !is_pivot[j])
        prob /= 2;
      mm_set(retval->transition_matrix, i, j, prob);
    }
  }

  /* derive probs between reverse states from probs between
     corresponding forward states, using the "reversibility"
     principle */
  for (i = hmm->nstates; i < nstates; i++) {
    for (j = hmm->nstates; j < nstates; j++) {
      if (vec_get(hmm->eq_freqs, -mapping[i]) == 0)
        prob = 0;               /* avoid div by zero */
      else
        prob = mm_get(hmm->transition_matrix, -mapping[j], -mapping[i]) *
          vec_get(hmm->eq_freqs, -mapping[j]) / 
          vec_get(hmm->eq_freqs, -mapping[i]);
      mm_set(retval->transition_matrix, i, j, prob);
    }
  }

  /* set probs between reverse states and pivot states, using the same
     principle.  */ 
  for (lidx = 0; lidx < lst_size(effective_pivot_states); lidx++) {
    i = lst_get_int(effective_pivot_states, lidx); /* pivot state */
    for (j = hmm->nstates; j < nstates; j++) { /* reverse state */
      /* from pivot state to reverse state */
      if (vec_get(hmm->eq_freqs, i) == 0)
        prob = 0;
      else 
        prob = mm_get(hmm->transition_matrix, -mapping[j], i) *
          vec_get(hmm->eq_freqs, -mapping[j]) / 
          (2 * vec_get(hmm->eq_freqs, i)); /* don't forget to halve */
      mm_set(retval->transition_matrix, i, j, prob);

      /* from reverse state to pivot state */
      if (vec_get(hmm->eq_freqs, -mapping[j]) == 0)
        prob = 0;
      else
        prob = mm_get(hmm->transition_matrix, i, -mapping[j]) *
          vec_get(hmm->eq_freqs, i) / 
          vec_get(hmm->eq_freqs, -mapping[j]);
      mm_set(retval->transition_matrix, j, i, prob);
    }    
  }

  /* (transition probs between forward and reverse states remain zero) */

  /* set new eq freqs */
  for (i = 0; i < nstates; i++) {
    prob = vec_get(hmm->eq_freqs, abs(mapping[i]));
    vec_set(retval->eq_freqs, i, 
                   is_pivot[i] ? prob : prob/2);
  }

  /* set begin and end vectors */
  if (hmm->begin_transitions != NULL) {
    for (i = 0; i < nstates; i++) {
      prob = vec_get(hmm->begin_transitions, abs(mapping[i]));
      vec_set(retval->begin_transitions, i, 
                     is_pivot[i] ? prob : prob/2);
    }
  }
  if (hmm->end_transitions != NULL) {
    for (i = 0; i < nstates; i++) {
      prob = vec_get(hmm->end_transitions, abs(mapping[i]));
      vec_set(retval->end_transitions, i, prob); 
                                /* NOTE: conditional prob */
    }
  }

  hmm_renormalize(retval);
  lst_free(effective_pivot_states);
  
  return retval;
}

/* Renormalize transition probabilities so that all rows sum to 1.
   Currently does not consider begin or end transitions.  Calls
   hmm_reset.  */
void hmm_renormalize(HMM *hmm) {
  mm_renormalize(hmm->transition_matrix);
  hmm_reset(hmm);
}

/* Sample a state path through a sequence using the stochastic traceback
   algorithm. */
void hmm_stochastic_traceback(HMM *hmm, double **forward_scores, int seqlen,
			      int *path) {
  int i, j, k, pass, maxidx, state;
  double max, score, z;
  List *predecessors;
  Vector *pv;
  
  /* Initialization */
  state = END_STATE;
  
  /* Recursion */
  for (i = seqlen; i > 0; i--) {
    checkInterruptN(i, 1000);
    max = -INFTY;
    maxidx=0;
    z = 1;
    predecessors = (state == END_STATE ? hmm->end_predecessors :
			  hmm->predecessors[state]);
    pv = vec_new(lst_size(predecessors));
    vec_zero(pv);
    /* To avoid underflows, normalization nust be done in log space before
       exponentiation of the probabililites. This requires three passes for
       each site. */
    for (pass = 0; pass < 3; pass++) { /* First pass just finds the max */
      for (j = 0; j < lst_size(predecessors); j++) {
	k = lst_get_int(predecessors, j);
	if (k == BEGIN_STATE)
	  continue;
	if (pass == 0) {
	  score = (forward_scores[k][i-1]
		   + hmm_get_transition_score(hmm, k, state));
	  if (score > max) {
	    max = score;
	    maxidx = j;
	  }
	} else if (pass == 1) { /* Second pass computes the summation portion
				   of the normalization factor */
	  if (j == maxidx)
	    continue;
	  z += exp2((forward_scores[k][i-1]
		    + hmm_get_transition_score(hmm, k, state)) - max);
	} else { /* Third pass finishes computation of the normalization factor
		    and performs the stochastic traceback recurrence */
	  if (j == 0) {
	    z = max + log2(z); /* Take log of summation (initilaized to 1, so
				  no need to add 1) and add l(max) to get
				  normalization factor. */
	  }
	  /* The core recurrence, with normalization in log space */
	  pv->data[j] = exp2(forward_scores[k][i-1]
		       + hmm_get_transition_score(hmm, k, state) - z);
	}
      }
    }
    /* Draw an index from the distribution and convert to a state */
    state = lst_get_int(predecessors, pv_draw_idx(pv));
    vec_free(pv); /* Cannot reuse this, as size may not be constant */
    path[i-1] = state;
  }
}

/* Set the transition_score_matrix in an hmm object. This must be done before
   calling any functions that use hmm_get_transition_score in a multithreaded
   context. */
void hmm_set_transition_score_matrix(HMM *hmm) {
  int i, j;
  double prob;
  Matrix *m = mat_new(hmm->nstates, hmm->nstates);

  /* "normal" transitions */
  for (i = 0; i < hmm->nstates; i++) {
    for (j = 0; j < hmm->nstates; j++) {
      prob = mm_get(hmm->transition_matrix, i, j);
      if (prob == 0) 
	mat_set(m, i, j, NEGINFTY);
      else
	mat_set(m, i, j, log2(prob));
    }
  }
  hmm->transition_score_matrix = m;

  /* end transitions, if used */
  if (hmm->end_transitions != NULL) {
    hmm->end_transition_scores =
      vec_new(hmm->end_transitions->size);
    for (i = 0; i < hmm->end_transitions->size; i++) {
      prob = vec_get(hmm->end_transitions, i);
      if (prob == 0)
	vec_set(hmm->end_transition_scores, i, NEGINFTY);
      else
	vec_set(hmm->end_transition_scores, i, log2(prob));
    }
  }
  
  /* begin transitions */
  hmm->begin_transition_scores =
    vec_new(hmm->begin_transitions->size);
  for (i = 0; i < hmm->begin_transitions->size; i++) {
    prob = vec_get(hmm->begin_transitions, i);
    if (prob == 0)
      vec_set(hmm->begin_transition_scores, i, NEGINFTY);
    else
      vec_set(hmm->begin_transition_scores, i, log2(prob));
  }
}

