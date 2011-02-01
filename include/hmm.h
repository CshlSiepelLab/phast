/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/**
   @file hmm.h
   Library of functions relating to Hidden Markov Models.  Includes
   simple reading and writing routines, as well as implementations of
   the Viterbi algorithm, the forward algorithm, and the backward
   algorithm.  Also includes a function to compute posterior
   probabilities.

   Emisson Scores - Probabilities of how likely a specific each residue is
   Transition Scores - Probabilities of transitioning to specific states from current one
   Very clear introduction to hmms - http://www.nature.com/nbt/journal/v22/n10/full/nbt1004-1315.html
   @ingroup hmm
*/

#ifndef HMM_H
#define HMM_H

#include <matrix.h>
#include <markov_matrix.h>
#include <misc.h>
#include <lists.h>
#include <vector.h>
#include <prob_vector.h>
/** Maximum number of states in the model */
#define MAXSTATES 1000
/** Used to identify starting state */
#define BEGIN_STATE -99
/** Used to identify finished state */
#define END_STATE -98

#define BEGIN_TRANSITIONS_TAG "BEGIN_TRANSITIONS:"
#define END_TRANSITIONS_TAG "END_TRANSITIONS:"
#define TRANSITION_MATRIX_TAG "TRANSITION_MATRIX:"
#define EQ_FREQS_TAG "EQUILIBRIUM_FREQUENCIES:"

/** Algorithm to use for solving HMMs */
typedef enum {VITERBI, /**< Viterbi for speedy HMM solving */
	      FORWARD, /**< Forward method of posterior decoding*/
	      BACKWARD /**< Backward method of posterior decoding*/
} hmm_mode;

/* NOTE: eventually need to be able to support an "adjacency list"
 * rather than an "adjacency matrix" representation of an HMM, for
 * better efficiency when there are many states and they are not fully
 * connected  */
/** Hidden Markov Model and meta data  */
typedef struct {
  int nstates;  /**< Number of current states in model */
  MarkovMatrix *transition_matrix; /**< Probabilities of moving from one state to another */
  Matrix *transition_score_matrix; /**< Entries are logs of entries
				     in transition matrix. */
  Vector *begin_transitions, /**< Beginning state transition probabilities */
  	 *end_transitions,   /**< Ending state transition probabilities */
    	 *begin_transition_scores,  /**< Entries are log of entries in vector */
	 *end_transition_scores,   /**< Entries are log of entries in vector */
	 *eq_freqs;		   /**< Equilibrium frequencies */
  List **predecessors,		/**< List of predecessor states in HMM, for each state i, the list of states that have a transition to that state */ 
  **successors;			/**< List of successor states in HMM, for each state i, the list of states that state i has a transition to */
  List *begin_successors, /**< List of states for which the begin state has a transition to */
 *end_predecessors;	  /**< List of states that have a transition to the end state */
} HMM;


/** Creates a new HMM object based on a Markov matrix of transition
   probabilities, a vector of transitions from the begin state, and a
   vector of transitions to the end state.  

   All transitions are expected to be described as probabilities, not 
   logs.  The vector of end-state transitions may be NULL, in which 
   case it will be assumed that the end state is not of interest (e.g.,
   the Viterbi and forward/backward algorithms will ignore the end 
   state).  NULL may also be specified for the begin_transitions; in 
   this case, they will be assumed to be uniform.  
   @param mm Matrix where next state depends only on current state
   @param eq_freqs Equilibrium frequencies
   @param begin_transitions Initial transition probabilities
   @param end_transitions Final transition probabilities
   @result Newly created HMM object based on Markov matrix transition probabilities.
 */
HMM* hmm_new(MarkovMatrix *mm, Vector *eq_freqs,
             Vector *begin_transitions, 
             Vector *end_transitions);
/** Create a new HMM object with the specified number of states.
    @param nstates Number of states new HMM object will need
    @param begin if == 1 HMM will include a begin_transaction vector initialized to zero
    @param end if == 1 HMM will include a begin_transaction vector initialized to zero  
    @result New empty hmm
*/
HMM *hmm_new_nstates(int nstates, int begin, int end);

/** Create a new copy of an existing HMM object.
    @param src HMM to copy
    @result clone of source HMM
 */
HMM *hmm_create_copy(HMM *src);

/** Free an HMM object. */
void hmm_free(HMM *hmm);

/** Get the score/probabilities of transitioning from state 'from_state' to state 'to_state'.
    @param hmm Model to get probabilities from
    @param from_state Beginning state to start gathering probabilities from
    @param to_state Final state, where to stop gathering probabilities
    @result Log probability associated with transition between from_state and to_state.
    @note Result includes special states BEGIN_STATE and END_STATE
*/
double hmm_get_transition_score(HMM *hmm, int from_state, int to_state);

/** Create a new HMM object from the contents of a file.
    @param F file containing HMM information

   @note File format is very simple; it includes specification of a transition matrix, a
   vector of transitions from the begin state, and (optionally) a
   vector of transitions to the end state.
   @note The transition matrix
   should be specified using the conventions described in
   markov_matrix.c; the others should simply be lists of
   whitespace-separated numbers.
   @note The description of each object must
   be preceded by a tag, as defined in hmm.h
   @see hmm.h
*/
HMM* hmm_new_from_file(FILE *F);

/** Save textual representation of HMM to a file
    @param F File descriptor to save HMM to
    @param hmm Hmm to get statistics from to write file
*/
void hmm_print(FILE *F, HMM *hmm);

/**
  Finds most probable path, according to the Viterbi algorithm.

  @param[in] hmm Model to use
  @param emission_scores Output scores, 2D array, hmm->nstates rows & columns
  @param[in] seqlen Length of path
  @param[out] path Array of integers indicating state numbers in the HMM
*/
void hmm_viterbi(HMM *hmm, double **emission_scores, int seqlen, int *path);

/** 
   Fills matrix of "forward" scores and returns total log probability
   of sequence. 
   @param[in] hmm Model to use
   @param emission_scores output scores, 2D array, hmm->nstates rows & seqlen columns
   @param[in] seqlen number of columns in emission_scores and forward_scores
   @param[out] foward_scores must be allocated to same size as emission_scores 
   @result total log probability of sequence
*/
double hmm_forward(HMM *hmm, double **emission_scores, int seqlen, 
                   double **forward_scores);
/** 
   Fills matrix of "backward" scores and returns total log probability
   of sequence.
   @param[in] hmm Model to use
   @param emission_scores Output scores, 2D array, hmm->nstates rows * seqlen columns
   @param[in] seqlen Number of columns in emission_scores and backward_scores
   @param[out] backward_scores Must be allocated to same size as emission_scores
   @result Total log probability of sequence
*/
double hmm_backward(HMM *hmm, double **emission_scores, int seqlen,
                    double **backward_scores);

/** Fills matrix of posterior probabilities.
   @param hmm Model to use
   @param emission_scores Output scores, 2D array, hmm->nstates rows & seqlen columns
   @param seqlen Number of columns in emission_scores and posterior_probs
   @param posterior_probs  (Optional) Must be allocated to same size as emission_scores
   @result Total log probability of sequence
*/
double hmm_posterior_probs(HMM *hmm, double **emission_scores, int seqlen,
                           double **posterior_probs);

void hmm_do_dp_forward(HMM *hmm, double **emission_scores, int seqlen, 
                       hmm_mode mode, double **full_scores, int **backptr);
void hmm_do_dp_backward(HMM *hmm, double **emission_scores, int seqlen, 
                        double **full_scores);
/** 
    Finds max or sum of score/transition combination over all previous
   states (max for Viterbi, sum for forward/backward).
   @param[in] hmm Model to use
   @param[in] full_scores Scores for all previous states
   @param[in] emission_scores Emission scores, 2D array
   @param backptr Pointer to predecessor, (Viterbi only)
   @param i Present state (may be END_STATE only if mode == FORWARD) or (BEGIN_STATE only if mode == BACKWARD)
   @param j Present column
   @result Max or sum of score/transition combination over all previous states
*/
double hmm_max_or_sum(HMM *hmm, double **full_scores, double **emission_scores,
                      int **backptr, int i, int j, hmm_mode mode);


void hmm_dump_matrices(HMM *hmm, double **emission_scores, int seqlen,
                       double **full_scores, int **backptr);
/**
   Reset arcs of HMM according to a matrix of counts and optionally,
   a specified matrix of pseudo-counts.

   @param[in] hmm Model to use
   @param[in] trans_counts Count of transitions for each part in matrix
   @param[in] trans_pseudocounts (Optional) Pseudo counts added to trans_counts
   @param[in] state_counts Count of states
   @param[in] beg_counts (Optional) Specify non-uniform bed transitions
   @param[in] bed_pseudocounts (Optional) Pseudo counts added to beg_counts   
*/
void hmm_train_from_counts(HMM *hmm, Matrix *trans_counts, 
                           Matrix *trans_pseudocounts,
                           Vector *state_counts,
                           Vector *state_pseudocounts,
                           Vector *beg_counts, 
                           Vector *beg_pseudocounts);

/* Retrain HMM according to observed path and, optionally, specified
   matrix of pseudo-counts.
   
   @param hmm Model to use
   @param[in] path Path to retrain HMM
   @param[in] npaths Amount of paths provided in parameter 'path'
   @param[in] trans_pseudocount (Optional) Transition pseudo count to adjust transition counts
   @param[in] state_pseudocounts (Optional) State pseudo counts to adjust for state counts
   @param[in] use_begin Whether or not to use beg counts
   @param[in] bed_pseudocounts (Optional) Beg pseudo counts to adjust beg counts
   @note Path indices are expected to range between 1 and hmm->nstates.  Each path is expected to be terminated by -1. 
*/
void hmm_train_from_paths(HMM *hmm, int **path, int npaths,
                          Matrix *trans_pseudocounts, 
                          Vector *state_pseudocounts,int use_begin,
                          Vector *beg_pseudocounts);

/** Update all counts according to a new path 
    @param trans_counts Transition counts to be updated
    @param state_counts State counts to be updated
    @param beg_counts Beg counts to be updated
    @param path New path used to calculate new counts
    @param len Length of path
    @param nstates Number of states
*/
void hmm_train_update_counts(Matrix *trans_counts, Vector *state_counts, 
                             Vector *beg_counts,
                             int *path, int len, int nstates);
/** Create a trivial HMM.
   Starts with a single state that transitions to itself with probability 1.
   @note This if useful for causing general HMM-based models to collapse to simpler models
   @result New trivial HMM
 */
HMM *hmm_create_trivial();

/**  
    Compute the total log likelihood of a specified path.
    @param hmm Model to use
    @param emission_scores Emission scores
    @param seqlen Length of path
    @param path Path to calculate total log likelihood for
 */
double hmm_path_likelihood(HMM *hmm, double **emission_scores, int seqlen, 
                           int *path);

/** Compute the total log likelihood of a subsequence of the input,
   using only the specified states.  Useful for scoring candidate
   predictions.
   @param hmm Model to use
   @param emission_scores Emission scores
   @param states List of indices of states making up the subsequence to compute log likelihood of
   @param begidx Beginning index
   @param len Overall number of states
   @result Log likelihood or subsequence identified by parameter 'states'
 */
double hmm_score_subset(HMM *hmm, double **emission_scores, List *states,
                        int begidx, int len);

/** 
   Compute the log odds score for a subsequence of the input, comparing
   likelihood based on a subset of states (test_states) against (null_states).  
   Useful for scoring candidate predictions. 
   @param hmm Model to use
   @param emission_scores Emission Scores
   @param test_states List of indices of states to compare against
   @param null_states List of indices of states to compare to
   @param begidx Beginning index of comparison
   @param len Overall number of states (must be same for both lists of states)
   return Log dds score between two subsets
   @see hmm_score_subset
*/
double hmm_log_odds_subset(HMM *hmm, double **emission_scores, 
                           List *test_states, List *null_states,
                           int begidx, int len);

/**
   Perform a cross product of two HMMs.  
   @param[out] dest Result of cross product
   @param[in] src1 First matrix in cross product
   @param[in] src2 Second matrix in cross product
   @note The pair of state i from src1 and state j from src2 is represented in the new
   HMM by state i*src2->nstates + j 
*/
void hmm_cross_product(HMM *dest, HMM *src1, HMM *src2);

/**
   Reset various attributes that are derived from the underlying
   matrix of transitions. 
   Should be called after the matrix is changed for any reason.  

   @param hmm Model to reset
   @note This routine assumes that hmm->nstates does not change.
*/
void hmm_reset(HMM *hmm);

/**
  Reverse and complement a sequence strand
Given an HMM, some of whose states represent strand-specific
   phenomena in DNA (e.g., coding regions, UTRs, introns), create a
   new HMM with a second version of all such states, corresponding to
   the reverse strand (assuming the original HMM describes the forward
   strand).  Transition probabilities between reverse-strand states
   will be a "reflection" of those for the forward-strand states,
   based on an assumption of reversibility.
   @param hmm Model to use
   @param pivot_states State indices indicating (in addition to 0) where Transitions could occur between forward- and reverse strand states
   @param mapping Defines correspondence between forward and reverse states.  Pre-allocated size of (hmm->nstates*2 - lst_size(pivot_states) - 1)
*/
HMM *hmm_reverse_compl(HMM *hmm, List *pivot_states, int *mapping);

/**
  Renormalizes transition probabilities so that all rows sum to 1.
  Currently does not consider begin or end transitions.  Calls
   hmm_reset. 
  @param hmm Model to renormalize
*/
void hmm_renormalize(HMM *hmm);

/** Sample a state path through a sequence using the stochastic traceback
   algorithm. 
   @param hmm Model to use
   @param forward_scores Scores using the forward algorithm
   @param seqlen Length of path
   @param path List of steps through the model
 */
void hmm_stochastic_traceback(HMM *hmm, double **forward_scores, 
			      int seqlen, int *path);

/** Set the transition_score_matrix in an hmm object. 
  @param hmm Model to prepare
  @warning This must be done before calling any functions that use hmm_get_transition_score in a multithreaded context.
*/
void hmm_set_transition_score_matrix(HMM *hmm);

#endif
