/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file em.h
    Function for training an phylo-HMM via EM including its models.
    @ingroup hmm
 */

#ifndef EM_H
#define EM_H

#include <hmm.h>
#include <tree_model.h>

//#define EM_CONVERGENCE_THRESHOLD 0.01
#define EM_CONVERGENCE_THRESHOLD 0.1 /* TEMPORARY! */

/** Train a Hidden Markov Model by using EM algorithm.
    @param hmm Hidden Markov Model to train
    @param models ???
    @param data Training data (the ith training sample in data must be of length 'sample_lens[i]' )
    @param nsamples Number of samples
    @param sample_lens Lengths of each data sample
    @param pseudocounts Pseudo counts of ???
    @param compute_emissions (Optional) Function to compute emissions. If NULL simply not called
    @param estimate_state_models (Optional) Function to estimate state models. If NULL fully general parameterization is assumed
    @param estimate_transitions (Optional) Function to estimate transitional probabilities (M step). If NULL fully general is assumed
    @param get_observation_index Function to get observation index
    @param log_function Function to use for logging statistics
    @param emissions_alloc (Optional) Used for emission probabilities (must be large enough for longest sample)
    @param logf Log to save statistics to
    @result Log likelihood of optimized model
    @note HMM and models must be initialized appropriately
    @note Must be one model for every state in the HMM
    @note If sample size is 1, emissions can be pre-computed
    @warning This function is experimental
*/
double hmm_train_by_em(HMM *hmm, void *models, void *data, int nsamples, 
                       int *sample_lens, Matrix *pseudocounts, 
                       void (*compute_emissions)(double**, void**, int, void*, 
                                                 int, int), 
                       void (*estimate_state_models)(TreeModel**, int, void*, 
                                                     double**, int, FILE*),
                       void (*estimate_transitions)(HMM*, void*, double**),
                       int (*get_observation_index)(void*, int, int),
                       void (*log_function)(FILE*, double, HMM*, void*, int),
                       double **emissions_alloc, FILE *logf);

#endif
