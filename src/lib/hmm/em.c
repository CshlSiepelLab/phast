/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: em.c,v 1.12 2008-11-12 02:07:59 acs Exp $ */

/* Experimental code for training a phylo-HMM by EM, including its
   phylogenetic models. */

#include <em.h>
#include <misc.h>
#include <sufficient_stats.h>
#include <fit_em.h>
#include <sys/time.h>

/* generic log function: show log likelihood and all HMM transitions
   probs */
void default_log_function(FILE *logf, double total_logl, HMM *hmm, 
                          void *data, int show_header) {
  int i, j;

  if (show_header) {
    fprintf(logf, "\nlogl\t");
    for (i = 0; i < hmm->nstates; i++) {
      for (j = 0; j < hmm->nstates; j++) {
        fprintf(logf, "(%d,%d)\t", i, j);
      }
    }
    fprintf(logf, "\n");
  }

  fprintf(logf, "%f\t", total_logl);
  for (i = 0; i < hmm->nstates; i++)
    for (j = 0; j < hmm->nstates; j++)
      fprintf(logf, "%f\t", mm_get(hmm->transition_matrix, i, j));
  fprintf(logf, "\n");
  fflush(logf);
}

/* hmm and models must be initialized appropriately */
/* must be one model for every state in the HMM */
/* the ith training sample in data must be of length 'sample_lens[i]' */
/* returns log likelihood of optimized model */
/* if sample size is one, emissions may be precomputed */
/* if estimate_state_models != NULL will be used to re-estimate on
   each iteration; otherwise will estimate transition probs only
   (compute_emissions and get_observation_index will be ignored) */
/* if estimate_transitions != NULL, it will be used for estimating
   transition probs (M step); otherwise a fully general
   parameterization will be assumed */
/* if emissions_alloc is non-NULL, it will be used for emission probs
   (must be large enough for longest sample) */
/* compute_emissions simply won't be called if NULL; this may make
   sense if estimate_state_models == NULL, nsamples == 1, and
   emissions are precomputed & passed in as emissions_alloc */
double hmm_train_by_em(HMM *hmm, void *models, void *data, int nsamples, 
                       int *sample_lens, Matrix *pseudocounts, 
                       void (*compute_emissions)(double**, void**, int, void*, 
                                                 int, int), 
                       void (*estimate_state_models)(void**, int, void*, 
                                                     double**, int, FILE*),
                       void (*estimate_transitions)(HMM*, void*, double**),
                       int (*get_observation_index)(void*, int, int),
                       void (*log_function)(FILE*, double, HMM*, void*, int),
		       double **emissions_alloc, FILE *logf) { 

  int i, k, l, s, obsidx, nobs=0, maxlen = 0, done, it;
  double **emissions, **forward_scores, **backward_scores, **E = NULL, **A;
  double *totalE = NULL, *totalA;
  double total_logl, prev_total_logl, val;
  List *val_list;

  struct timeval start_time, end_time;

  if (estimate_state_models != NULL && 
      (get_observation_index == NULL || compute_emissions == NULL))
    die("ERROR: (hmm_train_by_em) If estimating state models, must pass in non-NULL functions get_observation_index and compute_emissions.\n");

  if (compute_emissions == NULL &&
      (estimate_state_models != NULL || nsamples > 1 || emissions_alloc == NULL))
    die("ERROR: (hmm_train_by_em) compute_emissions function required.\n");
      
  if (logf != NULL)
    gettimeofday(&start_time, NULL);

  for (s = 0; s < nsamples; s++)
    if (sample_lens[s] > maxlen) maxlen = sample_lens[s];

  forward_scores = (double**)smalloc(hmm->nstates * sizeof(double*));
  backward_scores = (double**)smalloc(hmm->nstates * sizeof(double*));

  if (emissions_alloc != NULL)
    emissions = emissions_alloc;
  else 
    emissions = (double**)smalloc(hmm->nstates * sizeof(double*));

  for (i = 0; i < hmm->nstates; i++){
    forward_scores[i] = (double*)smalloc(maxlen * sizeof(double));
    backward_scores[i] = (double*)smalloc(maxlen * sizeof(double));
    if (emissions_alloc == NULL) 
      emissions[i] = (double*)smalloc(maxlen * sizeof(double));
  }
  A = (double**)smalloc(hmm->nstates * sizeof(double*));
  totalA = (double*)smalloc(hmm->nstates * sizeof(double));
  for (k = 0; k < hmm->nstates; k++) 
    A[k] = (double*)smalloc(hmm->nstates * sizeof(double));

  if (estimate_state_models) {
    nobs = get_observation_index(data, -1, -1); /* this is a bit
                                                   clumsy, but will do
                                                   for now */
    E = (double**)smalloc(hmm->nstates * sizeof(double*));
    totalE = (double*)smalloc(hmm->nstates * sizeof(double));
    for (k = 0; k < hmm->nstates; k++) 
      E[k] = (double*)smalloc(nobs * sizeof(double));
  }

  val_list = lst_new_dbl(hmm->nstates);

  prev_total_logl = NEGINFTY;
  done = FALSE;

  for (it = 1; !done; it++) {
    checkInterrupt();
    total_logl = 0;

    /* initialize 'A' and 'E' counts (see below) */
    for (k = 0; k < hmm->nstates; k++) {
      for (l = 0; l < hmm->nstates; l++) A[k][l] = 0;
      totalA[k] = 0;
    }
    if (estimate_state_models != NULL) {
      for (k = 0; k < hmm->nstates; k++) {
        for (obsidx = 0; obsidx < nobs; obsidx++) E[k][obsidx] = 0;
        totalE[k] = 0;
      }
    }

    for (s = 0; s < nsamples; s++) {
      double logp_fw, logp_bw;
      
      if (compute_emissions == NULL || 
	  (estimate_state_models == NULL && nsamples == 1 && it > 1))
	;			/* no need to compute emissions */
      else
	compute_emissions(emissions, models, hmm->nstates, data, 
			  s, sample_lens[s]);

      logp_fw = hmm_forward(hmm, emissions, sample_lens[s], 
                            forward_scores);
      logp_bw = hmm_backward(hmm, emissions, sample_lens[s], 
                             backward_scores);

      if (fabs(logp_fw - logp_bw) > 1.0)
        if (logf != NULL) 
          fprintf(logf, "WARNING: forward and backward algorithms returned different total log\nprobabilities (%f and %f, respectively).\n", logp_fw, logp_bw);

      total_logl += logp_fw;

      for (i = 0; i < sample_lens[s]; i++) {
        double this_logp;

        /* to avoid rounding errors, estimate total log prob
           separately for each column */
        if (estimate_state_models != NULL) {
          lst_clear(val_list);
          for (l = 0; l < hmm->nstates; l++) 
            lst_push_dbl(val_list, (forward_scores[l][i] + 
                                    backward_scores[l][i]));
          this_logp = log_sum(val_list);
          obsidx = get_observation_index(data, s, i);
          for (k = 0; k < hmm->nstates; k++) {
            /* compute expected number of times each state emits each
               distinct observation ('E' in Durbin et al.'s notation; see
               pp. 63-64) */
            val = exp2(forward_scores[k][i] + backward_scores[k][i] - 
                       this_logp);
            E[k][obsidx] += val;
            totalE[k] += val;
          }
        }
           
        /* compute expected number of transitions from each state to
           each other ('A' in Durbin et al.'s notation, pp. 63-64) */
	if (i != sample_lens[s]-1) {
        for (k = 0; k < hmm->nstates; k++) {
          for (l = 0; l < hmm->nstates; l++) {
            val = exp2(forward_scores[k][i] + 
                       hmm_get_transition_score(hmm, k, l) + 
                       emissions[l][i+1] + backward_scores[l][i+1] - 
                       logp_fw);
                                /* FIXME: begin and end states? start
                                   and end idx */
            A[k][l] += val;
            totalA[k] += val;
          }
        }
	}
      }
    }

    if (logf != NULL) {         /* do this before updating params;
                                   otherwise you're outputting the current
                                   likelihood with the new params,
                                   which is confusing */
      if (log_function != NULL)
        log_function(logf, total_logl, hmm, data, it == 1);
      else 
        default_log_function(logf, total_logl, hmm, NULL, it == 1);
    }

    /* check convergence */
    if (fabs(total_logl - prev_total_logl) <= EM_CONVERGENCE_THRESHOLD)
      done = TRUE;              /* no param update */

    else {
      prev_total_logl = total_logl;

      /* update transitions; use special function if given, otherwise
         assume fully general parameterization  */
      if (estimate_transitions != NULL)
        estimate_transitions(hmm, data, A);
      else 
        for (k = 0; k < hmm->nstates; k++)
          for (l = 0; l < hmm->nstates; l++) 
            mm_set(hmm->transition_matrix, k, l, A[k][l] / totalA[k]);

      /* FIXME: begin and end */

      hmm_reset(hmm);

      /* FIXME: need to use pseudocounts here */  

      /* re-estimate state models */
      if (estimate_state_models  != NULL)
        estimate_state_models(models, hmm->nstates, data, E, nobs, logf);
    }
  }

  if (logf != NULL) {
    gettimeofday(&end_time, NULL);
    fprintf(logf, "\nNumber of iterations: %d\nTotal time: %.4f sec.\n", it, 
            end_time.tv_sec - start_time.tv_sec + 
            (end_time.tv_usec - start_time.tv_usec)/1.0e6);
  }

  for (i = 0; i < hmm->nstates; i++) {
    sfree(forward_scores[i]);
    sfree(backward_scores[i]);
    if (emissions_alloc == NULL) sfree(emissions[i]);
    sfree(A[i]);
    if (estimate_state_models != NULL) sfree(E[i]);
  }
  sfree(forward_scores);
  sfree(backward_scores);
  if (emissions_alloc == NULL) sfree(emissions);
  sfree(A);
  sfree(totalA);
  if (estimate_state_models != NULL) {
    sfree(E);
    sfree(totalE);
  }
  lst_free(val_list);

  return total_logl;
}

