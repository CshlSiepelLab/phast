/* $Id: em.c,v 1.1.1.1 2004-06-03 22:43:11 acs Exp $
   Written by Adam Siepel, 2003
   Copyright 2003, Adam Siepel, University of California */

/* Experimental code for training a phylo-HMM by EM, including its
   phylogenetic models. */

#include <em.h>
#include <assert.h>
#include <sufficient_stats.h>
#include <fit_em.h>

/* hmm and models must be initialized appropriately */
/* must be one model for every state in the HMM */
/* the ith training sample in data must be of length 'sample_lens[i]' */
/* returns log likelihood of optimized model */
double hmm_train_by_em(HMM *hmm, void *models, void *data, int nsamples, 
                       int *sample_lens, gsl_matrix *pseudocounts, 
                       void (*compute_emissions)(double**, void**, int, void*, 
                                                 int, int), 
                       void (*estimate_state_models)(void**, int, void*, 
                                                     double**, int),
                       int (*get_observation_index)(void*, int, int)) { 

  int i, k, l, s, obsidx, nobs, maxlen = 0;
  double **emissions, **forward_scores, **backward_scores, **E, **A;
  double *totalE, *totalA;
  double total_logl, prev_total_logl, val;
  List *val_list;

  /* FIXME: rate categories!  currently ignored  (connection with pooling) */

  for (s = 0; s < nsamples; s++)
    if (sample_lens[s] > maxlen) maxlen = sample_lens[s];

  nobs = get_observation_index(data, -1, -1); /* FIXME */

  emissions = (double**)smalloc(hmm->nstates * sizeof(double*));
  forward_scores = (double**)smalloc(hmm->nstates * sizeof(double*));
  backward_scores = (double**)smalloc(hmm->nstates * sizeof(double*));
  for (i = 0; i < hmm->nstates; i++){
    emissions[i] = (double*)smalloc(maxlen * sizeof(double));
    forward_scores[i] = (double*)smalloc(maxlen * sizeof(double));
    backward_scores[i] = (double*)smalloc(maxlen * sizeof(double));
  }
  E = (double**)smalloc(hmm->nstates * sizeof(double*));
  A = (double**)smalloc(hmm->nstates * sizeof(double*));
  totalE = (double*)smalloc(hmm->nstates * sizeof(double));
  totalA = (double*)smalloc(hmm->nstates * sizeof(double));
  for (k = 0; k < hmm->nstates; k++) {
    E[k] = (double*)smalloc(nobs * sizeof(double));
    A[k] = (double*)smalloc(hmm->nstates * sizeof(double));
  }

  val_list = lst_new_dbl(hmm->nstates);

  prev_total_logl = NEGINFTY;
  while (1) {
    total_logl = 0;

    /* initialize 'A' and 'E' counts (see below) */
    for (k = 0; k < hmm->nstates; k++) {
      for (obsidx = 0; obsidx < nobs; obsidx++)
        E[k][obsidx] = 0;
      for (l = 0; l < hmm->nstates; l++)
        A[k][l] = 0;
      totalE[k] = totalA[k] = 0;
    }

    for (s = 0; s < nsamples; s++) {
      double logp_fw, logp_bw;

      compute_emissions(emissions, models, hmm->nstates, data, s, sample_lens[s]);

      logp_fw = hmm_forward(hmm, emissions, sample_lens[s], 
                            forward_scores);
      logp_bw = hmm_backward(hmm, emissions, sample_lens[s], 
                             backward_scores);

      if (abs(logp_fw - logp_bw) > 0.01)
        fprintf(stderr, "WARNING: forward and backward algorithms returned different total log\nprobabilities (%f and %f, respectively).\n", logp_fw, logp_bw);

      total_logl += logp_fw;

      for (i = 0; i < sample_lens[s]; i++) {
        double this_logp;

        /* to avoid rounding errors, estimate total log prob
           separately for each column */
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
          
          /* compute expected number of transitions from each state to
             each other ('A' in Durbin et al.'s notation, pp. 63-64) */
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

    /* normalize As and Es */
    for (k = 0; k < hmm->nstates; k++) {
/*       for (obsidx = 0; obsidx < nobs; obsidx++)  */
/*         E[k][obsidx] /= totalE[k]; */
      for (l = 0; l < hmm->nstates; l++)
        A[k][l] /= totalA[k];
    }

    /* check convergence */
    fprintf(stderr, "Training likelihood: %f\n", total_logl);
    if (total_logl - prev_total_logl <= EM_CONVERGENCE_THRESHOLD)
      break;                    /* FIXME: do this *before* computing
                                   expected counts  */
    prev_total_logl = total_logl;

    /* TODO: update transitions of HMM according to A values, using
       pseudocounts; combine with above? */ 
    for (k = 0; k < hmm->nstates; k++)
      for (l = 0; l < hmm->nstates; l++) 
        mm_set(hmm->transition_matrix, k, l, A[k][l]);
                                /* should all of this be done "in
                                   place" from the beginning?  */
                                /* FIXME: begin and end? */
    /* FIXME: pseudocounts! */
    hmm_reset(hmm);

    /* re-estimate state models */
    estimate_state_models(models, hmm->nstates, data, E, nobs);
    /* TODO: allow for "pools" of models */
  }

  for (i = 0; i < hmm->nstates; i++) {
    free(emissions[i]);
    free(forward_scores[i]);
    free(backward_scores[i]);
    free(E[i]);
    free(A[i]);
  }
  free(emissions);
  free(forward_scores);
  free(backward_scores);
  free(E);
  free(A);
  free(totalE);
  free(totalA);
  lst_free(val_list);

  return total_logl;
}



/* THESE FUNCTIONS SHOULD BE MOVED */

#include <tree_likelihoods.h>

/* implementation of generic routine required by EM code to fill out
   matrix of emissions */
void compute_emissions_phyhmm(double **emissions, void **models, int nmodels,
                              void *data, int sample, int length) {
  int k;
  PooledMSA *pmsa = (PooledMSA*)data;
  for (k = 0; k < nmodels; k++) 
    tl_compute_log_likelihood((TreeModel*)models[k], 
                              lst_get_ptr(pmsa->source_msas, sample),
                              emissions[k], -1, NULL);
}

void estimate_state_models_phyhmm(void **models, int nmodels, void *data, 
                                  double **E, int nobs) {
  int k, obsidx;
  gsl_vector *params;
  PooledMSA *pmsa = (PooledMSA*)data;

  for (k = 0; k < nmodels; k++) {
    TreeModel *tm = (TreeModel*)models[k];
    params = tm_params_init_from_model(tm);
    for (obsidx = 0; obsidx < nobs; obsidx++) 
      pmsa->pooled_msa->ss->cat_counts[k][obsidx] = E[k][obsidx];
/*     msa_get_base_freqs_tuples(pmsa->pooled_msa, tm->backgd_freqs,  */
/*                               tm->order+1, k); */
                                /* need to reestimate background
                                   freqs, using new category counts */
    tm_fit(tm, pmsa->pooled_msa, params, k, -1, OPT_HIGH_PREC, NULL);
                                /* FIXME: should use analytical gradients */
    gsl_vector_free(params); 
  }
}

int get_observation_index_phyhmm(void *data, int sample, int position) {
  PooledMSA *pmsa = (PooledMSA*)data;
  MSA *msa;
  if (sample == -1 || position == -1) 
    return pmsa->pooled_msa->ss->ntuples;
  msa = (MSA*)lst_get_ptr(pmsa->source_msas, sample);
  return msa->ss->tuple_idx[position];
}


