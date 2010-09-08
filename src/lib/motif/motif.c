/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: motif.c,v 1.7 2008-11-12 02:07:59 acs Exp $ */

/* Motif-finding routines */

#include "motif.h"
#include "em.h"
#include "msa.h"
#include "time.h"
#include "sufficient_stats.h"
#include "tree_model.h"
#include "tree_likelihoods.h"
#include "stacks.h"
#include "ctype.h"
#include "external_libs.h"
#include "misc.h"

#define DERIV_EPSILON 1e-6

/* this can be used to override analytical computation of derivatives */
#define NUMERICAL_DERIVS 0
/* #define NUMERICAL_DERIVS 1 */

/* functions passed to general EM code (see below) */
/* (single sequences) */
void mn_compute_emissions(double **emissions, void **models, int nmodels,
                          void *data, int sample, int length);
void mn_estim_mods(void **models, int nmodels, void *data, double **E, 
                   int nobs);
int mn_get_obs_idx(void *data, int sample, int position);

/* (multi-sequences) */
void phy_compute_emissions(double **emissions, void **models, int nmodels,
                           void *data, int sample, int length);
void phy_estim_mods(void **models, int nmodels, void *data, double **E, 
                    int nobs);
int phy_get_obs_idx(void *data, int sample, int position);

/* function passed to opt_bfgs in discriminative training (see below) */
double mtf_compute_conditional(Vector *params, void *data);
void mtf_compute_conditional_grad(Vector *grad, Vector *params, 
                                  void *data, Vector *lb, Vector *ub);

/* comparison function for sorting motifs in descending order by log
   likelihood (see below) */
int score_compare(const void* ptr1, const void* ptr2) {
  double val1 = (*((Motif**)ptr1))->score;
  double val2 = (*((Motif**)ptr2))->score;
  if (val2 > val1) return 1;
  else if (val2 < val1) return -1;
  return 0;
}

/* Find motifs in a collection individual sequences or multiple
   alignments, either using EM or discriminative training.  If
   'multiseq' == 1 then 'data' must be a PooledMSA object; otherwise,
   it must be a SeqSet object.  Return value is a list of Motif
   objects.  The best 'nmotif' motifs of size 'motif_size' are
   selected.  The 'backgd' argument indicates a background model,
   either a TreeModel (if multiseq == 1) or a Vector (if multiseq
   == 0).  If 'has_motif' is non-NULL, then motifs are learned by
   discriminative training, with 'has_motif' interpreted as an array
   of values indicating whether each training example should be
   considered a positive (1) or negative (0) example, or somewhere in
   between.  In this case, 'backgd' is only used for initialization.
   The 'prior' argument indicates an initial value for the prior
   probability that a motif instance appears in each sequence (used
   with EM only).  See calling code in phast_motif.c regarding
   'init_list,' 'sample_parms,' and 'npseudocounts.' */
List* mtf_find(void *data, int multiseq, int motif_size, int nmotifs, 
               TreeNode *tree, void *backgd, double *has_motif, double prior, 
               int nrestarts, List *init_list, int sample_parms, 
               int npseudocounts) {

  int i, j, k, cons, trial, alph_size, nparams = -1;
  double *alpha;
  List *motifs = lst_new_ptr(nrestarts * 
                             (init_list != NULL ? lst_size(init_list) : 1));
  List *tmpl;
  char *cons_str = smalloc((motif_size + 1) * sizeof(char));
  SeqSet *seqset = !multiseq ? data : NULL;
  PooledMSA *pmsa = multiseq ? data : NULL;
  Vector **freqs = smalloc((motif_size + 1) * sizeof(void*));
  Vector *params = NULL, *lower_bounds = NULL, *upper_bounds = NULL;
  int *inv_alphabet = multiseq ? pmsa->pooled_msa->inv_alphabet :
    seqset->set->inv_alphabet;
  Hashtable *hash;

  cons_str[motif_size] = '\0';
  alph_size = multiseq ? strlen(pmsa->pooled_msa->alphabet) : 
    strlen(seqset->set->alphabet);
  alpha = smalloc(alph_size * sizeof(double));
  for (i = 0; i < alph_size; i++) alpha[i] = 1;      
  for (i = 0; i <= motif_size; i++) freqs[i] = vec_new(alph_size);
          
  if (has_motif == NULL) {      /* non-discriminative case only */
    if (multiseq) 
      vec_copy(freqs[0], ((TreeModel*)backgd)->backgd_freqs);
    else 
      vec_copy(freqs[0], backgd);
  }

  for (cons = 0; 
       cons < (init_list == NULL ? 1 : lst_size(init_list)); 
       cons++) {              /* (loop only once if no init_list) */

    String *initstr = init_list == NULL ? NULL : 
      lst_get_ptr(init_list, cons);

    for (trial = 0; trial < nrestarts; trial++) {
      Motif *m;

      if (nrestarts == 1)
        fprintf(stderr, "Trying candidate %d ... ", cons+1);
      else 
        fprintf(stderr, "Trying candidate %d, trial %d ... ", 
                cons+1, trial+1);

      if (initstr == NULL)
        for (i = 1; i <= motif_size; i++) 
          mtf_draw_multinomial(freqs[i], alpha);

      else
        mtf_init_from_consensus(initstr, freqs, inv_alphabet,
                                npseudocounts, sample_parms, motif_size);

      /* create a new motif object */
      m = multiseq ? 
        mtf_new(motif_size, 1, freqs, pmsa, backgd, 0.25) :
        mtf_new(motif_size, 0, freqs, seqset, NULL, 0);

      /* now train */
      if (has_motif == NULL) {  /* EM training */
        if (multiseq)           
          m->score = mtf_em(m->ph_mods, pmsa, lst_size(pmsa->source_msas), 
                           pmsa->lens, m->motif_size, prior, 
                           phy_compute_emissions, phy_estim_mods, 
                           phy_get_obs_idx, m->postprob, m->bestposition);
        else 
          m->score = mtf_em(m->freqs, seqset, seqset->set->nseqs, seqset->lens, 
                           m->motif_size, prior, mn_compute_emissions, 
                           mn_estim_mods, mn_get_obs_idx, m->postprob,
                           m->bestposition);
      }
      else {                    /* discriminative training */
        double retval;
        m->has_motif = has_motif;

        if (nparams == -1) {    /* first time only */
          int params_per_model = multiseq ? 
            tm_get_nparams(m->ph_mods[1]) : /* assume all are the same */
            m->alph_size;
    
          nparams = params_per_model * m->motif_size + 1;
                                /* one more for motif threshold */
          params = vec_new(nparams);
          lower_bounds = vec_new(nparams); 
          vec_set_all(lower_bounds, 0.00001);
          upper_bounds = vec_new(nparams);
          vec_set_all(upper_bounds, 1);
          vec_set(lower_bounds, 0, NEGINFTY); /* threshold */
          vec_set(upper_bounds, 0, INFTY);
          /* no upper bounds */
        }

        /* initialize params */
        j = 0;
        vec_set(params, j++, 2 * motif_size);
                                /* approx 2 nats per model seems to be
                                   a reasonable initialization for the
                                   threshold */
        for (i = 1; i <= m->motif_size; i++) {
          if (multiseq) {
            Vector *tm_params = tm_params_new_init_from_model(m->ph_mods[i]);
/*             vec_set(upper_bounds, j, 20); */ /* FIXME: have to relax upper bound for rate constant */
/*             vec_set(lower_bounds, j, .25); */ /* FIXME: avoid degenerate case */
            for (k = 0; k < tm_params->size; k++)
              vec_set(params, j++, vec_get(tm_params, k));
            vec_free(tm_params);
          }
          else 
            for (k = 0; k < m->alph_size; k++)
              vec_set(params, j++, vec_get(m->freqs[i], k));
        }
	if (j != nparams)
	  die("ERROR mtf_find j (%i) != nparams (%i)\n", j, nparams);
          
        retval = opt_bfgs(mtf_compute_conditional, params, m, &m->score, 
                          lower_bounds, upper_bounds, NULL,
                          NUMERICAL_DERIVS ? NULL : 
                          mtf_compute_conditional_grad, 
                          OPT_LOW_PREC, NULL);

        m->score *= -1;

        if (retval != 0) 
          /* (the opt_bfgs code produces an error message) */
          fprintf(stderr, " ... continuing ... ");
      }      

      mtf_get_consensus(m, cons_str);
      fprintf(stderr, "(consensus = '%s', score = %.3f)\n", cons_str, m->score);

      mtf_predict(m, m->training_data, m->bestposition, m->samplescore, 
                  has_motif);   /* predict and score best motif */

      lst_push_ptr(motifs, m);
    }
  }

  lst_qsort(motifs, score_compare);

  /* for now, report no more than one motif per consensus
     (best-scoring version) */
  hash = hsh_new(lst_size(motifs));
  tmpl = lst_new_ptr(nmotifs);
  for (i = 0; i < lst_size(motifs); i++) {
    Motif *m = lst_get_ptr(motifs, i);
    mtf_get_consensus(m, cons_str);
    if (hsh_get_int(hash, cons_str) == -1 && lst_size(tmpl) < nmotifs) {
      lst_push_ptr(tmpl, m);
      hsh_put_int(hash, cons_str, 1);
    }
    else mtf_free(m);
  }
  hsh_free(hash);
  lst_free(motifs);
  motifs = tmpl;

  for (i = 0; i <= motif_size; i++) vec_free(freqs[i]);
  free(freqs);

  if (params != NULL) vec_free(params);
  if (lower_bounds != NULL) vec_free(lower_bounds);
  if (upper_bounds != NULL) vec_free(upper_bounds);
  free(cons_str);
  free(alpha);

  return motifs;
}

/* implementations of functions passed to EM code */

/* multinomial: compute emissions for all models for the given sample */
void mn_compute_emissions(double **emissions, void **models, int nmodels,
                          void *data, int sample, int length) {
  int i, j;
  MSA *set = ((SeqSet*)data)->set;
  Vector *logmod = vec_new(((Vector*)models[0])->size);
  for (i = 0; i < nmodels; i++) {
    Vector *thismod = (Vector*)models[i];
    if (emissions[i] == NULL) continue;
                                /* can force certain models to be
                                   ignored by setting to NULL */
    if (thismod->size != ((Vector*)models[0])->size)
      die("ERROR mn_compute_emissions bad dimensions\n");
    for (j = 0; j < thismod->size; j++) {
      double val = vec_get(thismod, j);
      vec_set(logmod, j, val == 0 ? NEGINFTY : log(val));
    }
    for (j = 0; j < length; j++) {
      emissions[i][j] = 
        vec_get(logmod, 
                       set->inv_alphabet[(int)set->seqs[sample][j]]);
    }
  }    
  vec_free(logmod);
}

/* multinomial: M step for models.  Based on expected counts, update
   model parameters.  In this case, we just base a new multinomial
   distrib. on the expected counts */
void mn_estim_mods(void **models, int nmodels, void *data, 
                   double **E, int nobs) {
  int i, j;
  for (i = 1; i < nmodels; i++) { /* skip the backgd model */
    double sum = 0;
    for (j = 0; j < nobs; j++) { 
      double val = max(E[i][j], MTF_EPSILON); /* don't let it go zero */
      vec_set(models[i], j, val);                                
      sum += val;
    }
    vec_scale(models[i], 1/sum);
  }
}

/* multinomial: return index for observation; in this case it's just
   the inv_alphabet for the specified sequence and position */
int mn_get_obs_idx(void *data, int sample, int position) {
  MSA *set = ((SeqSet*)data)->set;
  if (sample < 0 || position < 0) return strlen(set->alphabet);
                                /* this hack lets the EM code figure
                                   out the size of state space */
  return set->inv_alphabet[(int)set->seqs[sample][position]];
}

/* phylogenetic: compute emissions for all models for the given sample */
void phy_compute_emissions(double **emissions, void **models, int nmodels,
                           void *data, int sample, int length) {
  int i, k;
  PooledMSA *pmsa = (PooledMSA*)data;
  double conv_factor = log(2);
  for (k = 0; k < nmodels; k++) {
    if (emissions[k] != NULL) { /* can force certain models to be
                                   ignored by setting to NULL */
      MSA *smsa = lst_get_ptr(pmsa->source_msas, sample);
      tl_compute_log_likelihood((TreeModel*)models[k], smsa, 
                                emissions[k], -1, NULL);
      for (i = 0; i < smsa->length; i++) emissions[k][i] *= conv_factor;
                                /* convert to natural log scale */
    }
  }
}

/* phylogenetic: based on expected counts, update model parameters (M step) */
void phy_estim_mods(void **models, int nmodels, void *data, 
                    double **E, int nobs) {
  int k, obsidx;
  Vector *params;
  PooledMSA *pmsa = (PooledMSA*)data;

  for (k = 1; k < nmodels; k++) {
    TreeModel *tm = (TreeModel*)models[k];
    params = tm_params_new_init_from_model(tm);
    for (obsidx = 0; obsidx < nobs; obsidx++) 
      pmsa->pooled_msa->ss->cat_counts[k][obsidx] = E[k][obsidx];
    tm_fit(tm, pmsa->pooled_msa, params, k, OPT_HIGH_PREC, NULL);
/*     fprintf(stderr, "%d: %f, %f\n", k, tm->scale, tr_total_len(tm->tree)); */
    vec_free(params); 
  }
}

/* phylogenetic: return index for observation; in this case it's the
   tuple index */
int phy_get_obs_idx(void *data, int sample, int position) {
  PooledMSA *pmsa = (PooledMSA*)data;
  MSA *msa;
  if (sample == -1 || position == -1) 
    return pmsa->pooled_msa->ss->ntuples;
  msa = (MSA*)lst_get_ptr(pmsa->source_msas, sample);
  return pmsa->tuple_idx_map[sample][msa->ss->tuple_idx[position]];
}

/* A little package of intermediate computations from
   mtf_compute_conditional that can be reused in
   mtf_compute_conditional_grad */
/* typedef struct { */
/*   double *l;   */                  /* l[s] is a kind of log motif probability
                                   for sample s: it's the log of the
                                   sum over all starting positions of
                                   the probability of a motif instance
                                   (see details below) */
  /* this will do for now; later, might try to store emissions per
     observation type per model */
/* } OptData; */

/* this is the function that is optimized in discriminative training;
   see Segal et al., RECOMB '02 */
double mtf_compute_conditional(Vector *params, void *data) {
  Motif *m = (Motif*)data;
  PooledMSA *pmsa = m->multiseq ? m->training_data : NULL;
  SeqSet *ss = m->multiseq ? NULL : m->training_data;
  int *sample_lens = m->multiseq ? pmsa->lens : ss->lens;
  int nsamples = m->multiseq ? lst_size(pmsa->source_msas) : ss->set->nseqs;
  int i, j, k, s;
  double threshold, tot_logprob; 
  double **emissions = smalloc((m->motif_size+1) * sizeof(double*));
  List *tmplist;
  int maxlen = 0;

  for (i = 0; i < nsamples; i++)
    if (sample_lens[i] > maxlen) maxlen = sample_lens[i];

  emissions[0] = NULL;          /* causes background emissions not to
                                   be computed */
  for (i = 1; i <= m->motif_size; i++)
    emissions[i] = (double*)smalloc(maxlen * sizeof(double));

  tmplist = lst_new_dbl(maxlen);

  /* unpack params */
  j = 0;
  threshold = vec_get(params, j++); 
  /* (threshold can be interpreted as log prior odds for a motif instance) */

  for (i = 1; i <= m->motif_size; i++) {
    if (m->multiseq) {
      tm_unpack_params(m->ph_mods[i], params, j);
      j += tm_get_nparams(m->ph_mods[i]);
    }
    else {
      double sum = 0;
      for (k = 0; k < m->alph_size; k++) sum += vec_get(params, j+k);
      for (k = 0; k < m->alph_size; k++)
        vec_set(m->freqs[i], k, vec_get(params, j++)/sum);
/*         vec_set(m->freqs[i], k, vec_get(params, j++)); */
    }
  }
  
  /* compute conditional log probability */
  tot_logprob = 0;
  for (s = 0; s < nsamples; s++) {
    double sample_logprob, lsum;

    /* first compute emissions */
    if (m->multiseq)
      phy_compute_emissions(emissions, (void**)m->ph_mods, m->motif_size+1, 
                            (void*)pmsa, s, pmsa->lens[s]);
    else
      mn_compute_emissions(emissions, (void**)m->freqs, m->motif_size+1, 
                           (void*)ss, s, ss->lens[s]);

    lst_clear(tmplist);
    for (i = 0; i < sample_lens[s] - m->motif_size; i++) { /* FIXME: +1 */
      double val = 0;
      for (j = 0; j < m->motif_size; j++) 
        val += emissions[j+1][i+j];
      lst_push_dbl(tmplist, val);
    }

    lsum = threshold - log(sample_lens[s] - m->motif_size) + log_sum_e(tmplist);
    /* this is the quantity inside the logit in equation 2, Segal et
       al., RECOMB '02.  But note that there's an error in that paper:
       the quantity 'v' either has to be defined as the prior odds
       (without the log) or has to appear outside the log in equation 2 */

    /* now, we want the log of the logit of the log sum (conditional log
       probability).  We also want to interpolate between the case in
       which has_motif[s] == 1 and has_motif[s] == 0. */

    /* the log prob of the sample is the log of a weighted sum of
       logit(logsum) and 1-logit(logsum) */ 

/*     sample_logprob = log((m->has_motif[s] + (1-m->has_motif[s]) * exp(-lsum))  */
/*                           (1 + exp(-lsum)));  */

    /* (do it this way to avoid underflow) */
    if (m->has_motif[s] == 0)
      sample_logprob = -lsum - log(1 + exp(-lsum));
    else if (m->has_motif[s] == 1)
      sample_logprob = -log(1 + exp(-lsum)); /* just for efficiency */
    else 
      sample_logprob = log(m->has_motif[s]) - log(1 + exp(-lsum)) + 
        log(1 + (1-m->has_motif[s])/m->has_motif[s] * exp(-lsum));

    if (!(finite(sample_logprob)))
      die("ERROR mtf_compute_conditional: sample_logprob not finite\n");

    tot_logprob += sample_logprob;
  }

  for (i = 1; i <= m->motif_size; i++) free(emissions[i]);
  free(emissions);
  lst_free(tmplist);

  return -tot_logprob;          /* negate because opt_bfgs minimizes */
}

/* functions used by mtf_compute_conditional_grad (see below) */
void mtf_compute_inner_derivs_phy(double **derivs, Motif *m, 
                                  Vector *params);
void mtf_compute_inner_derivs_mn(double **derivs, Motif *m, 
                                 Vector *params);

/* compute gradients for mtf_compute_conditional */
void mtf_compute_conditional_grad(Vector *grad, Vector *params, 
                                  void *data, Vector *lb, Vector *ub) {
  Motif *m = (Motif*)data;
  PooledMSA *pmsa = m->multiseq ? m->training_data : NULL;
  SeqSet *ss = m->multiseq ? NULL : m->training_data;
  int *sample_lens = m->multiseq ? pmsa->lens : ss->lens;
  int nsamples = m->multiseq ? lst_size(pmsa->source_msas) : ss->set->nseqs;
  int i, j, k, s, params_per_model, nobs, maxlen = 0;
  double **emissions = smalloc((m->motif_size+1) * sizeof(double*));
  double **Dp;
  int *param_to_model, *pos_to_type;
  double *motif_lprob;
  List *tmplist;
  double threshold;

  for (i = 0; i < nsamples; i++)
    if (sample_lens[i] > maxlen) maxlen = sample_lens[i];

  emissions[0] = NULL;          /* causes background emissions not to
                                   be computed */
  for (i = 1; i <= m->motif_size; i++)
    emissions[i] = (double*)smalloc(maxlen * sizeof(double));
  pos_to_type = smalloc(maxlen * sizeof(int));
  motif_lprob = smalloc(maxlen * sizeof(double));
  tmplist = lst_new_dbl(maxlen);

  /* note: we're assuming the motif models are already consistent with
     the parameters vector */

  /* construct a mapping from parameters to models */
  params_per_model = m->multiseq ? 
    tm_get_nparams(m->ph_mods[1]) : 
    m->alph_size;  
  param_to_model = smalloc((m->motif_size * params_per_model + 1) * sizeof(int));
  param_to_model[0] = -1;
  for (i = 0; i < m->motif_size; i++)
    for (j = 0; j < params_per_model; j++)
      param_to_model[i*params_per_model + j + 1] = i + 1;

  nobs = m->multiseq ? pmsa->pooled_msa->ss->ntuples : m->alph_size;

  /* compute inner derivatives.  Let Dp[i][j] be the partial
     derivative wrt parameter i of log P(X_j|theta(i)), where X_j is
     the jth *type* of observation, and theta(i) is the motif model
     associated with param i (a many-to-one relationship is assumed
     between parameters and models).  The threshold parameter (index
     0) is a special case; it's not associated with any motif model
     and is ignored here.  */
  threshold = vec_get(params, 0);
  Dp = smalloc(params->size * sizeof(void*));
  for (i = 1; i < params->size; i++) Dp[i] = smalloc(nobs * sizeof(double));
  if (m->multiseq) mtf_compute_inner_derivs_phy(Dp, m, params);
  else mtf_compute_inner_derivs_mn(Dp, m, params);

  /* now compute the gradient, sample by sample and parameter by
     parameter */
  vec_zero(grad);
  for (s = 0; s < nsamples; s++) {
    double l, c, factor;

    /* compute emissions; also define mapping of observations to
       observation types */
    if (m->multiseq) {
      MSA *source_msa;
      phy_compute_emissions(emissions, (void**)m->ph_mods, m->motif_size+1, 
                            (void*)pmsa, s, pmsa->lens[s]);
      source_msa = lst_get_ptr(pmsa->source_msas, s);
      for (j = 0; j < sample_lens[s]; j++)
        pos_to_type[j] = pmsa->tuple_idx_map[s][source_msa->ss->tuple_idx[j]];
    }
    else {
      mn_compute_emissions(emissions, (void**)m->freqs, m->motif_size+1, 
                           (void*)ss, s, ss->lens[s]);
      for (j = 0; j < sample_lens[s]; j++)
        pos_to_type[j] = ss->set->inv_alphabet[(int)ss->set->seqs[s][j]];
    }

    /* we need the (log) prob of a motif instance starting at each
       position in the sample */
    lst_clear(tmplist);
    for (j = 0; j < sample_lens[s] - m->motif_size; j++) { /* FIXME: +1 */
      double logprob = 0;
      for (k = 0; k < m->motif_size; k++) 
        logprob += emissions[k+1][j+k];
      lst_push_dbl(tmplist, logprob);
      motif_lprob[j] = logprob;
    }

    /* intermediate quantities used in derivative calculations */
    l = log_sum_e(tmplist);
    c = -log(sample_lens[s] - m->motif_size);
    factor = (1 - logit(threshold + c + l)) * (2 * m->has_motif[s] - 1) / 
      (m->has_motif[s] + (1 - m->has_motif[s]) * exp(-threshold - c - l));
                                /* this factor appears in all partial
                                   derivatives */

    /* first take care of the threshold parameter (index 0).  In this
       case, the partial derivative is just the factor computed above */
    vec_set(grad, 0, vec_get(grad, 0) + factor);

    /* now compute the partial derivatives for the other parameters:
       each one is essentially a sum over all starting positions of
       the inner derivatives weighted by the post prob of a motif */
    for (i = 1; i < params->size; i++) {
      double deriv = 0;         /* partial derivative for this 
                                   parameter and this sample */
      int model = param_to_model[i];

      for (j = 0; j < sample_lens[s] - m->motif_size; j++) /* FIXME: +1 */
        if (motif_lprob[j] > NEGINFTY) /* avoid 0 * infty situation */
          deriv += exp(motif_lprob[j] - l) * Dp[i][pos_to_type[j+model-1]];
      deriv *= factor;

      vec_set(grad, i, vec_get(grad, i) + deriv);
    }
  }

  vec_scale(grad, -1);   /* because function is negated */

  free(param_to_model);
  free(pos_to_type);
  free(motif_lprob);
  lst_free(tmplist);
  for (i = 1; i < params->size; i++) free(Dp[i]);
  free(Dp);
  for (i = 1; i <= m->motif_size; i++) free(emissions[i]);
  free(emissions);
}

/* compute "inner derivatives" for a phylogenetic modif model (see
   mtf_compute_conditional_grad).  derivs[i][j] will be set to the
   (numerically estimated) partial derivative of the log probability
   of the jth type of observation given the model associated with the
   ith parameter */
void mtf_compute_inner_derivs_phy(double **derivs, Motif *m, 
                                  Vector *params) {
  int i, j, k, mod, base_param_idx;
  PooledMSA *pmsa = m->training_data;
  MSA *dummy_msa;
  double *col_ll, *col_ll_tweak;
  double conv_factor = log(2);  /* need to use natural log scale */

  /* create a dummy alignment based on the pmsa, so that the
     existing function for likelihood computation can be used */
  dummy_msa = msa_create_copy(pmsa->pooled_msa, 1);
  dummy_msa->length = dummy_msa->ss->ntuples;
  dummy_msa->ss->tuple_idx = smalloc(dummy_msa->length * sizeof(int));
  for (i = 0; i < dummy_msa->length; i++) {
    dummy_msa->ss->counts[i] = 1; /* one version of each tuple */
    dummy_msa->ss->tuple_idx[i] = i;
  }

  col_ll = smalloc(dummy_msa->length * sizeof(double));
  col_ll_tweak = smalloc(dummy_msa->length * sizeof(double));
  base_param_idx = 1;           /* base param for first model (skip param 0) */
  for (mod = 1; mod <= m->motif_size; mod++) {
    int nparams_this_model = tm_get_nparams(m->ph_mods[mod]);
    tm_unpack_params(m->ph_mods[mod], params, base_param_idx);
    tl_compute_log_likelihood(m->ph_mods[mod], dummy_msa, col_ll, -1, NULL);
    for (k = 0; k < nparams_this_model; k++) {
      double origparm;
      i = base_param_idx + k;   /* index of parameter of interest */
      origparm = vec_get(params, i);
      vec_set(params, i, origparm + DERIV_EPSILON);
      tm_unpack_params(m->ph_mods[mod], params, base_param_idx);
      tl_compute_log_likelihood(m->ph_mods[mod], dummy_msa, col_ll_tweak, 
                                -1, NULL);

      for (j = 0; j < dummy_msa->ss->ntuples; j++)
        derivs[i][j] = conv_factor * (col_ll_tweak[j] - col_ll[j]) / 
          DERIV_EPSILON;
      
      vec_set(params, i, origparm);
    }
    base_param_idx += nparams_this_model; /* base param for next model */
  }

  msa_free(dummy_msa);
  free(col_ll);
  free(col_ll_tweak);
}

/* compute "inner derivatives" for a multinomial modif model (see
   mtf_compute_conditional_grad).  derivs[i][j] will be set to the
   partial derivative of the log probability
   of the jth type of observation given the model associated with the
   ith parameter */
void mtf_compute_inner_derivs_mn(double **derivs, Motif *m, 
                                 Vector *params) {
  int i, j, model;
  for (model = 1; model <= m->motif_size; model++) {
    double sum = 0;
    int base_idx = (model-1) * m->alph_size + 1;
    for (i = 0; i < m->alph_size; i++) 
      sum += vec_get(params, base_idx + i);

    for (i = 0; i < m->alph_size; i++) {
      double thisparam = vec_get(params, base_idx + i);
      for (j = 0; j < m->alph_size; j++) {
        if (i == j)
          derivs[base_idx+i][j] = (thisparam == 0 ?
                                   INFTY : 
                                   1/thisparam - 1/sum);
        else 
          derivs[base_idx+i][j] = -1/sum;
      }
    }
  }
}

/* estimate a (multinomial) background model from a set of sequences */
void mtf_estim_backgd_mn(SeqSet *s, Vector *model) {
  double count = 0;
  int i, j;
  for (i = 0; i < s->set->nseqs; i++) {
    for (j = 0; j < s->lens[i]; j++) {
      int idx = s->set->inv_alphabet[(int)s->set->seqs[i][j]];
      if (idx < 0) continue;
      vec_set(model, idx, vec_get(model, idx) + 1);
      count++;
    }
  }
  vec_scale(model, 1.0/count);
}

/* find a single motif by EM, given a pre-initialized set of models.
   Functions must be provided for computing "emission" probabilities
   under all models, for updating model parameters given posterior
   expected counts, and for indexing observations.  It is assumed that
   the motif appears at most once in each sequence of observations.
   This function can be used with phylogenetic models or ordinary
   multinomial models.  The models that are passed in are updated, and
   at convergence, represent the (apparent) m.l.e.  The return value
   is the maximized log likelihood.  The number of models is assumed
   to be width+1, where width is the motif width.  The first model is
   assumed to represent the background distribution and its parameters
   will not be updated.  A prior probability that each sequence has a
   motif must be specified.  If non-NULL, postprob and bestposition
   are expected to be arrays of size nsamples; they will be populated
   with values indicating, respectively for each sample, the posterior
   prob. that a motif appears, and the starting position of the best
   instance of the motif.  */
double mtf_em(void *models, void *data, int nsamples, 
              int *sample_lens, int width, double motif_prior,
              void (*compute_emissions)(double**, void**, int, void*, 
                                        int, int), 
              void (*estimate_state_models)(void**, int, void*, 
                                            double**, int),
              int (*get_observation_index)(void*, int, int),
              double *postprob, int *bestposition) {
  
  int i, j, k, s, obsidx, nobs, maxlen = 0;
  double **emissions, **E;
  double *logpY, *postpY;
  double total_logl, prev_total_logl, expected_nmotifs, max=0, window_sum;
  List *tmplst;

  for (s = 0; s < nsamples; s++)
    if (sample_lens[s] > maxlen) maxlen = sample_lens[s];

  nobs = get_observation_index(data, -1, -1); /* convention is to
                                                 return total number
                                                 in this case */

  emissions = (double**)smalloc((width+1) * sizeof(double*));
  for (i = 0; i <= width; i++)
    emissions[i] = (double*)smalloc(maxlen * sizeof(double));
  E = (double**)smalloc((width+1) * sizeof(double*));
  for (k = 1; k <= width; k++) 
    E[k] = (double*)smalloc(nobs * sizeof(double));
  tmplst = lst_new_dbl(maxlen);
  logpY = smalloc(maxlen * sizeof(double));
  postpY = smalloc(maxlen * sizeof(double));

  prev_total_logl = NEGINFTY;
  while (1) {
    total_logl = expected_nmotifs = 0;
    for (k = 1; k <= width; k++) 
      for (obsidx = 0; obsidx < nobs; obsidx++)
        E[k][obsidx] = 0;

    for (s = 0; s < nsamples; s++) {
      double tot_ll_backgd, tot_ll_motif, sample_logl, postpZ;

      compute_emissions(emissions, models, width+1, data, s, sample_lens[s]);

      tot_ll_backgd = 0;
      for (i = 0; i < sample_lens[s]; i++)
        tot_ll_backgd += emissions[0][i]; /* log likelihood of backgd
                                             model (for this sample);
                                             for the moment, leave out
                                             the prior */

      /* let the ith element of logpY be the log of the joint (prior)
         probability that there is a motif and it starts at position i
         in the current sample */
      lst_clear(tmplst);
      for (i = 0; i < sample_lens[s] - width; i++) {
        logpY[i] = log(motif_prior/(sample_lens[s] - width)) + tot_ll_backgd;
        /* (use backgd model for all positions but motif; motif
           positions factored out below) */
        for (j = 0; j < width; j++)
          logpY[i] += emissions[j+1][i+j] - emissions[0][i+j];
        lst_push_dbl(tmplst, logpY[i]);
      }
      tot_ll_motif = log_sum_e(tmplst); /* log likelihood of motif model
                                         (sum over all starting points
                                         for the motif) */

      /* now put in the prior for the backgd model */
      tot_ll_backgd += log(1-motif_prior);

      if (tot_ll_motif < NEGINFTY) 
        sample_logl = tot_ll_backgd;
      else
        sample_logl = tot_ll_motif + log(1 + exp(tot_ll_backgd - tot_ll_motif));
                                /* do it this way to avoid underflow */
      
      if (!(finite(sample_logl)))
	die("ERROR mtf_em sample_logl not finite\n");

      /* now let postpY[i] be the posterior probability that there is a
         motif and it starts at position i */
      if (bestposition != NULL) { bestposition[s] = -1; max = 0; }
      window_sum = 0;
      for (i = 0; i < sample_lens[s] - width; i++) {
        postpY[i] = exp(logpY[i] - sample_logl);

        /* this is a hack used by MEME to avoid giving preference to
           repetitive motifs: force sum to be at most one within each 
           window of size width */
        window_sum += postpY[i];
        if (i >= width) window_sum -= postpY[i-width];
        if (window_sum > 1) {
          for (j = max(0, i-width+1); j <= i; j++)
            postpY[j] /= window_sum;
          window_sum = 1;
        }
      }

      /* have to do this on a separate pass because of the
         scaling hack */
      for (i = 0; i < sample_lens[s] - width; i++) {
        if (bestposition != NULL && postpY[i] > max) { 
          bestposition[s] = i;
          max = postpY[i];
        }
      }

      /* postpZ is the posterior probability that there is a motif in
         this sample (a sum over all postpYs) */
      postpZ = exp(tot_ll_motif - sample_logl);
      expected_nmotifs += postpZ;
      if (postprob != NULL) postprob[s] = postpZ;
      
      total_logl += sample_logl; /* running total across samples */

      /* now update expected numbers of each type of character
         generated by each motif state */
      for (i = 0; i < sample_lens[s] - width; i++) {
        for (k = 0; k < width; k++) {
          obsidx = get_observation_index(data, s, i+k);
          E[k+1][obsidx] += postpY[i];
        }
      }
    }

    /* check convergence */
/*     fprintf(stderr, "Training likelihood: %f\n", total_logl); */

    if (total_logl - prev_total_logl <= MTF_EM_CONVERGENCE_THRESHOLD)
      break;              

    prev_total_logl = total_logl;

    /* re-estimate state models */
    estimate_state_models(models, width+1, data, E, nobs);

    /* update motif prior */
    motif_prior = min(1-MTF_EPSILON, expected_nmotifs/nsamples);
                                /* don't let it go quite to 1 */
  }

  for (i = 0; i <= width; i++) {
    free(emissions[i]);
    if (i > 0) free(E[i]);
  }
  free(emissions);
  free(E);
  free(postpY);
  free(logpY);
  lst_free(tmplst);

  return total_logl;
}

/* draw the parameters of a multinomial distribution from a Dirichlet
   distrib. */
void mtf_draw_multinomial(Vector *v, double *alpha) {
  double theta[v->size];
  int i;
  dirichlet_draw(v->size, alpha, theta);
  for (i = 0; i < v->size; i++) vec_set(v, i, theta[i]);
}

/* scan a set of DNA sequences for the most prevalent n-tuples of bases */
void mtf_get_common_ntuples(SeqSet *s, List *tuples, int tuple_size, 
                            int number) {
  int alph_size = strlen(s->set->alphabet);
  int idx, cutoff, i, j, k;
  int ntuples = int_pow(alph_size, tuple_size), high = ntuples/alph_size;
  int *counts = smalloc(ntuples * sizeof(int));
  char str[tuple_size+1];
  List *tmplst, *aux_tuples;

  str[tuple_size] = '\0';

  for (i = 0; i < ntuples; i++) counts[i] = 0;
 
  for (i = 0; i < s->set->nseqs; i++) {
    idx = 0;
    for (k = 0; k < tuple_size-1; k++)
      idx += s->set->inv_alphabet[(int)s->set->seqs[i][tuple_size-k-2]] * 
        int_pow(alph_size, k);
    
    for (j = tuple_size-1; j < s->lens[i]; j++) {
      idx = (idx % high) * alph_size; /* shift left */
      idx += s->set->inv_alphabet[(int)s->set->seqs[i][j]]; /* add next base */
      counts[idx]++;
    }
  }

  tmplst = lst_new_int(ntuples);
  for (i = 0; i < ntuples; i++) 
    if (counts[i] > 0) lst_push_int(tmplst, counts[i]);
  lst_qsort_int(tmplst, DESCENDING);
  cutoff = lst_get_int(tmplst, number-1);
  lst_clear(tuples);
  aux_tuples = lst_new_ptr(number);
  for (i = 0; i < ntuples; i++) {
    if (counts[i] >= cutoff) {
      get_tuple_str(str, i, tuple_size, s->set->alphabet);
      if (counts[i] == cutoff) 
        lst_push_ptr(aux_tuples, str_new_charstr(str));
      else 
        lst_push_ptr(tuples, str_new_charstr(str));
    }
  }

  /* this is necessary to get exactly 'number' tuples */
  for (i = 0; i < lst_size(aux_tuples); i++) {
    if (lst_size(tuples) < number) 
      lst_push_ptr(tuples, lst_get_ptr(aux_tuples, i));
    else str_free(lst_get_ptr(aux_tuples, i));
  }

  lst_free(tmplst);
  lst_free(aux_tuples);
  free(counts);
}

/* randomly sample n-tuples from sequence set */
void mtf_sample_ntuples(SeqSet *s, List *tuples, int tuple_size, int number) {
  int n, i, j;
  char tuple[tuple_size + 1];
  tuple[tuple_size] = '\0';
  /* for simplicity, we'll assume they're all about the same length,
     and sample number/seqset->nseqs from each sequence */
  n = ceil(1.0 * number/s->set->nseqs);
  for (i = 0; i < s->set->nseqs && lst_size(tuples) < number; i++) {
    for (j = 0; j < n && lst_size(tuples) < number; j++) {
      int start = rint(1.0 * (s->lens[i] - tuple_size) * unif_rand());
      strncpy(tuple, &s->set->seqs[i][start], tuple_size);
      lst_push_ptr(tuples, str_new_charstr(tuple));
    }
  }
}

/* perform a "soft" initialization of a series of multinomial models
   from a consensus sequence.  Initialization can be deterministic
   (treat pseudocounts as counts) or probabilistic (draw parameters
   from Dirichlet distribution).  If target size is larger than
   consensus size, flanking models will be added with uniform
   distributions (or draws from uniform Dirichlet) */
void mtf_init_from_consensus(String *consensus, Vector **mods, 
                             int *inv_alph, int npseudocounts, 
                             int probabilistic, int target_size) {
  int i, j, sum, offset, size/* , Nfactor */;

  size = mods[0]->size;         /* mods[0] assumed backgd, not altered */

/*   Nfactor = (inv_alph[(int)'N'] >= 0 ? 1 : 0); */

  for (i = 0; i < (target_size - consensus->length) / 2; i++) {
    for (j = 0; j < size; j++) {
/*       if (j == inv_alph[(int)'N']) vec_set(mods[i+1], j, 0); */
/*       else */ 
      vec_set(mods[i+1], j, 1.0/(size/* -Nfactor */));
    }
  }
  offset = i;

  sum = size /* - Nfactor */ - 1 + npseudocounts;
  for (i = 0; i < consensus->length; i++) {
    for (j = 0; j < size; j++) {
      if (j == inv_alph[(int)consensus->chars[i]])
        vec_set(mods[i+offset+1], j, 1.0 * npseudocounts / sum);
/*       else if (j == inv_alph[(int)'N']) */
/*         vec_set(mods[i+offset+1], j, 0); */
      else
        vec_set(mods[i+offset+1], j, 1.0/sum);
    }
  }

  for (; i+offset < target_size; i++) {
    for (j = 0; j < size; j++) {
/*       if (j == inv_alph[(int)'N']) vec_set(mods[i+offset+1], j, 0); */
/*       else  */
      vec_set(mods[i+offset+1], j, 1.0/(size/* -Nfactor */));
    }
  }
}

/* given a possibly large set of consensus sequences representing
   starting models, return a subset that look promising, based on
   unmaximized scores.  (This is similar to the heuristic used
   by MEME.) */
void mtf_winnow_starts(void *data, List *origseqs, int ntochoose, 
                       List *bestseqs, int multiseq, int motif_size, 
                       TreeNode *tree, void *backgd, double *has_motif) {

  SeqSet *ss = !multiseq ? data : NULL;
  PooledMSA *pmsa = multiseq ? data : NULL;
  Vector **freqs = smalloc((motif_size + 1) * sizeof(void*));
  int alph_size = multiseq ? strlen(pmsa->pooled_msa->alphabet) : 
    strlen(ss->set->alphabet);
  int *inv_alphabet = multiseq ? pmsa->pooled_msa->inv_alphabet :
    ss->set->inv_alphabet;
  int nsamples = multiseq ? lst_size(pmsa->source_msas) : ss->set->nseqs;
  int *sample_lens = multiseq ? pmsa->lens : ss->lens;
  List *score_lst = lst_new_dbl(lst_size(origseqs));
  double *score = smalloc(lst_size(origseqs) * sizeof(double));
  List *tmplst;
  int i, j, k, s, maxlen;
  double **emissions;
  double threshold;

  for (i = 0; i <= motif_size; i++) freqs[i] = vec_new(alph_size);
          
  if (has_motif == NULL) {      /* non-discriminative case only */
    if (multiseq) 
      vec_copy(freqs[0], ((TreeModel*)backgd)->backgd_freqs);
    else 
      vec_copy(freqs[0], backgd);
  }

  maxlen = 0;
  for (s = 0; s < nsamples; s++)
    if (sample_lens[s] > maxlen) maxlen = sample_lens[s];
  emissions = (double**)smalloc((motif_size+1) * sizeof(double*));
  if (has_motif != NULL) emissions[0] = NULL;
  else emissions[0] = (double*)smalloc(maxlen * sizeof(double));
  for (j = 1; j <= motif_size; j++) 
    emissions[j] = (double*)smalloc(maxlen * sizeof(double));
  tmplst = lst_new_dbl(maxlen);

  for (i = 0; i < lst_size(origseqs); i++) {
    String *seq = lst_get_ptr(origseqs, i);
    mtf_init_from_consensus(seq, freqs, inv_alphabet, 7, 0, motif_size);

    /* compute score */
    score[i] = 0;
    for (s = 0; s < nsamples; s++) {
      Motif *m;
      double lsum;

      if (multiseq) { 
        m = mtf_new(motif_size, 1, freqs, pmsa, backgd, 0.5);
        phy_compute_emissions(emissions, (void**)m->ph_mods, motif_size+1, 
                              data, s, sample_lens[s]);
      }
      else {
        m = mtf_new(motif_size, 0, freqs, ss, NULL, 0);
        mn_compute_emissions(emissions, (void**)m->freqs, motif_size+1, 
                             data, s, sample_lens[s]);
      }

      lst_clear(tmplst);
      for (j = 0; j < sample_lens[s] - m->motif_size; j++) { /* FIXME: +1 */
        double val = 0;
        for (k = 0; k < m->motif_size; k++) {
          val += emissions[k+1][j+k];
          if (has_motif == NULL) val -= emissions[0][j+k];
        }
        lst_push_dbl(tmplst, val);
      }

      lsum = -log(sample_lens[s] - m->motif_size) + log_sum_e(tmplst);
      if (has_motif != NULL) lsum += 2 * motif_size; 
                                /* approx of threshold */

      if (has_motif == NULL) 
        score[i] += lsum;
      else if (has_motif[s] == 0)
        score[i] += -lsum - log(1 + exp(-lsum));
      else if (has_motif[s] == 1)
        score[i] += -log(1 + exp(-lsum));
      else 
        score[i] += log(m->has_motif[s]) - log(1 + exp(-lsum)) + 
          log(1 + (1-m->has_motif[s])/m->has_motif[s] * exp(-lsum));

      mtf_free(m);
    }

    lst_push_dbl(score_lst, score[i]);
  }
  lst_qsort_dbl(score_lst, DESCENDING);
  threshold = lst_get_dbl(score_lst, ntochoose-1);
  lst_clear(bestseqs);
  for (i = 0; i < lst_size(origseqs); i++) 
    if (score[i] >= threshold) {
      lst_push_ptr(bestseqs, lst_get_ptr(origseqs, i));
/*       fprintf(stderr, "Keeping %s (%.3f) \n", ((String*)lst_get_ptr(origseqs, i))->chars, score[i]); */
    }
/*     else */
/*       fprintf(stderr, "Discarding %s (%.3f) \n", ((String*)lst_get_ptr(origseqs, i))->chars, score[i]); */

  free(score);
  lst_free(score_lst);
  lst_free(tmplst);
  for (j = 0; j <= motif_size; j++) {
    if (emissions[j] != NULL) free(emissions[j]);
    vec_free(freqs[j]);
  }
  free(emissions);
  free(freqs);
}

/* create a new seqset */
SeqSet *mn_new_seqset(int nseqs) {
  SeqSet * retval = smalloc(sizeof(SeqSet));
  retval->set = msa_new(smalloc(nseqs * sizeof(void*)), 
                        smalloc(nseqs * sizeof(void*)), 
                        nseqs, -1, NULL);
  retval->lens = smalloc(nseqs * sizeof(int));
  return retval;
}

/* free a seqset */
void mn_free_seqset(SeqSet *s) {
  msa_free(s->set);
  free(s->lens);
  free(s);
}

/* create a set of individual (gapless) sequences from a set of
   multiple alignments.  The new set can either include only the
   reference sequence from each alignment (set 'refseq' to desired
   one-based index) or all sequences (set 'refseq' to -1).  The return
   value is a new SeqSet object. */
SeqSet *mtf_get_seqset(List *msas, int refseq, int min_allowable_size) {
  int i, j, k, nseqs = 0;
  SeqSet *retval;

  if (refseq < 0)              /* scan for total number */
    for (i = 0; i < lst_size(msas); i++) 
      nseqs += ((MSA*)lst_get_ptr(msas, i))->nseqs;
  else 
    nseqs = lst_size(msas);

  retval = mn_new_seqset(nseqs);

  for (i = 0, k = 0; i < lst_size(msas); i++) {
    MSA *msa = lst_get_ptr(msas, i);
    j = refseq >= 1 ? refseq-1 : 0; /* if refseq, only loop once;
                                       otherwise, go through all
                                       seqs  */
    for (; j < msa->nseqs; j++) {
      int l, m = 0;
      /* copy seq w/o gaps */
      retval->set->seqs[k] = smalloc((msa->length+1) * sizeof(char));
      for (l = 0; l < msa->length; l++)
        if (msa->seqs[j][l] != GAP_CHAR)
          retval->set->seqs[k][m++] = msa->seqs[j][l];

      if (m < min_allowable_size) {
        fprintf(stderr, "WARNING: ignoring sequence '%s'.\n", msa->names[j]);
        free(retval->set->seqs[k]);
        continue;
      }

      retval->set->seqs[k][m] = '\0';
      retval->lens[k] = m;
      if (m > retval->set->length) retval->set->length = m;

      /* also copy name */
      retval->set->names[k] = strdup(msa->names[j]);
      k++;

      if (refseq >= 0) break;
    }
  }
  if (k < nseqs) retval->set->nseqs = k;
  return retval;
}

/* this hack avoids problems with non-invertible matrices arising from
   simple initialization strategies; it's used in mtf_new, below */
/* void fuzz(Vector *v) { */
/*   double sum = 0; */
/*   int i; */
/*   srand(time(NULL)); */
/*   for (i = 0; i < v->size; i++) { */
/*     add a small amount of noise to each value */ 
/*     double val = vec_get(v, i); */
/*     val += (rand() / (1000.0 * RAND_MAX));  */
/*     vec_set(v, i, val); */
/*     sum += val; */
/*   } */
/*   vec_scale(v, 1.0/sum); */
/* } */

/* create a new Motif object.  The last two parameters are used only
   if multiseq == 1.  New copies are created of 'freqs' and (if
   appropriate) backgd_phmod. */
Motif* mtf_new(int motif_size, int multiseq, Vector **freqs, 
               void *training_data, TreeModel *backgd_phmod, 
               double scale_factor) {
  int i, n;
  Motif *m = smalloc(sizeof(Motif));
  m->motif_size = motif_size;
  m->multiseq = multiseq;

  m->freqs = smalloc((motif_size+1) * sizeof(void*));
  m->training_data = training_data;
  
  if (multiseq) {
    m->alphabet = ((PooledMSA*)training_data)->pooled_msa->alphabet;
    n = lst_size(((PooledMSA*)training_data)->source_msas);
    m->refseq = 1;

    if (backgd_phmod != NULL) {
      m->ph_mods = smalloc((motif_size+1) * sizeof(void*));
      m->ph_mods[0] = tm_create_copy(backgd_phmod);
      for (i = 1; i <= motif_size; i++) {
        m->ph_mods[i] = tm_create_copy(backgd_phmod);
        tm_free_rmp(m->ph_mods[i]); /* FIXME: better way? */
        m->ph_mods[i]->estimate_branchlens = TM_BRANCHLENS_NONE;
        m->ph_mods[i]->estimate_backgd = 1;
        m->ph_mods[i]->allow_gaps = 0;
        tm_init_rmp(m->ph_mods[i]);
        vec_copy(m->ph_mods[i]->backgd_freqs, freqs[i]);
        m->freqs[i] = m->ph_mods[i]->backgd_freqs;
                                /* in this case, point to the same
                                   object (otherwise problem keeping
                                   them synched) */
/*         fuzz(m->ph_mods[i]->backgd_freqs);  */
        tm_scale_branchlens(m->ph_mods[i], scale_factor, 0);
      }
    }
  }
  else {
    for (i = 0; i <= motif_size; i++) {
      m->freqs[i] = vec_new(freqs[i]->size);
      vec_copy(m->freqs[i], freqs[i]);
    }

    m->alphabet = ((SeqSet*)training_data)->set->alphabet;
    n = ((SeqSet*)training_data)->set->nseqs;
    m->refseq = -1;
  }

  m->alph_size = strlen(m->alphabet);
  m->postprob = smalloc(n * sizeof(double));
  m->bestposition = smalloc(n * sizeof(int));
  m->samplescore = smalloc(n * sizeof(double));
  m->score = NEGINFTY;
  m->has_motif = NULL;
  m->coord_maps = NULL;
  return m;  
}

/* free a Motif object */
void mtf_free(Motif *m) {
  int i;
  if (m->multiseq && m->ph_mods != NULL) {
    for (i = 0; i <= m->motif_size; i++) 
      tm_free(m->ph_mods[i]);
    free(m->ph_mods);
  }
  else 
    for (i = 0; i <= m->motif_size; i++) 
      vec_free(m->freqs[i]);
  if (m->coord_maps != NULL) {
    int nobs = lst_size(((PooledMSA*)m->training_data)->source_msas);
    for (i  = 0; i < nobs; i++)
      msa_map_free(m->coord_maps[i]);
    free(m->coord_maps);
  }
  free(m->freqs);
  free(m->postprob);
  free(m->bestposition);
  free(m->samplescore);
  free(m);
}

/* derive a consensus sequence from a motif (majority each position)
   Assumes DNA rather than amino-acid alphabet. */
void mtf_get_consensus(Motif *m, char *consensus) {
  int i, j;
  for (i = 1; i <= m->motif_size; i++) {
    double pur_freq = 0, pyr_freq = 0;
    consensus[i-1] = 'N';
    for (j = 0; j < m->alph_size; j++) {
      double f = vec_get(m->freqs[i], j);
      if (f > 0.5) {
        consensus[i-1] = m->alphabet[j];
        break;
      }
      else if (IS_PURINE(m->alphabet[j]))
        pur_freq += f;
      else 
        pyr_freq += f;
    }
    if (consensus[i-1] == 'N') {
      if (pur_freq > 0.75) consensus[i-1] = 'R';
      else if (pyr_freq > 0.75) consensus[i-1] = 'Y';
    }
  }  
}

/* build coordinate maps for the multiple alignments in the training
   set for a multi-sequence motif; allows efficient transformation to
   coordinate system of reference sequence */
void mtf_build_coord_maps(Motif *m) {
  PooledMSA *pmsa = m->training_data;
  int i, nobs = lst_size(pmsa->source_msas);
  m->coord_maps = smalloc(nobs * sizeof(void*));
  for (i = 0; i < nobs; i++)
    m->coord_maps[i] = msa_build_coord_map(lst_get_ptr(pmsa->source_msas, i),
                                           m->refseq);
}

/* used in various functions below for mapping alignment coordinates
   to coords in reference sequence */
static PHAST_INLINE 
int safe_map(msa_coord_map *map, int pos) {
  if (pos <= 0) return pos;
  if (pos > map->msa_len) pos = map->msa_len;
  return msa_map_msa_to_seq(map, pos);
}

/* print a human-readable description of a motif and instances of it
   found in a training data set.  Can be used for single-sequence or
   multi-sequence motifs.  Assumes there is a reference
   sequence in the case of multi-sequence motifs. */
void mtf_print(FILE *F, Motif *m) {
  int i, j, k, nobs;
  List *notfound;
  PooledMSA *pmsa = m->multiseq ? (PooledMSA*)m->training_data : NULL;
  SeqSet *ss = m->multiseq ? NULL : (SeqSet*)m->training_data;

  /* consensus */
  char *cons = smalloc((m->motif_size+1) * sizeof(char));
  cons[m->motif_size] = '\0';
  mtf_get_consensus(m, cons);
  fprintf(F, "Consensus: %s\nscore: %.3f\n\n", cons, m->score);
  free(cons);

  /* position-specific estimates */
  fprintf(F, "Motif summary:\n\n%6s", "pos.");
  for (i = 0; i < m->alph_size; i++) fprintf(F, " %6c", m->alphabet[i]);
  if (m->multiseq) fprintf(F, " %6c", 't');
  fprintf(F, "\n");
  for (i = 1; i <= m->motif_size; i++) {
    fprintf(F, "%6d", i);
    for (j = 0; j < m->alph_size; j++) 
      fprintf(F, " %6.4f", vec_get(m->freqs[i], j));
    if (m->multiseq && m->ph_mods[i] != NULL && m->ph_mods[0] != NULL) 
      fprintf(F, " %6.4f", 
              tr_total_len(m->ph_mods[i]->tree) * m->ph_mods[i]->scale / 
              tr_total_len(m->ph_mods[0]->tree));
    fprintf(F, "\n");
  }
  
  /* show best motif instances */
  nobs = m->multiseq ? lst_size(pmsa->source_msas) :
    ss->set->nseqs;
  if (m->multiseq && m->coord_maps == NULL)
    mtf_build_coord_maps(m);
  notfound = lst_new_ptr(nobs);
  if (m->multiseq)
    fprintf(F, "\nMotif instances in reference sequence:\n\n");
  else 
    fprintf(F, "\nMotif instances:\n\n");
  /* first show only instances in reference sequence */
  for (i = 0; i < nobs; i++) {
    char *name, *seq;
    int len;
    if (m->multiseq) {
      MSA *smsa = lst_get_ptr(pmsa->source_msas, i);
      name = smsa->names[m->refseq-1];
      seq = smsa->seqs[m->refseq-1];
      len = pmsa->lens[i];             
    }
    else {
      name = ((SeqSet*)m->training_data)->set->names[i];
      seq = ((SeqSet*)m->training_data)->set->seqs[i];
      len = ((SeqSet*)m->training_data)->lens[i];
    }    
    if (m->has_motif != NULL && !m->has_motif[i])
        continue;
    else if (m->has_motif == NULL && m->postprob[i] < 0.5) {
      lst_push_ptr(notfound, name);
      continue;
    }
    fprintf(F, "%-15s %6d ", name, m->multiseq ? 
            safe_map(m->coord_maps[i], m->bestposition[i] - 10 + 1) :
            m->bestposition[i] - 10 + 1);
    for (j = m->bestposition[i] - 10; 
         j < m->bestposition[i] + m->motif_size + 10; j++) {
      if (j < 0 || j >= len) printf("*");
      else if (j < m->bestposition[i] || 
               j >= m->bestposition[i] + m->motif_size) 
        fprintf(F, "%c", tolower(seq[j]));
      else 
        fprintf(F, "%c", seq[j]);
    }
    fprintf(F, "%6d (%.2f)\n", m->multiseq ?
            safe_map(m->coord_maps[i], m->bestposition[i] + m->motif_size + 10 - 1) :
            m->bestposition[i] + m->motif_size + 10 - 1, 
            m->samplescore[i]);
  }
  if (lst_size(notfound) > 0) {
    fprintf(F, "\nInstances were not found in the following sequences (posterior prob. < 50%%):\n");
    for (i = 0; i < lst_size(notfound); i++)
      fprintf(F, "\t%s\n", (char*)lst_get_ptr(notfound, i));
  }

  /* now show complete alignments (if m->multiseq) */
  if (m->multiseq) {
    fprintf(F, "\nFull alignments:\n\n");
    for (i = 0; i < lst_size(pmsa->source_msas); i++) {
      MSA *msa = lst_get_ptr(pmsa->source_msas, i);
      if ((m->has_motif != NULL && !m->has_motif[i]) ||
          (m->has_motif == NULL && m->postprob[i] < 0.5)) 
      continue;
      for (k = 0; k < msa->nseqs; k++) {
        if (k == m->refseq-1)
          fprintf(F, "%-15s %6d ", msa->names[k], m->multiseq ?
                  safe_map(m->coord_maps[i], m->bestposition[i] - 10 + 1) :
                  m->bestposition[i] - 10 + 1);
        else 
          fprintf(F, "%-15s %6s ", msa->names[k], ""); /* space */
        for (j = m->bestposition[i] - 10; 
             j < m->bestposition[i] + m->motif_size + 10; j++) {
          if (j < 0 || j >= msa->length) fprintf(F, "*");
          else if (j < m->bestposition[i] || 
                   j >= m->bestposition[i] + m->motif_size) 
            fprintf(F, "%c", tolower(msa->seqs[k][j]));
          else 
            fprintf(F, "%c", msa->seqs[k][j]);
        }
        if (k == 0) 
          fprintf(F, "%6d (%.2f)", m->multiseq ?
                  safe_map(m->coord_maps[i], m->bestposition[i] + m->motif_size + 10 - 1) :
                  m->bestposition[i] + m->motif_size + 10 - 1, 
                  m->samplescore[i]);
        fprintf(F, "\n");
      }
      fprintf(F, "\n");
    }
  }
  lst_free(notfound);
}

/* print an html-formatted version of a motif and instances of it
   found in a training data set.  Can be used for single-sequence or
   multi-sequence motifs.  Assumes there is a reference
   sequence in the case of multi-sequence motifs. */
void mtf_print_html(FILE *F, Motif *m) {
  int i, j, k, nobs;
  List *notfound;
  PooledMSA *pmsa = m->multiseq ? (PooledMSA*)m->training_data : NULL;
  SeqSet *ss = m->multiseq ? NULL : (SeqSet*)m->training_data;
  int *inv_alphabet = m->multiseq ? pmsa->pooled_msa->inv_alphabet :
    ss->set->inv_alphabet;
  int *motstart, *motend;
  char **color = smalloc(m->alph_size * sizeof(char*));
  GFF_Feature **posfeat;

  /* consensus */
  char *cons = smalloc((m->motif_size+1) * sizeof(char));
  cons[m->motif_size] = '\0';
  mtf_get_consensus(m, cons);

  /* temporary */
  color[0] = "blue"; color[1] = "red"; color[2] = "green"; color[3] = "purple";

  fprintf(F, "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n\
<html>\n\
<head><title>Motif details</title></head>\n\
<body>\n\
<h1>Motif details</h1>\n");

  fprintf(F, "Consensus: %s<br>\nscore: %.3f<br>\n\n", cons, m->score);
  free(cons);

  /* position-specific estimates */
  fprintf(F, "<br>Position-specific distributions:<br><br>\n\n");
  fprintf(F, "<table border=\"1\" rules=\"all\">\n");
  fprintf(F, "<tr><td align=\"center\">pos.</td>");
  for (i = 0; i < m->alph_size; i++) fprintf(F, "<td align=\"center\">%c</td>", m->alphabet[i]);
  if (m->multiseq) fprintf(F, "<td align=\"center\">%c</td>", 't');
  fprintf(F, "</tr>\n");
  for (i = 1; i <= m->motif_size; i++) {
    fprintf(F, "<tr><td align=\"center\">%d</td>", i);
    for (j = 0; j < m->alph_size; j++) 
      fprintf(F, "<td align=\"center\">%.4f</td>", vec_get(m->freqs[i], j));
    if (m->multiseq && m->ph_mods[i] != NULL && m->ph_mods[0] != NULL) 
      fprintf(F, "<td align=\"center\">%.4f</td>", 
              tr_total_len(m->ph_mods[i]->tree) * m->ph_mods[i]->scale / 
              tr_total_len(m->ph_mods[0]->tree));
    fprintf(F, "</tr>\n");
  }
  fprintf(F, "</table>\n");
  
  /* show best motif instances */
  nobs = m->multiseq ? lst_size(pmsa->source_msas) :
    ss->set->nseqs;
  if (m->multiseq && m->coord_maps == NULL)
    mtf_build_coord_maps(m);
  posfeat = smalloc(nobs * sizeof(void*));
  for (i = 0; i < nobs; i++) posfeat[i] = NULL;
  motstart = smalloc(nobs * sizeof(int));
  motend = smalloc(nobs * sizeof(int));

  notfound = lst_new_ptr(nobs);
  if (m->multiseq)
    fprintf(F, "\n<br><br>Best instances found (reference sequence version):<br><br>\n\n");
  else 
    fprintf(F, "\n<br><br>Best instances found:<br><br>\n\n");
  /* first show only instances in reference sequence */
  fprintf(F, "<table>\n");
  for (i = 0; i < nobs; i++) {
    char *name, *seq;
    int len;
    if (m->multiseq) {
      name = ((MSA*)lst_get_ptr(pmsa->source_msas, i))->names[m->refseq-1];
      seq = ((MSA*)lst_get_ptr(pmsa->source_msas, i))->seqs[m->refseq-1];
      len = pmsa->lens[i];             
    }
    else {
      name = ((SeqSet*)m->training_data)->set->names[i];
      seq = ((SeqSet*)m->training_data)->set->seqs[i];
      len = ((SeqSet*)m->training_data)->lens[i];
    }    
    if (m->has_motif != NULL && !m->has_motif[i]) 
        continue;
    else if (m->has_motif == NULL && m->postprob[i] < 0.5) {
      lst_push_ptr(notfound, name);
      continue;
    }

    /* get genomic coords of the motif and its surrounding region, if
       possible */
    posfeat[i] = 
      gff_new_feature_genomic_pos(str_new_charstr(name), str_new(1), 
                                  str_new(1), 0, GFF_NULL_FRAME, 
                                  str_new(1), 1);
                                /* returns NULL if name is not in the
                                   format of genomic coords */
    if (posfeat[i] != NULL) {
      if (posfeat[i]->strand == '-') { 
        motend[i] = posfeat[i]->end - (m->multiseq ? 
                 safe_map(m->coord_maps[i], m->bestposition[i] + 1) - 1:
                 m->bestposition[i]);
        motstart[i] = motend[i] - m->motif_size + 1;
      }
      else {
        motstart[i] = posfeat[i]->start + (m->multiseq ?
                   safe_map(m->coord_maps[i], m->bestposition[i] + 1) - 1:
                   m->bestposition[i]);
        motend[i] = motstart[i] + m->motif_size - 1;
      }
      motstart[i] -= m->motif_size / 2;
      motend[i] += m->motif_size / 2;
    }

    /* now print row in table */
    fprintf(F, "<tr><td>");
    if (posfeat[i] != NULL)
      fprintf(F, "<a href=\"%s&position=%s:%d-%d\" TARGET=_blank>", 
              HGTRACKS_URL, posfeat[i]->seqname->chars, 
              posfeat[i]->start, posfeat[i]->end);
    fprintf(F, "%s", name);
    if (posfeat[i] != NULL) fprintf(F, "</a>");
    fprintf(F, "</td><td align=\"right\">%d</td><td><font face=\"Courier\">", 
            m->multiseq ? 
            safe_map(m->coord_maps[i], m->bestposition[i] - 10 + 1) :
            m->bestposition[i] - 10 + 1);
    for (j = m->bestposition[i] - 10; 
         j < m->bestposition[i] + m->motif_size + 10; j++) {
      if (j < 0 || j >= len) fprintf(F, "*");
      else if (j < m->bestposition[i] || 
               j >= m->bestposition[i] + m->motif_size) 
        fprintf(F, "%c", tolower(seq[j]));
      else {
        if (j == m->bestposition[i] && posfeat[i] != NULL)
          fprintf(F, "<a href=\"%s&position=%s:%d-%d\" TARGET=_blank>", 
                  HGTRACKS_URL, posfeat[i]->seqname->chars, motstart[i], 
                  motend[i]);
        fprintf(F, "<font color=\"%s\">%c</font>", 
                inv_alphabet[(int)seq[j]] >= 0 ? 
                color[inv_alphabet[(int)seq[j]]] : "gray", seq[j]);
        if (j == m->bestposition[i] + m->motif_size - 1 && posfeat[i] != NULL)
          fprintf(F, "</a>");
      }
    }
    fprintf(F, "</font></td><td align=\"right\">%d</td><td align=\"right\">(%.2f)</td></tr>\n", 
            m->multiseq ?
            safe_map(m->coord_maps[i], m->bestposition[i] + m->motif_size + 10 - 1) :
            m->bestposition[i] + m->motif_size + 10 - 1, m->samplescore[i]);
  }
  fprintf(F, "</table>\n");

  if (lst_size(notfound) > 0) {
    fprintf(F, "\n<br>Instances were not found in the following sequences (posterior prob. < 50%%):\n<ul>");
    for (i = 0; i < lst_size(notfound); i++)
      fprintf(F, "<li>%s\n", (char*)lst_get_ptr(notfound, i));
    fprintf(F, "</ul>\n");
  }

  /* now show complete alignments (if m->multiseq) */
  if (m->multiseq) {
    fprintf(F, "\n<br><br>Full alignments:<br><br>\n\n");
    fprintf(F, "<table>\n");
    for (i = 0; i < lst_size(pmsa->source_msas); i++) {
      MSA *msa = lst_get_ptr(pmsa->source_msas, i);
      if ((m->has_motif != NULL && !m->has_motif[i]) ||
          (m->has_motif == NULL && m->postprob[i] < 0.5)) 
      continue;
      for (k = 0; k < msa->nseqs; k++) {
        if (k == m->refseq-1) {
          fprintf(F, "<tr><td>");
          if (posfeat[i] != NULL)
            fprintf(F, "<a href=\"%s&position=%s:%d-%d\" TARGET=_blank>", 
                    HGTRACKS_URL, posfeat[i]->seqname->chars, 
                    posfeat[i]->start, posfeat[i]->end);
          fprintf(F, "%s", msa->names[k]);
          if (posfeat != NULL) fprintf(F, "</a>");          
          fprintf(F, "</td><td align=\"right\">%d</td>", 
                  m->multiseq ?
                  safe_map(m->coord_maps[i], m->bestposition[i] - 10 + 1) :
                  m->bestposition[i] - 10 + 1);
        }
        else 
          fprintf(F, "<tr><td>%s</td><td></td>", msa->names[k]);
        fprintf(F, "<td><font face=\"Courier\">");
        for (j = m->bestposition[i] - 10; 
             j < m->bestposition[i] + m->motif_size + 10; j++) {
          if (j < 0 || j >= msa->length) fprintf(F, "*");
          else if (j < m->bestposition[i] || 
                   j >= m->bestposition[i] + m->motif_size) 
            fprintf(F, "%c", tolower(msa->seqs[k][j]));
          else {
            if (j == m->bestposition[i] && posfeat[i] != NULL)
              fprintf(F, "<a href=\"%s&position=%s:%d-%d\" TARGET=_blank>", 
                      HGTRACKS_URL, posfeat[i]->seqname->chars, motstart[i], 
                      motend[i]);
            fprintf(F, "<font color=\"%s\">%c</font>", 
                    inv_alphabet[(int)msa->seqs[k][j]] >= 0 ?
                    color[inv_alphabet[(int)msa->seqs[k][j]]] : "gray", 
                    msa->seqs[k][j]);
            if (j == m->bestposition[i] + m->motif_size - 1 && posfeat[i] != NULL)
              fprintf(F, "</a>");
          }
        }
        fprintf(F, "</font></td>");
        if (k == 0) 
          fprintf(F, "<td align=\"right\">%d</td><td align=\"right\">(%.2f)</td></tr>\n", 
                  m->multiseq ?
                  safe_map(m->coord_maps[i], m->bestposition[i] + m->motif_size + 10 - 1) :
                  m->bestposition[i] + m->motif_size + 10 - 1, 
                  m->samplescore[i]);
        else 
          fprintf(F, "<td></td><td></td></tr>\n");
      }
      fprintf(F, "\n");
    }
    fprintf(F, "</table>\n");
  }
  fprintf(F, "</body>\n</html>\n");
  lst_free(notfound);
  for (i = 0; i < nobs; i++)
    if (posfeat[i] != NULL) gff_free_feature(posfeat[i]);
  free(posfeat); free(motstart); free(motend);
}

void mtf_print_summary_html(FILE *F, List *motifs, String *prefix) {

  char cons[STR_SHORT_LEN];
  int i;
  String *link_fname = str_new(STR_SHORT_LEN);

  fprintf(F, "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n\
<html><head><title>Best-scoring motifs</title></head><body><center><p>\n\
<table border=\"1\" rules=\"all\">\n\
<caption><em>Best-scoring motifs</em></caption>\n\
<tr><td align=\"center\">No.</td>\n\
<td align=\"center\">Consensus</td>\n\
<td align=\"center\">Score</td></tr>\n");

  for (i = 0; i < lst_size(motifs); i++) {
    Motif *m = lst_get_ptr(motifs, i);
    str_cpy(link_fname, prefix);
    str_remove_path(link_fname);
    str_append_int(link_fname, i+1);
    str_append_charstr(link_fname, ".html");

    cons[m->motif_size] = '\0';
    mtf_get_consensus(m, cons);
    fprintf(F, "<td align=\"center\">%d</td><td align=\"center\"><a href=\"%s\"><font face=\"Courier\">%s</font></a></td><td align=\"right\">%.3f</td></tr>\n", i+1, link_fname->chars, cons, m->score);
  }

  fprintf(F, "</table></body></html>\n");
  str_free(link_fname);
}

/* Predict the best motif in each sample of a data set.  If
   m->multiseq, data should be a PooledMSA, otherwise it should be a
   SeqSet.  If has_motif != NULL, then predictions will be made only
   for samples i such that has_motif[i] >= 0.5.  The starting position
   of the best motifs for each sample i is stored in bestposition[i]
   and a score is stored in score[i]. */
void mtf_predict(Motif *m, void *data, int *bestposition, double *score, 
                 double *has_motif) {
  PooledMSA *pmsa = m->multiseq ? m->training_data : NULL;
  SeqSet *ss = m->multiseq ? NULL : m->training_data;
  int nsamples = m->multiseq ? lst_size(pmsa->source_msas) : ss->set->nseqs;
  int *sample_lens = m->multiseq ? pmsa->lens : ss->lens;
  int i, j, s, maxlen = 0;
  double **emissions;

  for (s = 0; s < nsamples; s++)
    if (sample_lens[s] > maxlen) maxlen = sample_lens[s];

  emissions = (double**)smalloc((m->motif_size+1) * sizeof(double*));
  emissions[0] = NULL;
  for (i = 1; i <= m->motif_size; i++)
    emissions[i] = (double*)smalloc(maxlen * sizeof(double));
                          
  for (s = 0; s < nsamples; s++) {
    if (has_motif != NULL && has_motif[s] < 0.5) { 
      bestposition[s] = -1;
      score[s] = 0; 
      continue; 
    }

    if (m->multiseq) 
      phy_compute_emissions(emissions, (void*)m->ph_mods, m->motif_size+1, 
                            pmsa, s, sample_lens[s]);
    else 
      mn_compute_emissions(emissions, (void*)m->freqs, m->motif_size+1, ss, 
                           s, sample_lens[s]);

    score[s] = NEGINFTY;      
    for (i = 0; i < sample_lens[s] - m->motif_size; i++) {
      double thisscore = 0;
      for (j = 0; j < m->motif_size; j++) 
        thisscore += emissions[j+1][i+j]; /* check: can ignore
                                             threshold/prior? */
      if (thisscore >= score[s]) {
        bestposition[s] = i;
        score[s] = thisscore;
      }
    }
  }

  for (i = 1; i <= m->motif_size; i++)
    free(emissions[i]);
  free(emissions);
}

/* add a feature to a gff for each predicted instance of a motif in
   the training set */
void mtf_add_features(Motif *m, GFF_Set *gff) {
  PooledMSA *pmsa = m->multiseq ? (PooledMSA*)m->training_data : NULL;
  SeqSet *ss = m->multiseq ? NULL : (SeqSet*)m->training_data;
  int i, nobs = m->multiseq ? lst_size(pmsa->source_msas) : ss->set->nseqs;

  if (m->multiseq && m->coord_maps == NULL)
    mtf_build_coord_maps(m);

  for (i = 0; i < nobs; i++) {
    GFF_Feature *f;
    char *name = m->multiseq ? 
      ((MSA*)lst_get_ptr(pmsa->source_msas, i))->names[m->refseq-1] :
      ((SeqSet*)m->training_data)->set->names[i];

    if ((m->has_motif != NULL && !m->has_motif[i]) ||
        (m->has_motif == NULL && m->postprob[i] < 0.5)) 
      continue;

    if ((f = gff_new_feature_genomic_pos(str_new_charstr(name), 
                                         str_new_charstr("phast_motif"),
                                         str_new_charstr("predicted_motif"), 
                                         m->samplescore[i], GFF_NULL_FRAME,
                                         str_new(1), 0)) == NULL) {
      fprintf(stderr, "WARNING: cannot create feature for sequence '%s'.\n",
              name);
      continue;
    }

    if (f->strand == '-') {
      f->end -= (m->multiseq ? 
                 safe_map(m->coord_maps[i], m->bestposition[i] + 1) - 1:
                 m->bestposition[i]);
      f->start = f->end - m->motif_size + 1;
    }
    else {
      f->start += (m->multiseq ?
                   safe_map(m->coord_maps[i], m->bestposition[i] + 1) - 1:
                   m->bestposition[i]);
      f->end = f->start + m->motif_size - 1;
    }

    lst_push_ptr(gff->features, f);
  }
}


/****************************OBVIATED***************************************/
/* Create an HMM that defines a motif model like the one used by MEME,
   with a single background state and a sequence of <size> states
   representing a motif.  The parameter 'prior' is the prior
   probability that the motif begins at any particular position.  For
   simplicity, we assume the motif does not begin in the first
   position or extend to the last position of a sequence.  If
   single_motif == 1, then the motif states can be visited at most
   once in a path through the HMM.  */
/* HMM* mtf_as_hmm(int size, double prior, int single_motif) { */
/*   int i; */
/*   HMM *hmm = single_motif ? hmm_new_nstates(size+2, 1, 1) :  */
/*     hmm_new_nstates(size+1, 1, 1); */ /* need an extra backgd state if
                                      single_motif */

/*   mm_set(hmm->transition_matrix, 0, 1, prior);  */
                                /* prob of entering motif */
/*   for (i = 1; i < size; i++) mm_set(hmm->transition_matrix, i, i+1, 1); */
                                /* probs of 1 through motif states */
/*   vec_set(hmm->begin_transitions, 0, 1); */
                                /* prob 1 of starting in backgd */
/*   mm_set(hmm->transition_matrix, 0, 0, 1-prior-MTF_EPSILON); */
                                /* prob of staying in backgd */
/*   vec_set(hmm->end_transitions, 0, MTF_EPSILON); */
                                /* ensures that you end in backgd */

/*   if (single_motif == 0)  */
/*     mm_set(hmm->transition_matrix, size, 0, 1); */
                                /* prob 1 of going to backgd after
                                   motif (implicit assumption that
                                   motifs cannot be adjacent) */
/*   else { */
/*     mm_set(hmm->transition_matrix, size, size+1, 1); */
                                /* prob 1 of going to 2nd backgd after
                                   motif */
/*     mm_set(hmm->transition_matrix, size+1, size+1, 1-prior-MTF_EPSILON); */
                                /* prob of staying in 2nd backgd */
                                /* NOTE: there's implicitly prob
                                   'prior' of trans. to a dummy
                                   state, w/ const. emission prob 0
                                   (necessary to avoid a bias favoring
                                   2nd backgd state).  WARNING:
                                   problem if transitions are
                                   normalized ... */
/*     vec_set(hmm->end_transitions, size+1, MTF_EPSILON); */
                                /* allows you to end in 2nd backgd also */
/*   } */

/*   hmm_reset(hmm); */

/*   return hmm; */
/* } */


/* Find a motif from a collection of multiple alignments.  Similar to
   above, but uses multiple alignments in place of single sequences,
   and phylogenetic models in place of multinomial distributions.
   Return value is an array of motif_size+1 pointers to
   TreeModels.  */
/* TreeModel** mtf_find_msa_em(List *msas, TreeNode *tree, int motif_size,  */
/*                             double prior, int nrestarts,  */
/*                             List *init_list, int sample_parms,  */
/*                             int nbest, int npseudocounts, */
/*                             double *postprob, int *bestposition, */
/*                             double *bestlogl) { */

/*   TreeModel **tmpmodels = smalloc((motif_size+1) * sizeof(void*)), */
/*     **models = smalloc((motif_size+1) * sizeof(void*)); */
/*   int i, cons, trial, alph_size; */
/*   double logl; */
/*   double *tmppostprob = smalloc(lst_size(msas) * sizeof(double)); */
/*   int *tmpbestposition = smalloc(lst_size(msas) * sizeof(int)), */
/*     *msa_lens = smalloc(lst_size(msas) * sizeof(int)); */
/*   double alpha[strlen(((MSA*)lst_get_ptr(msas, 0))->alphabet)]; */
/*   List *tmpl = NULL; */
/*   PooledMSA *pmsa; */
/*   Vector **mult = smalloc((motif_size+1) * sizeof(void*)); */

/*   *bestlogl = NEGINFTY; */

  /* estimate background model */
/*   pmsa = ss_pooled_from_msas(msas, 1, motif_size, NULL); */

/*   msa_remove_N_from_alph(pmsa->pooled_msa); */
/*   for (i = 0; i < lst_size(pmsa->source_msas); i++) */
/*     msa_remove_N_from_alph(lst_get_ptr(pmsa->source_msas, i)); */
/*   alph_size = strlen(pmsa->pooled_msa->alphabet); */

/*   tmpmodels[0] = tm_new(tr_create_copy(tree), NULL, NULL, F81,  */
/*                         pmsa->pooled_msa->alphabet, 1, 0, -1); */
/*   tm_fit(tmpmodels[0], pmsa->pooled_msa,  */
/*          tm_params_init(tmpmodels[0], .1, 5, 0),  */
/*          -1, -1, 0, OPT_MED_PREC, NULL); */
/*   models[0] = tmpmodels[0]; */

  /* various initializations */
/*   mult[0] = NULL; */
/*   for (i = 1; i <= motif_size; i++) { */
/*     tmpmodels[i] = models[i] = NULL; */
/*     mult[i] = vec_new(alph_size); */
/*   } */
/*   for (i = 0; i < alph_size; i++) alpha[i] = 1;       */
/*   for (i = 0; i < lst_size(msas); i++) */
/*     msa_lens[i] = ((MSA*)lst_get_ptr(msas, i))->length; */

  /* select subset of init strings, if necessary */
  /* better way to do this?  use conservation? */
  /* FIXME: need to derive a seqset from the msas to use this */
/*   if (nbest > 0 && init_list != NULL) { */
/*     fprintf(stderr, "Winnowing candidate start strings ...\n"); */
/*     tmpl = lst_new_ptr(nbest); */
/*     mtf_winnow_starts(msas, msa_lens, tmpmodels[0], motif_size,  */
/*                       nbest, tmpl, init_list); */
/*     init_list = tmpl; */
/*   } */

/*   for (cons = 0;  */
/*        cons < (init_list == NULL ? 1 : lst_size(init_list));  */
/*        cons++) {  */             /* (loop only once if no init_list) */

/*     String *initstr = init_list == NULL ? NULL :  */
/*       lst_get_ptr(init_list, cons); */

/*     for (trial = 0; trial < nrestarts; trial++) { */

/*       if (nrestarts == 1) */
/*         fprintf(stderr, "Trying candidate %d ... ", cons+1); */
/*       else  */
/*         fprintf(stderr, "Trying candidate %d, trial %d ... ",  */
/*                 cons+1, trial+1); */

      /* first set up multionomial models, possibly based on
         consensus seqs; reuse routines from above */
/*       if (initstr == NULL) */
/*         for (i = 1; i <= motif_size; i++)  */
/*           mtf_draw_multinomial(mult[i], alpha); */
/*       else */
/*         mtf_init_from_consensus(initstr, mult, pmsa->pooled_msa->inv_alphabet, */
/*                                 npseudocounts, sample_parms, motif_size); */

      /* now derive tree models from them  */
/*       for (i = 1; i <= motif_size; i++) { */
        /* FIXME: avoid copying each time? avoid reinit? */
/*         if (tmpmodels[i] != NULL) tm_free(tmpmodels[i]); */ /* FIXME: init */
/*         tmpmodels[i] = tm_create_copy(tmpmodels[0]); */
/*         tm_free_rmp(tmpmodels[i]); */
/*         tmpmodels[i]->estimate_scale_only = 1; */
/*         tmpmodels[i]->estimate_backgd = 1; */
/*         tm_init_rmp(tmpmodels[i]); */
/*         vec_copy(tmpmodels[i]->backgd_freqs, mult[i]); */
/*         tm_scale_branchlens(tmpmodels[i], 0.5, 0); */
/*       } */

/*       logl = mtf_find_em(tmpmodels, pmsa, lst_size(msas), msa_lens,  */
/*                          motif_size, prior, phy_compute_emissions,  */
/*                          phy_estim_mods, phy_get_obs_idx, tmppostprob, */
/*                          tmpbestposition); */

/*       fprintf(stderr, "(logL = %.1f)\n", logl); */

/*       if (logl > *bestlogl) { */
/*         *bestlogl = logl; */
/*         for (i = 1; i <= motif_size; i++) { */
/*           if (models[i] != NULL) tm_free(models[i]);  */
/*           models[i] = tm_create_copy(tmpmodels[i]); */
/*         } */
/*         if (postprob != NULL) */
/*           for (i = 0; i < lst_size(msas); i++) */
/*             postprob[i] = tmppostprob[i]; */
/*         if (bestposition != NULL) */
/*           for (i = 0; i < lst_size(msas); i++) */
/*             bestposition[i] = tmpbestposition[i]; */
/*       } */
/*     } */
/*   } */

/*   for (i = 1; i <= motif_size; i++) { */
/*     tm_free(tmpmodels[i]); */
/*     vec_free(mult[i]); */
/*   } */
/*   free(tmpmodels); */
/*   free(mult); */
/*   free(tmpbestposition); */
/*   free(tmppostprob); */
/*   free(msa_lens); */
/*   if (tmpl != NULL) lst_free(tmpl); */
/*   ss_free_pooled_msa(pmsa); */
/*   return models; */
/* } */


