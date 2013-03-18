/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: fit_em.c,v 1.8.2.1 2009-03-19 17:19:07 mt269 Exp $ */

/* Functions for fitting tree models by EM */


#include <fit_em.h>
#include <misc.h>
#include <tree_model.h>
#include <stringsplus.h>
#include <ctype.h>
#include <numerical_opt.h>
#include <tree_likelihoods.h>
#include <subst_mods.h>
#include <time.h>
#include <sys/time.h>
#include <sufficient_stats.h>
#include <matrix.h>
#include <sys/types.h>
#include <unistd.h>
#include <complex.h>
#include <math.h>
#include <dgamma.h>

#define DERIV_EPSILON 1e-5      
                                /* used for numerical est. of derivatives */

/* internal functions */
double tm_partial_ll_wrapper(Vector *params, void *data);
double tm_partial_ll_wrapper_fast(Vector *params, void *data);
double tm_likelihood_wrapper(Vector *params, void *data);
void tm_log_em(FILE *logf, int header_only, double val, Vector *params);
void compute_grad_em_approx(Vector *grad, Vector *params, void *data, 
                            Vector *lb, Vector *ub);
void compute_grad_em_exact(Vector *grad, Vector *params, void *data, 
                           Vector *lb, Vector *ub);
void get_neighbors(int *neighbors, int state, int order, int alph_size);


void remove_ratevar_from_param_map(TreeModel *mod, Vector *params);

/* fit a tree model using EM */
int tm_fit_em(TreeModel *mod, MSA *msa, Vector *params, int cat, 
              opt_precision_type precision, int max_its, FILE *logf,
	      FILE *error_file) {
  double ll, improvement;
  Vector *lower_bounds, *upper_bounds, *opt_params;
  int retval = 0, it, i, home_stretch = 0, nratecats, npar;
  double lastll = NEGINFTY, alpha=0, rK0=0, freqK0=0, branchlen_scale;
  struct timeval start_time, end_time, post_prob_start, post_prob_end;
  TreeModel *proj_mod = NULL;
  char tmp_mod_fname[STR_SHORT_LEN];
  FILE *F;
  void (*grad_func)(Vector*, Vector*, void*, Vector*, 
                    Vector*);
  double (*likelihood_func)(Vector *, void*);
  Matrix *H;
  int opt_ratevar_freqs=0;
  opt_precision_type bfgs_prec = OPT_LOW_PREC;
                                /* will be adjusted as necessary */

  /* obtain sufficient statistics for MSA, if necessary */
  if (msa->ss == NULL) {
    if (msa->seqs == NULL)
      die("ERROR tm_fit_em: msa->seqs == NULL\n");
    ss_from_msas(msa, mod->order+1, 0, NULL, NULL, NULL, -1, 
		 subst_mod_is_codon_model(mod->subst_mod));
  }

  if (mod->backgd_freqs == NULL) {
    tm_init_backgd(mod, msa, cat);
    for (i=0; i<mod->backgd_freqs->size; i++)
      vec_set(params, mod->backgd_idx+i, vec_get(mod->backgd_freqs, i));
  }

  if (mod->tree == NULL) {      /* weight matrix */
    mod->lnL = tl_compute_log_likelihood(mod, msa, NULL, NULL, cat, NULL) * 
      log(2);
    return 0;
  }

  /* if using an expensive model, set up filename for temporary files */
  if (mod->order >= 2) {
    char tmpstr[20];

    strcpy(tmpstr, "phast");

    sprintf(tmp_mod_fname, "fit_em.%s.%d.mod", tmpstr, getpid());
  }

  if (logf != NULL)
    gettimeofday(&start_time, NULL);

  /* package with mod any data needed to compute likelihoods */
  mod->msa = msa;               
  mod->tree_posteriors = tl_new_tree_posteriors(mod, msa, 0, 0, 0, 1, 0, 0,
                                                mod->empirical_rates ? 1 : 0);
  mod->category = cat;

  /* in the case of rate variation, start by ignoring then reinstate
     when close to convergence */
  nratecats = mod->nratecats;
  if (nratecats > 1 && mod->param_map[mod->ratevar_idx] != -1) {
    remove_ratevar_from_param_map(mod, params);
    opt_ratevar_freqs = 1;
    if (mod->alt_subst_mods == NULL) {
      alpha = mod->alpha;
      mod->alpha = -nratecats;    /* code to ignore rate variation
  				   temporarily */
      mod->nratecats = 1;
      freqK0 = mod->freqK[0];
      rK0 = mod->rK[0];
      mod->freqK[0] = mod->rK[0] = 1;
    }
  }

  if (logf != NULL) tm_log_em(logf, 1, 0, params);
  grad_func = (proj_mod == NULL && 
               mod->estimate_branchlens == TM_BRANCHLENS_ALL && 
               mod->subst_mod != JC69 && mod->subst_mod != F81 &&
               !mod->estimate_backgd && mod->alt_subst_mods==NULL &&
	       mod->selection_idx < 0 ? 
	       compute_grad_em_approx : NULL);
                                /* functions for analytical gradients
                                   don't yet know about estimating
                                   scale or backgd freqs, also require
                                   diagonalization  */

  /* routines for derivative computation assume complex numbers */
  if (grad_func != NULL)
    mm_set_eigentype(mod->rate_matrix, COMPLEX_NUM);

  if (mod->estimate_backgd)
    likelihood_func = tm_likelihood_wrapper;
  else likelihood_func = tm_partial_ll_wrapper;

  npar=0;
  for (i=0; i<params->size; i++) 
    if (mod->param_map[i] >= npar) 
      npar = mod->param_map[i]+1;
  opt_params = vec_new(npar);
  for (i=0; i<npar; i++) 
    vec_set(opt_params, i, 0);  // in some cases some parameters are not used but need to be initialized
  for (i=0; i<params->size; i++) {
    if (mod->param_map[i] >= 0) 
      vec_set(opt_params, mod->param_map[i], vec_get(params, i));
    vec_set(mod->all_params, i, vec_get(params, i));
  }

  tm_new_boundaries(&lower_bounds, &upper_bounds, npar, mod, 0);

  H = mat_new(npar, npar);
  mat_set_identity(H);

  if (mod->estimate_branchlens == TM_BRANCHLENS_NONE ||
      mod->alt_subst_mods != NULL ||
      mod->selection_idx >= 0)
    mod->scale_during_opt = 1;
  if (mod->estimate_branchlens == TM_BRANCHLENS_ALL) {
    for (i=0; i < mod->tree->nnodes; i++) {
      TreeNode *n = lst_get_ptr(mod->tree->nodes, i);
      if (n != mod->tree && n->hold_constant)
	mod->scale_during_opt = 1;
    }
  }


  for (it = 1;  ; it++) {
    double tmp;
    checkInterrupt();

    tm_unpack_params(mod, opt_params, -1);
    
    /* if appropriate, dump intermediate version of model for inspection */
    if (mod->order >= 2) {
      F = phast_fopen(tmp_mod_fname, "w+"); tm_print(F, mod); phast_fclose(F);
    }

    /* obtain posterior probabilities and likelihood */
    if (logf != NULL) 
      gettimeofday(&post_prob_start, NULL);

    ll = tl_compute_log_likelihood(mod, msa, NULL, NULL, cat, mod->tree_posteriors) 
      * log(2); 

    if (logf != NULL) {
      gettimeofday(&post_prob_end, NULL);
      fprintf(logf, "\nTime to collect posterior probabilities: %.4f sec.\n", 
              post_prob_end.tv_sec - post_prob_start.tv_sec + 
              (post_prob_end.tv_usec - post_prob_start.tv_usec)/1.0e6);      
      tm_log_em(logf, 0, ll, mod->all_params);
    }

    improvement = fabs((lastll - ll)/ll);
    lastll = ll;

    /* check convergence */
    if ((max_its > 0 && it > max_its) || 
	(improvement < TM_EM_CONV(precision) &&
	 bfgs_prec == precision && mod->nratecats == nratecats))
                                /* don't exit unless precision is
                                   already at its max and rate
                                   variation has been reintroduced (if
                                   necessary) */
      break;

    /* adjust inner optimization strategy as necessary */
    if (improvement < TM_EM_CONV(OPT_CRUDE_PREC)) {
      /* change gradient function first (if necessary), and use
         slightly better precision with BFGS; then on a subsequent
         pass maximize BFGS precision.  A big jump in likelihood
         usually occurs when the gradient function is first changed. */
      if (grad_func == compute_grad_em_approx) {
        if (logf != NULL) fprintf(logf, "Switching to exact gradients.\n");
        grad_func = compute_grad_em_exact;
        if (bfgs_prec == OPT_LOW_PREC && bfgs_prec != precision) 
          bfgs_prec = OPT_MED_PREC;
      }
      else {
        home_stretch = 1;
        if (bfgs_prec != precision) {
          if (logf != NULL) 
            fprintf(logf, "Switching to higher precision with BFGS.\n");
          bfgs_prec = precision;
        }
      }
    }
    /* NOTE: medium precision case could possibly be improved
       by switching to high precision for BFGS approx when improvement <
       2e-5.  Could save some iterations in the outer loop (collecting
       post probs is expensive).  Will require some testing and tuning to
       get it right.  Probably try with U2S, R2, U2. */
    
   /* in case of empirical rate variation, also reestimate mixing
       proportions (rate weights).  The m.l.e. for these is a simple
       function of the posterior probs. of the rate categories and the
       rate constants (the maximization step of EM decomposes into two
       separate problems) */
    if (mod->nratecats > 1 && mod->empirical_rates && opt_ratevar_freqs) {
      if (mod->site_model) {
	tm_site_model_set_ml_weights(mod, params, mod->tree_posteriors->rcat_expected_nsites);
      } else {
	double sum = 0;
	for (i = 0; i < mod->nratecats; i++) 
	  sum += mod->tree_posteriors->rcat_expected_nsites[i];
	for (i = 0; i < mod->nratecats; i++) {
	  vec_set(params, mod->ratevar_idx+i, 
		  mod->tree_posteriors->rcat_expected_nsites[i] / sum);
	}
      }
      for (i=0; i < tm_get_nratevarparams(mod); i++)
	vec_set(mod->all_params, mod->ratevar_idx+i, vec_get(params, mod->ratevar_idx+i));
    }

    opt_bfgs(likelihood_func, opt_params, (void*)mod, &tmp, lower_bounds,
             upper_bounds, logf, grad_func, bfgs_prec, H, NULL); 

    if (mod->nratecats != nratecats && 
        improvement < TM_EM_CONV(OPT_CRUDE_PREC) && home_stretch) {
      int nrateparams = tm_get_nratevarparams(mod);
      int old_npar = npar;
      if (logf != NULL) fprintf(logf, "Introducing rate variation.\n");
      mod->nratecats = nratecats;
      mod->alpha  = alpha;
      mod->rK[0] = rK0;
      mod->freqK[0] = freqK0;
      if (opt_ratevar_freqs && !mod->empirical_rates) {
	npar += nrateparams;
	vec_realloc(opt_params, npar);
	for (i=0; i < nrateparams; i++) {
	  mod->param_map[mod->ratevar_idx+i] = i + old_npar;
	  vec_set(opt_params, i + old_npar, 
		  vec_get(mod->all_params, mod->ratevar_idx + i));
	  
	}
	if (lower_bounds != NULL) vec_free(lower_bounds);
	if (upper_bounds != NULL) vec_free(upper_bounds);
	tm_new_boundaries(&lower_bounds, &upper_bounds, npar, mod, 0);
	mat_free(H);
	H = mat_new(npar, npar);
	mat_set_identity(H);
      }
    }
  }

  mod->lnL = ll;

  if (error_file != NULL)
    tm_variance(error_file, mod, msa, mod->all_params, cat);

  /* take care of final scaling of rate matrix and branch lengths */
  branchlen_scale = 1;
  if (mod->scale_during_opt == 0 && mod->subst_mod != JC69 && 
      mod->subst_mod != F81) {
    branchlen_scale *= tm_scale_rate_matrix(mod); 
    tm_scale_params(mod, params, branchlen_scale); 
    /* FIXME: this inserted so params reflect model on exit, but won't
       be correct if either condition below holds -- not currently an
       issue because they aren't supported in phyloBoot anyway */
  }
  if (mod->estimate_branchlens == TM_SCALE_ONLY) { 
    branchlen_scale *= mod->scale;
    mod->scale = 1;
  }
  if (mod->nratecats > 1 && mod->empirical_rates) { 
    double rate_scale = 0;
    for (i = 0; i < mod->nratecats; i++) /* scale for expected rate */
      rate_scale += mod->rK[i] * mod->freqK[i];
    branchlen_scale *= rate_scale;
  }
  if (branchlen_scale != 1)
    tm_scale_branchlens(mod, branchlen_scale, 0); 

  /* close log */
  if (logf != NULL) {
    gettimeofday(&end_time, NULL);
    fprintf(logf, "\nNumber of iterations: %d\nTotal time: %.4f sec.\n", it, 
            end_time.tv_sec - start_time.tv_sec + 
            (end_time.tv_usec - start_time.tv_usec)/1.0e6);
  }

  vec_free(lower_bounds);
  tl_free_tree_posteriors(mod, msa, mod->tree_posteriors);
  mod->tree_posteriors = NULL;

  vec_free(opt_params);
  return retval;
}


void remove_ratevar_from_param_map(TreeModel *mod, Vector *params) {
  int i, j;
  if (mod->nratecats == 1) return;

  for (i=mod->ratevar_idx; i < mod->ratevar_idx + tm_get_nratevarparams(mod); i++) {
    if (mod->param_map[i] >= 0) {
      for (j=0; j < params->size; j++)
	if (mod->param_map[j] > mod->param_map[i])
	  mod->param_map[j]--;
      mod->param_map[i] = -1;
    }
  }
}


double tm_partial_ll_wrapper(Vector *params, void *data) {
  TreeModel *mod = (TreeModel*)data;
  TreePosteriors *post = mod->tree_posteriors;
  tm_unpack_params(mod, params, -1);
  return -tl_compute_partial_ll_suff_stats(mod, post) * log(2);
}

/* Print a line to a log file that describes the state of the
   optimization procedure on a given iteration.  The value of the
   function is output, along with the values of all parameters.  If
   "header_only == 1", an appropriate header is printed. */
void tm_log_em(FILE *logf, int header_only, double val, Vector *params) {
  int i;
  char tmp[30];
  if (header_only) {
    fprintf(logf, "%15s ", "f(x)");
    for (i = 0; i < params->size; i++) {
      sprintf(tmp, "x_%d", i);
      fprintf(logf, "%15s ", tmp);
    }
    fprintf(logf, "\n");
  }
  else {
    fprintf(logf, "%15.6f ", val);
    for (i = 0; i < params->size; i++) 
      fprintf(logf, "%15.6f ", vec_get(params, i));
    fprintf(logf, "\n");
  }
  fflush(logf);
}

/* (used in compute_grad_em_approx) given model info and a state number,
   obtain the "neighbors" of the state -- that is, states
   corresponding to all character tuples that differ from it by no
   more than one character */
void get_neighbors(int *neighbors, int state, int order, int alph_size) {
  int place, j, k;
  for (place = 0; place <= order; place++) {
    int state_digit = (state % int_pow(alph_size, place+1)) / 
      int_pow(alph_size, place);   
    int refval = state - state_digit * int_pow(alph_size, place);
    for (j = 0, k = 0; j < alph_size; j++) {
      if (j == state_digit) continue;
      neighbors[place * (alph_size-1) + k] = refval + j *
        int_pow(alph_size, place);
      k++;
    }
  }        

  neighbors[(order+1) * (alph_size-1)] = state;  
                                /* every state is a neighbor of itself */
}

/* Compute gradient for tree model using (approximate) analytical rate
   matrix derivs.  NOTE: the tree model is assumed to be up to date
   wrt the parameters, including the eigenvalues and eigenvectors, and
   the exponentiated matrices associated with each edge (okay if
   tm_unpack_params called since last parameter update) */
void compute_grad_em_approx(Vector *grad, Vector *params, void *data, 
                          Vector *lb, Vector *ub) {

  TreeModel *mod = (TreeModel*)data;
  MarkovMatrix *P, *Q = mod->rate_matrix;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int ndigits = mod->order + 1;
  int nneighbors = ((alph_size-1) * ndigits + 1);
  int neighbors[nstates][nneighbors];
  int mark_col[nstates];
  int i, j, k, l, m, rcat, params_idx, node, lidx, lidx2, orig_size;
  int grad_idx, idx, numpar;
  TreeNode *n;
  List *traversal;
  double t;
  double freqK[mod->nratecats], rK_tweak[mod->nratecats];

  List *erows = lst_new_int(4), *ecols = lst_new_int(4), 
    *distinct_rows = lst_new_int(2), *distinct_cols = lst_new_int(4);

  static double **q = NULL, **q2 = NULL, **q3 = NULL, 
    **dq = NULL, **dqq = NULL, **qdq = NULL, **dqq2 = NULL, **qdqq = NULL, 
    **q2dq = NULL, **dqq3 = NULL, **qdqq2 = NULL, **q2dqq = NULL, 
    **q3dq = NULL;
  static Complex *diag = NULL;

  if  (Q->evals_z == NULL || Q->evec_matrix_z == NULL || Q->evec_matrix_inv_z == NULL)
    die("ERRROR: compute_grad_em_approx got NULL value in eigensystem; error diagonalizing matrix.");

  if (diag == NULL) {
    diag = (Complex*)smalloc(nstates * sizeof(Complex));
    set_static_var((void**)&diag);
  }

  /* init memory (first time only) */
  if (q == NULL) {
    q = (double**)smalloc(nstates * sizeof(double*));
    set_static_var((void**)&q);
    q2 = (double**)smalloc(nstates * sizeof(double*));
    q3 = (double**)smalloc(nstates * sizeof(double*));
    dq = (double**)smalloc(nstates * sizeof(double*));
    dqq = (double**)smalloc(nstates * sizeof(double*));
    qdq = (double**)smalloc(nstates * sizeof(double*));
    dqq2 = (double**)smalloc(nstates * sizeof(double*));
    qdqq = (double**)smalloc(nstates * sizeof(double*));
    q2dq = (double**)smalloc(nstates * sizeof(double*));
    dqq3 = (double**)smalloc(nstates * sizeof(double*));
    qdqq2 = (double**)smalloc(nstates * sizeof(double*));
    q2dqq = (double**)smalloc(nstates * sizeof(double*));
    q3dq = (double**)smalloc(nstates * sizeof(double*));

    for (i = 0; i < nstates; i++) {
      q[i] = (double*)smalloc(nstates * sizeof(double));
      q2[i] = (double*)smalloc(nstates * sizeof(double));
      q3[i] = (double*)smalloc(nstates * sizeof(double));
      dq[i] = (double*)smalloc(nstates * sizeof(double));
      dqq[i] = (double*)smalloc(nstates * sizeof(double));
      qdq[i] = (double*)smalloc(nstates * sizeof(double));
      dqq2[i] = (double*)smalloc(nstates * sizeof(double));
      qdqq[i] = (double*)smalloc(nstates * sizeof(double));
      q2dq[i] = (double*)smalloc(nstates * sizeof(double));
      q3dq[i] = (double*)smalloc(nstates * sizeof(double));
      q2dqq[i] = (double*)smalloc(nstates * sizeof(double));
      qdqq2[i] = (double*)smalloc(nstates * sizeof(double));
      dqq3[i] = (double*)smalloc(nstates * sizeof(double));
    }
  }
  
  /* set Q, zero Q^2 and Q^3 */
  for (i = 0; i < nstates; i++) {
    for (j = 0; j < nstates; j++) {
        q[i][j] = mm_get(mod->rate_matrix, i, j); /* just for
                                                     convenience
                                                     below */
        q2[i][j] = q3[i][j] = 0;
    }
  }

  /* obtain the neighbors of each state, used throughout to improve
     efficiency of multiplications */
  for (i = 0; i < nstates; i++) 
    get_neighbors(neighbors[i], i, mod->order, alph_size);

  /* compute Q^2 */
  for (i = 0; i < nstates; i++) {
    for (j = 0; j < nstates; j++) 
      for (k = 0; k < nneighbors; k++) 
        q2[i][j] += q[i][neighbors[i][k]] * q[neighbors[i][k]][j];
  }
  /* you can do this a bit more efficiently, but I don't think it's
     worth the trouble */

  /* compute Q^3 */
  for (i = 0; i < nstates; i++) {
    for (j = 0; j < nstates; j++) 
      for (k = 0; k < nneighbors; k++) 
        q3[i][j] += q[i][neighbors[i][k]] * q2[neighbors[i][k]][j];
  }

  vec_zero(grad);

  /* compute partial derivs for branch length params */
  traversal = tr_preorder(mod->tree); /* branch-length parameters
                                         correspond to pre-order
                                         traversal of tree */
  idx = 0;
  for (j = 0; j < lst_size(traversal); j++) {

    n = lst_get_ptr(traversal, j);
    if (n == mod->tree) continue;
    params_idx = mod->bl_idx + idx++;
    grad_idx = mod->param_map[params_idx];
    if (grad_idx < 0) continue;

    for (rcat = 0; rcat < mod->nratecats; rcat++) {
      P = mod->P[n->id][rcat];
      t = n->dparent * mod->rK[rcat]; /* the factor of 1/2 is taken
                                         care of here, in the def. of
                                         n->dparent */

      /* main diagonal of matrix of eigenvalues * exponentials of
         eigenvalues for branch length t*/
      for (i = 0; i < nstates; i++)
        diag[i] = z_mul_real(z_mul(z_exp(z_mul_real(zvec_get(Q->evals_z, i), t)), zvec_get(Q->evals_z, i)), mod->rK[rcat]);

      /* save time by only using complex numbers in the inner loop if
         necessary (each complex mult equivalent to four real mults and
         two real adds) */
      if (tm_node_is_reversible(mod, n)) {
        for (k = 0; k < nstates; k++) {
          for (l = 0; l < nstates; l++) {
            double p = mm_get(P, k, l);
            double dp_dt = 0;
            double dp_dt_div_p;

            for (i = 0; i < nstates; i++) 
              dp_dt += (zmat_get(Q->evec_matrix_z, k, i)).x * (diag[i]).x * 
                (zmat_get(Q->evec_matrix_inv_z, i, l)).x;

            /* have to handle case of p == 0 carefully -- want contrib
               to derivative to be zero if dp_dt == 0 or
               expected_nsubst_tot == 0 (as will normally be the case)
               and want to avoid a true inf value */
            if (p == 0) {
              if (dp_dt == 0) dp_dt_div_p = 0;
              else if (dp_dt < 0) dp_dt_div_p = NEGINFTY;
              else dp_dt_div_p = INFTY;
            }
            else dp_dt_div_p = dp_dt / p;
          
            vec_set(grad, grad_idx, vec_get(grad, grad_idx) +
                           dp_dt_div_p * 
                           mod->tree_posteriors->expected_nsubst_tot[rcat][k][l][n->id]); 

          }
        }
      }
      else {                      /* non-reversible model -- need to
                                     allow for complex numbers */
        for (k = 0; k < nstates; k++) {
          for (l = 0; l < nstates; l++) {
            double p = mm_get(P, k, l);
            double partial_p_div_p;
            Complex partial_p = z_set(0, 0);

            for (i = 0; i < nstates; i++) 
              partial_p = z_add(partial_p, z_mul(z_mul(zmat_get(Q->evec_matrix_z, k, i), diag[i]), zmat_get(Q->evec_matrix_inv_z, i, l)));
          
            if (! (fabs(partial_p.y) <= TM_IMAG_EPS))
	      die("ERROR compute_grad_em_approx: fabs(partial_p.y=%e) should be <= %e\n",
		  partial_p.y, TM_IMAG_EPS);

            /* see comments for real case (above) */
            if (p == 0) {
              if (partial_p.x == 0) partial_p_div_p = 0;
              else if (partial_p.x < 0) partial_p_div_p = NEGINFTY;
              else partial_p_div_p = INFTY;
            }
            else partial_p_div_p = partial_p.x / p;

            vec_set(grad, grad_idx, vec_get(grad, grad_idx) +
                           partial_p_div_p * 
                           mod->tree_posteriors->expected_nsubst_tot[rcat][k][l][n->id]);
          }
        }
      }
    }
  }

  /* compute partial deriv for alpha (if dgamma) */
  if (mod->nratecats > 1 && !mod->empirical_rates) {
    grad_idx = mod->param_map[mod->ratevar_idx];
    if (grad_idx < 0)
      die("ERROR compute_grad_em_approx: grad_idx=%i should be >= 0\n",
	  grad_idx);
    /* for numerical est. of derivatives of rate consts wrt alpha */
    DiscreteGamma(freqK, rK_tweak, mod->alpha + DERIV_EPSILON, 
                  mod->alpha + DERIV_EPSILON, mod->nratecats, 0);

    for (rcat = 0; rcat < mod->nratecats; rcat++) {
      double dr_da = (rK_tweak[rcat] - mod->rK[rcat]) / DERIV_EPSILON;

      for (j = 0; j < mod->tree->nnodes; j++) {
        n = lst_get_ptr(mod->tree->nodes, j);
        if (n->parent == NULL) continue;
        t = n->dparent * mod->rK[rcat];
	if (n->id == mod->root_leaf_id) {
	  if (t != 0.0)
	    die("ERROR compute_grad_em_approx: t != 0\n");
	  continue;
	}
        P = mod->P[n->id][rcat];

        for (i = 0; i < nstates; i++)
          diag[i] = z_mul_real(z_mul(z_exp(z_mul_real(zvec_get(Q->evals_z, i), t)), zvec_get(Q->evals_z, i)), n->dparent * dr_da);

        /* only use complex numbers if necessary (as above) */
        if (tm_node_is_reversible(mod, n)) {
          for (k = 0; k < nstates; k++) {
            for (l = 0; l < nstates; l++) {
              double p = mm_get(P, k, l);
              double dp_da = 0; 
              double dp_da_div_p;

              for (i = 0; i < nstates; i++) 
                dp_da += 
                  (zmat_get(Q->evec_matrix_z, k, i)).x * diag[i].x * 
                  (zmat_get(Q->evec_matrix_inv_z, i, l)).x;
 
              if (p == 0) {
                if (dp_da == 0) dp_da_div_p = 0;
                else if (dp_da < 0) dp_da_div_p = NEGINFTY;
                else dp_da_div_p = INFTY;
              }
              else dp_da_div_p = dp_da / p;

              vec_set(grad, grad_idx, 
		      vec_get(grad, grad_idx) +
		      dp_da_div_p * 
		      mod->tree_posteriors->expected_nsubst_tot[rcat][k][l][n->id]); 
            }
          }
        }
        else {                  /* non-reversible model -- use complex
                                   numbers */
          for (k = 0; k < nstates; k++) {
            for (l = 0; l < nstates; l++) {
              double p = mm_get(P, k, l);
              double dp_da_div_p;
              Complex dp_da = z_set(0, 0);
              for (i = 0; i < nstates; i++) 
                dp_da = z_add(dp_da, z_mul(z_mul(zmat_get(Q->evec_matrix_z, k, i), diag[i]), zmat_get(Q->evec_matrix_inv_z, i, l)));

              if (!(fabs(dp_da.y) <= TM_IMAG_EPS))
		die("ERROR compute_grad_em_approx: fabs(dp_da.y=%i) should be <= %e\n", dp_da.y, TM_IMAG_EPS);
	      
              if (p == 0) {
                if (dp_da.x == 0) dp_da_div_p = 0;
                else if (dp_da.x < 0) dp_da_div_p = NEGINFTY;
                else dp_da_div_p = INFTY;
              }
              else dp_da_div_p = dp_da.x / p;

              vec_set(grad, grad_idx, vec_get(grad, grad_idx) +
		      dp_da_div_p * 
		      mod->tree_posteriors->expected_nsubst_tot[rcat][k][l][n->id]); 
            }
          }
        }
      }
    }
  }
  else if (mod->empirical_rates && (mod->nratecats > 1 || mod->alpha < 0)) {
                                /* empirical rates -- gradient is zero
                                   wrt rate weights (they're already
                                   incorporated into the post
                                   probs) */
    int nrc = mod->nratecats > 1 ? mod->nratecats : -mod->alpha;
    for (; nrc >= 1; nrc--) {
      grad_idx = mod->param_map[mod->ratevar_idx + nrc - 1];
      if (grad_idx >= 0)
	vec_set(grad, grad_idx, 0);
    }
  }
  else if (mod->alpha < 0) {     /* dgamma temporarily disabled -- grad
                                   for alpha is zero */
    grad_idx = mod->param_map[mod->ratevar_idx];
    if (grad_idx >= 0)
      vec_set(grad, grad_idx, 0);
  }

  /* compute partial derivs for rate matrix params */

  /* FIXME: temporary */
  if (mod->subst_mod == JC69 || mod->subst_mod == K80 ||
      mod->subst_mod == UNDEF_MOD)
    die("ERROR compute_grad_em_approx: bad model\n");

  numpar = tm_get_nratematparams(mod);
  for (idx=0; idx<numpar; idx++) {
    params_idx = mod->ratematrix_idx + idx;
    grad_idx = mod->param_map[params_idx];
    if (grad_idx < 0) continue;
    /* zero all matrices */
    for (i = 0; i < nstates; i++)
      for (j = 0; j < nstates; j++)
       dq[i][j] = dqq[i][j] = qdq[i][j] = dqq2[i][j] = 
         qdqq[i][j] = q2dq[i][j] = dqq3[i][j] = 
         qdqq2[i][j] = q2dqq[i][j] = q3dq[i][j] = 0;

    /* element coords (rows/col pairs) at which current param appears in Q */
    lst_cpy(erows, mod->rate_matrix_param_row[params_idx]);
    lst_cpy(ecols, mod->rate_matrix_param_col[params_idx]);
    if (lst_size(erows) != lst_size(ecols))
      die("ERROR compute_grad_em_approx: size of erows (%i) does not match size of ecols (%i)\n", erows, ecols);

    /* set up dQ, the partial deriv of Q wrt the current param */
    for (i = 0; i < nstates; i++) mark_col[i] = 0;
    lst_clear(distinct_rows);
    lst_clear(distinct_cols);
    for (i = 0, orig_size = lst_size(erows); i < orig_size; i++) {
      l = lst_get_int(erows, i); 
      m = lst_get_int(ecols, i);

      if (dq[l][m] != 0)    /* row/col pairs should be unique */
	die("ERROR compute_grad_em_approx: dq[%i][%i] should be zero but is %e\n", l, m, dq[l][m]);

      dq[l][m] = subst_mod_is_reversible(mod->subst_mod) ? 
        vec_get(mod->backgd_freqs, m) : 1;
                                /* FIXME: may need to generalize */
      
      if (dq[l][m] == 0) continue; 
      /* possible if reversible with zero eq freq */

      /* keep track of distinct rows and cols with non-zero entries */
      /* also add diagonal elements to 'rows' and 'cols' lists, as
         necessary */
      if (dq[l][l] == 0) {      /* new row */
        lst_push_int(distinct_rows, l);
        lst_push_int(erows, l);
        lst_push_int(ecols, l);
      }
      if (!mark_col[m]) {       /* new col */
        lst_push_int(distinct_cols, m); 
        mark_col[m] = 1; 
      }
      if (!mark_col[l]) {       /* row also col, because of diag elment */
        lst_push_int(distinct_cols, l); 
        mark_col[l] = 1; 
      }

      dq[l][l] -= dq[l][m];     /* note that a param can appear
                                   multiple times in a row */
    }

    /* compute (dQ)Q */
    for (lidx = 0; lidx < lst_size(erows); lidx++) {
      i = lst_get_int(erows, lidx);
      k = lst_get_int(ecols, lidx);
      for (j = 0; j < nstates; j++)
        dqq[i][j] += dq[i][k] * q[k][j];
    }

    /* compute Q(dQ) */
    for (lidx = 0; lidx < lst_size(erows); lidx++) {
      k = lst_get_int(erows, lidx);
      j = lst_get_int(ecols, lidx);
      for (i = 0; i < nstates; i++)
        qdq[i][j] += q[i][k] * dq[k][j];
    }

    /* compute (dQ)Q^2 */
    for (lidx = 0; lidx < lst_size(erows); lidx++) {
      i = lst_get_int(erows, lidx);
      k = lst_get_int(ecols, lidx);
      for (j = 0; j < nstates; j++)
        dqq2[i][j] += dq[i][k] * q2[k][j];
    }

    /* compute Q(dQ)Q */
    for (lidx = 0; lidx < lst_size(distinct_rows); lidx++) {
      k = lst_get_int(distinct_rows, lidx);
      for (lidx2 = 0; lidx2 < nneighbors; lidx2++) {
        i = neighbors[k][lidx2];
        for (j = 0; j < nstates; j++)
          qdqq[i][j] += q[i][k] * dqq[k][j];
      }
    }

    /* compute Q^2(dQ) */
    for (lidx = 0; lidx < lst_size(erows); lidx++) {
      k = lst_get_int(erows, lidx);
      j = lst_get_int(ecols, lidx);
      for (i = 0; i < nstates; i++)
        q2dq[i][j] += q2[i][k] * dq[k][j];
    }

    /* compute (dQ)Q^3 */
    for (lidx = 0; lidx < lst_size(erows); lidx++) {
      i = lst_get_int(erows, lidx);
      k = lst_get_int(ecols, lidx);
      for (j = 0; j < nstates; j++)
        dqq3[i][j] += dq[i][k] * q3[k][j];
    }

    /* compute Q(dQ)Q^2 */
    for (lidx = 0; lidx < lst_size(distinct_rows); lidx++) {
      k = lst_get_int(distinct_rows, lidx);
      for (lidx2 = 0; lidx2 < nneighbors; lidx2++) {
        i = neighbors[k][lidx2];
        for (j = 0; j < nstates; j++)
          qdqq2[i][j] += q[i][k] * dqq2[k][j];
      }
    }

    /* compute Q^2(dQ)Q */
    for (lidx = 0; lidx < lst_size(distinct_cols); lidx++) {
      k = lst_get_int(distinct_cols, lidx);
      for (lidx2 = 0; lidx2 < nneighbors; lidx2++) {
        j = neighbors[k][lidx2];
        for (i = 0; i < nstates; i++)
          q2dqq[i][j] += q2dq[i][k] * q[k][j];
      }
    }

    /* compute Q^3(dQ) */
    for (lidx = 0; lidx < lst_size(erows); lidx++) {
      k = lst_get_int(erows, lidx);
      j = lst_get_int(ecols, lidx);
      for (i = 0; i < nstates; i++)
        q3dq[i][j] += q3[i][k] * dq[k][j];
    }

    for (rcat = 0; rcat < mod->nratecats; rcat++) {
      for (node = 0; node < mod->tree->nnodes; node++) {
        double taylor2, taylor3, taylor4;

        if (node == mod->tree->id)
          continue; 
        n = lst_get_ptr(mod->tree->nodes, node);
        t = n->dparent * mod->rK[rcat];
	if (n->id == mod->root_leaf_id) {
	  if (t != 0.0)
	    die("ERROR compute_grad_em_approx: t should be 0.0 but is %e\n", t);
	  continue;
	}

        P = mod->P[n->id][rcat];

        /* this allows for fewer mults in intensive loop below */
        taylor2 = t*t/2;
        taylor3 = t*t*t/6;
        taylor4 = t*t*t*t/24;

        for (i = 0; i < nstates; i++) {
          for (j = 0; j < nstates; j++) {
            double partial_p = t * dq[i][j] + 
              taylor2 * (dqq[i][j] + qdq[i][j]) +
              taylor3 * (dqq2[i][j] + qdqq[i][j] + q2dq[i][j]) +
              taylor4 * (dqq3[i][j] + qdqq2[i][j] + q2dqq[i][j] + 
                         q3dq[i][j]);
            double p = mm_get(P, i, j);
            double partial_p_div_p;

            /* handle case of p == 0 carefully, as described above */
            if (p == 0) {
              if (partial_p == 0) partial_p_div_p = 0;
              else if (partial_p < 0) partial_p_div_p = NEGINFTY;
              else partial_p_div_p = INFTY;
            }
            else partial_p_div_p = partial_p / p;
            vec_set(grad, grad_idx, vec_get(grad, grad_idx) + 
		    partial_p_div_p * 
		    mod->tree_posteriors->expected_nsubst_tot[rcat][i][j][node]);
          }
        }
      }
    }
  }
  vec_scale(grad, -1);
  lst_free(erows); lst_free(ecols); lst_free(distinct_rows); 
  lst_free(distinct_cols);
}

/* Like above, but using the approach outlined by Schadt and Lange for
   computing the partial derivatives wrt rate matrix parameters.
   Slower, but gives exact results */
void compute_grad_em_exact(Vector *grad, Vector *params, void *data, 
                           Vector *lb, Vector *ub) {

  TreeModel *mod = (TreeModel*)data;
/*   int alph_size = strlen(mod->rate_matrix->states); */
  int nstates = mod->rate_matrix->size;
  int i, j, k, l, m, rcat, params_idx, node, lidx, orig_size;
  int idx, grad_idx, numpar;
  TreeNode *n;
  MarkovMatrix *P, *Q;
  List *traversal;
  List *erows = lst_new_int(4), *ecols = lst_new_int(4), 
    *distinct_rows = lst_new_int(2);
  double t;
  double freqK[mod->nratecats], rK_tweak[mod->nratecats];

  static double **dq = NULL;
  static Complex **f = NULL, **tmpmat = NULL, **sinv_dq_s = NULL;
  static Complex *diag = NULL;

  if (diag == NULL) {
    diag = (Complex*)smalloc(nstates * sizeof(Complex));
    set_static_var((void**)&diag);
  }

  Q = mod->rate_matrix;
  if (Q->evals_z == NULL || Q->evec_matrix_z == NULL ||
      Q->evec_matrix_inv_z == NULL)
    die("ERROR compute_grade_em_exact got NULL value in eigensystem; error diagonalizing rate matrix\n");

  /* init memory (first time only) */
  if (dq == NULL) {
    dq = (double**)smalloc(nstates * sizeof(double*));
    set_static_var((void**)&dq);
    f = (Complex**)smalloc(nstates * sizeof(Complex*));
    tmpmat = (Complex**)smalloc(nstates * sizeof(Complex*));
    sinv_dq_s = (Complex**)smalloc(nstates * sizeof(Complex*));

    for (i = 0; i < nstates; i++) {
      dq[i] = (double*)smalloc(nstates * sizeof(double));
      f[i] = (Complex*)smalloc(nstates * sizeof(Complex));
      tmpmat[i] = (Complex*)smalloc(nstates * sizeof(Complex));
      sinv_dq_s[i] = (Complex*)smalloc(nstates * sizeof(Complex));
    }
  }
  
  vec_zero(grad);

  /* compute partial derivs for branch length params */
  traversal = tr_preorder(mod->tree); /* branch-length parameters
                                         correspond to pre-order
                                         traversal of tree */
  idx = 0;
  for (j = 0; j < lst_size(traversal); j++) {

    n = lst_get_ptr(traversal, j);
    if (n == mod->tree) continue;
    params_idx = mod->bl_idx + idx++;
    grad_idx = mod->param_map[params_idx];
    if (grad_idx < 0) continue;
    if (n->id == mod->root_leaf_id)
      die("ERROR compute_grad_em_exact: n->id == mod->root_leaf_id = %i\n",
	  n->id);

    for (rcat = 0; rcat < mod->nratecats; rcat++) {
      P = mod->P[n->id][rcat];
      t = n->dparent * mod->rK[rcat]; /* the factor of 1/2 is taken
                                         care of here, in the def. of
                                         n->dparent */

      /* main diagonal of matrix of eigenvalues * exponentials of
         eigenvalues for branch length t*/
      for (i = 0; i < nstates; i++)
        diag[i] = z_mul_real(z_mul(z_exp(z_mul_real(zvec_get(Q->evals_z, i), t)), zvec_get(Q->evals_z, i)), mod->rK[rcat]);

      /* save time by only using complex numbers in the inner loop if
         necessary (each complex mult equivalent to four real mults and
         two real adds) */
      if (tm_node_is_reversible(mod, n)) {
        for (k = 0; k < nstates; k++) {
          for (l = 0; l < nstates; l++) {
            double p = mm_get(P, k, l);
            double dp_dt = 0;
            double dp_dt_div_p;

            for (i = 0; i < nstates; i++) 
              dp_dt += (zmat_get(Q->evec_matrix_z, k, i)).x * diag[i].x * 
                (zmat_get(Q->evec_matrix_inv_z, i, l)).x;

            /* have to handle case of p == 0 carefully -- want contrib
               to derivative to be zero if dp_dt == 0 or
               expected_nsubst_tot == 0 (as will normally be the case)
               and want to avoid a true inf value */
            if (p == 0) {
              if (dp_dt == 0) dp_dt_div_p = 0;
              else if (dp_dt < 0) dp_dt_div_p = NEGINFTY;
              else dp_dt_div_p = INFTY;
            }
            else dp_dt_div_p = dp_dt / p;
          
            vec_set(grad, grad_idx, vec_get(grad, grad_idx) +
		    dp_dt_div_p *
		    mod->tree_posteriors->expected_nsubst_tot[rcat][k][l][n->id]);

          }
        }
      }
      else {                      /* non-reversible model -- need to
                                     allow for complex numbers */
        for (k = 0; k < nstates; k++) {
          for (l = 0; l < nstates; l++) {
            double p = mm_get(P, k, l);
            double dp_dt_div_p;
            Complex dp_dt = z_set(0, 0);

            for (i = 0; i < nstates; i++) 
              dp_dt = z_add(dp_dt, z_mul(z_mul(zmat_get(Q->evec_matrix_z, k, i), diag[i]), zmat_get(Q->evec_matrix_inv_z, i, l)));

            /* see comments for real case (above) */
            if (p == 0) {
              if (dp_dt.x == 0) dp_dt_div_p = 0;
              else if (dp_dt.x < 0) dp_dt_div_p = NEGINFTY;
              else dp_dt_div_p = INFTY;
            }
            else dp_dt_div_p = dp_dt.x / p;

            if (!(fabs(dp_dt.y) <= TM_IMAG_EPS))
	      die("ERROR compute_grad_exact: fabs(dp_dt.y=%e) should be <= %e\n",
		  dp_dt.y, TM_IMAG_EPS);
            vec_set(grad, grad_idx, vec_get(grad, grad_idx) +
                           dp_dt_div_p *
                           mod->tree_posteriors->expected_nsubst_tot[rcat][k][l][n->id]);
          }
        }
      }
    }
  }

  /* compute partial deriv for alpha (if dgamma) */
  if (mod->nratecats > 1 && !mod->empirical_rates) {
    grad_idx = mod->param_map[mod->ratevar_idx];
    if (grad_idx >= 0) {
      /* for numerical est. of derivatives of rate consts wrt alpha */
      DiscreteGamma(freqK, rK_tweak, mod->alpha + DERIV_EPSILON, 
		    mod->alpha + DERIV_EPSILON, mod->nratecats, 0);
      
      for (rcat = 0; rcat < mod->nratecats; rcat++) {
	double dr_da = (rK_tweak[rcat] - mod->rK[rcat]) / DERIV_EPSILON;
	for (j = 0; j < mod->tree->nnodes; j++) {
	  n = lst_get_ptr(mod->tree->nodes, j);
	  if (n->parent == NULL) continue;
	  t = n->dparent * mod->rK[rcat];
	  if (n->id == mod->root_leaf_id) {
	    if (t != 0.0)
	      die("ERROR compute_grad_exact: t should be zero but is %e\n", t);
	    continue;
	  }
	  P = mod->P[n->id][rcat];
	  
	  for (i = 0; i < nstates; i++)
	    diag[i] = z_mul_real(z_mul(z_exp(z_mul_real(zvec_get(Q->evals_z, i), t)), zvec_get(Q->evals_z, i)), n->dparent * dr_da);
	  
	  /* only use complex numbers if necessary (as above) */
	  if (tm_node_is_reversible(mod, n)) {
	    for (k = 0; k < nstates; k++) {
	      for (l = 0; l < nstates; l++) {
		double p = mm_get(P, k, l);
		double dp_da = 0; 
		double dp_da_div_p;
		
		for (i = 0; i < nstates; i++) 
		  dp_da += 
		    (zmat_get(Q->evec_matrix_z, k, i)).x * diag[i].x *
		    (zmat_get(Q->evec_matrix_inv_z, i, l)).x;
		
		if (p == 0) {
		  if (dp_da == 0) dp_da_div_p = 0;
		  else if (dp_da < 0) dp_da_div_p = NEGINFTY;
		  else dp_da_div_p = INFTY;
		}
		else dp_da_div_p = dp_da / p;
		
		vec_set(grad, grad_idx, 
			vec_get(grad, grad_idx) +
			dp_da_div_p * 
			mod->tree_posteriors->expected_nsubst_tot[rcat][k][l][n->id]); 
	      }
	    }
	  }
	  else {                  /* non-reversible model -- use complex
				     numbers */
	    for (k = 0; k < nstates; k++) {
	      for (l = 0; l < nstates; l++) {
		double p = mm_get(P, k, l);
		double dp_da_div_p;
		Complex dp_da = z_set(0, 0);
		for (i = 0; i < nstates; i++) 
		  dp_da = z_add(dp_da, z_mul(z_mul(zmat_get(Q->evec_matrix_z, k, i), diag[i]), zmat_get(Q->evec_matrix_inv_z, i, l)));
		
		if (!(fabs(dp_da.y) <= TM_IMAG_EPS))
		  die("ERROR compute_grad_exact: fabs(dp_da.y=%e) should be <= %e\n",
		      dp_da.y, TM_IMAG_EPS);
		
		if (p == 0) {
		  if (dp_da.x == 0) dp_da_div_p = 0;
		  else if (dp_da.x < 0) dp_da_div_p = NEGINFTY;
		  else dp_da_div_p = INFTY;
		}
		else dp_da_div_p = dp_da.x / p;
		
		vec_set(grad, grad_idx, vec_get(grad, grad_idx) +
			dp_da_div_p * 
			mod->tree_posteriors->expected_nsubst_tot[rcat][k][l][n->id]); 
	      }
	    }
	  }
	}
      }
    }
  }
  else if (mod->empirical_rates && (mod->nratecats > 1 || mod->alpha < 0)) {
                                /* empirical rates -- gradient is zero
                                   wrt rate weights (they're already
                                   incorporated into the post
                                   probs) */
    int nrc = mod->nratecats > 1 ? mod->nratecats : -mod->alpha;
    for (; nrc >= 1; nrc--) {
      grad_idx = mod->param_map[mod->ratevar_idx + nrc - 1];
      if (grad_idx >= 0)
	vec_set(grad, grad_idx, 0);
    }
  }
  else if (mod->alpha < 0) {     /* dgamma temporarily disabled -- grad
				    for alpha is zero */
    grad_idx = mod->param_map[mod->ratevar_idx];
    if (grad_idx >= 0)
      vec_set(grad, grad_idx, 0);
  }

  /* compute partial derivs for rate matrix params */

  /* FIXME: temporary */
  if (mod->subst_mod == JC69 || mod->subst_mod == K80 ||
      mod->subst_mod == UNDEF_MOD)
    die("ERROR compute_grad_exact: bad subst mod\n");

  numpar = tm_get_nratematparams(mod);
  for (idx=0; idx < numpar; idx++) {
    params_idx = mod->ratematrix_idx + idx;
    grad_idx = mod->param_map[params_idx];
    if (grad_idx < 0) continue;

    for (i = 0; i < nstates; i++) {
      for (j = 0; j < nstates; j++) {
        dq[i][j] = 0;
	tmpmat[i][j] = z_set(0, 0);
	sinv_dq_s[i][j] = z_set(0, 0);
      }
    }

    /* element coords (rows/col pairs) at which current param appears in Q */
    lst_cpy(erows, mod->rate_matrix_param_row[params_idx]);
    lst_cpy(ecols, mod->rate_matrix_param_col[params_idx]);
    if (lst_size(erows) != lst_size(ecols))
      die("ERROR compute_grad_exact: size of erows (%i) does not match size of ecols (%i)\n", lst_size(erows), lst_size(ecols));

    /* set up dQ, the partial deriv of Q wrt the current param */
    lst_clear(distinct_rows);
    for (i = 0, orig_size = lst_size(erows); i < orig_size; i++) {
      l = lst_get_int(erows, i); 
      m = lst_get_int(ecols, i);

      if (dq[l][m] != 0)    /* row/col pairs should be unique */
	die("ERROR compute_grad_exact dq[%i][%i] should be zero but is %e\n",
	    l, m, dq[l][m]);

      dq[l][m] = subst_mod_is_reversible(mod->subst_mod) ? 
        vec_get(mod->backgd_freqs, m) : 1;
                                /* FIXME: may need to generalize */
      
      if (dq[l][m] == 0) continue; 
      /* possible if reversible with zero eq freq */

      /* keep track of distinct rows and cols with non-zero entries */
      /* also add diagonal elements to 'rows' and 'cols' lists, as
         necessary */
      if (dq[l][l] == 0) {      /* new row */
        lst_push_int(distinct_rows, l);
        lst_push_int(erows, l);
        lst_push_int(ecols, l);
      }

      dq[l][l] -= dq[l][m];     /* note that a param can appear
                                   multiple times in a row */
    }

    /* compute S^-1 dQ S */
    for (lidx = 0; lidx < lst_size(erows); lidx++) {
      i = lst_get_int(erows, lidx);
      k = lst_get_int(ecols, lidx);
      for (j = 0; j < nstates; j++)
        tmpmat[i][j] = z_add(tmpmat[i][j], z_mul_real(zmat_get(Q->evec_matrix_z, k, j), dq[i][k]));
    }

    for (lidx = 0; lidx < lst_size(distinct_rows); lidx++) {
      k = lst_get_int(distinct_rows, lidx);
      for (i = 0; i < nstates; i++) {
        for (j = 0; j < nstates; j++) {
          sinv_dq_s[i][j] =
            z_add(sinv_dq_s[i][j], z_mul(zmat_get(Q->evec_matrix_inv_z, i, k), tmpmat[k][j]));
        }
      }
    }

    for (rcat = 0; rcat < mod->nratecats; rcat++) {
      for (node = 0; node < mod->tree->nnodes; node++) {
        if (node == mod->tree->id)
          continue; 

        n = lst_get_ptr(mod->tree->nodes, node);
        t = n->dparent * mod->rK[rcat];
	if (n->id == mod->root_leaf_id) {
	  if (t != 0.0)
	    die("ERROR compute_grad_exact expected t to be zero but was %e\n",
		t);
	  continue;
	}
        P = mod->P[n->id][rcat];

        /* as above, it's worth it to have separate versions of the
           computations below for the real and complex cases */

        if (tm_node_is_reversible(mod, n)) { /* real case */
          /* build the matrix F */
          for (i = 0; i < nstates; i++) {
            for (j = 0; j < nstates; j++) {
              if ((zvec_get(Q->evals_z, i)).x ==
                  (zvec_get(Q->evals_z, j)).x)
		f[i][j].x = exp((zvec_get(Q->evals_z, i)).x * t) * t;
              else
                f[i][j].x = (exp((zvec_get(Q->evals_z, i)).x * t) 
			     - exp((zvec_get(Q->evals_z, j)).x * t)) /
		  ((zvec_get(Q->evals_z, i)).x - (zvec_get(Q->evals_z, j)).x);
            }
          }

          /* compute (F o S^-1 dQ S) S^-1 */
          for (i = 0; i < nstates; i++) {
            for (j = 0; j < nstates; j++) {
              tmpmat[i][j].x = 0;
              for (k = 0; k < nstates; k++) 
                tmpmat[i][j].x += f[i][k].x * sinv_dq_s[i][k].x *
		  (zmat_get(Q->evec_matrix_inv_z, k, j)).x;
            }
          }

          /* compute S (F o S^-1 dQ S) S^-1; simultaneously compute
             gradient elements */
          for (i = 0; i < nstates; i++) {
            for (j = 0; j < nstates; j++) {
              double partial_p = 0;
              double p = mm_get(P, i, j);
              double partial_p_div_p;

              for (k = 0; k < nstates; k++) 
                partial_p += (zmat_get(Q->evec_matrix_z, i, k)).x * 
                  tmpmat[k][j].x;

              /* handle case of p == 0 carefully, as described above */
              if (p == 0) {
                if (partial_p == 0) partial_p_div_p = 0;
                else if (partial_p < 0) partial_p_div_p = NEGINFTY;
                else partial_p_div_p = INFTY;
              }
              else partial_p_div_p = partial_p / p;

              vec_set(grad, grad_idx, vec_get(grad, grad_idx) + 
		      partial_p_div_p *
		      mod->tree_posteriors->expected_nsubst_tot[rcat][i][j][node]);
            }
          }
        }
        else {                    /* complex case */
          /* build the matrix F */
          for (i = 0; i < nstates; i++) {
            for (j = 0; j < nstates; j++) {
              if (z_eq(zvec_get(Q->evals_z, i),
		       zvec_get(Q->evals_z, j)))
                f[i][j] = z_mul_real(z_exp(z_mul_real(zvec_get(Q->evals_z, i), t)), t);
              else
                f[i][j] = z_div(z_sub(z_exp(z_mul_real(zvec_get(Q->evals_z, i), t)), z_exp(z_mul_real(zvec_get(Q->evals_z, j), t))), z_sub(zvec_get(Q->evals_z, i), zvec_get(Q->evals_z, j)));
	      
            }
          }

          /* compute (F o S^-1 dQ S) S^-1 */
          for (i = 0; i < nstates; i++) {
            for (j = 0; j < nstates; j++) {
	      tmpmat[i][j] = z_set(0, 0);
              for (k = 0; k < nstates; k++) 
                tmpmat[i][j] = z_add(tmpmat[i][j], z_mul(f[i][k], z_mul(sinv_dq_s[i][k], zmat_get(Q->evec_matrix_inv_z, k, j))));
            }
          }

          /* compute S (F o S^-1 dQ S) S^-1; simultaneously compute
             gradient elements */
          for (i = 0; i < nstates; i++) {
            for (j = 0; j < nstates; j++) {
              double p = mm_get(P, i, j);
              double partial_p_div_p;
              Complex partial_p = z_set(0, 0);

              for (k = 0; k < nstates; k++) 
                partial_p = z_add(partial_p, z_mul(zmat_get(Q->evec_matrix_z, i, k), tmpmat[k][j]));

              if (!(fabs(partial_p.y) <= TM_IMAG_EPS))
		die("ERROR compute_grad_exact: fabs(partial_p.y=%e) should be <= %e\n", partial_p.y, TM_IMAG_EPS);

              /* handle case of p == 0 carefully, as described above */
              if (p == 0) {
                if (partial_p.x == 0) partial_p_div_p = 0;
                else if (partial_p.x < 0) partial_p_div_p = NEGINFTY;
                else partial_p_div_p = INFTY;
              }
              else partial_p_div_p = partial_p.x / p;

              vec_set(grad, grad_idx, vec_get(grad, grad_idx) + 
		      partial_p_div_p * 
		      mod->tree_posteriors->expected_nsubst_tot[rcat][i][j][node]);
            }
          }
        }
      }
    }
  }
  vec_scale(grad, -1);
  lst_free(erows); lst_free(ecols); lst_free(distinct_rows); 
}

