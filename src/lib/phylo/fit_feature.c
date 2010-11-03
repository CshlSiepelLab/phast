/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: fit_feature.c,v 1.5 2008-12-06 16:46:31 acs Exp $ */

/* Functions to compute likelihoods, estimate scale factors, perform
   LRTs, score tests, etc. for multi-column features.  Generalization
   of fit_columm.c to features.  Makes use of several of the
   single-column functions. */

#include <fit_feature.h>
#include <fit_column.h>
#include <sufficient_stats.h>
#include <tree_likelihoods.h>

#define SIGFIGS 4
/* number of significant figures to which to estimate scale
   parameters (currently affects 1d parameter estimation only) */

/* Compute and return the log likelihood of a tree model with respect
   to a given feature.  Calls col_compute_log_likelihood for each
   column */ 
double ff_compute_log_likelihood(TreeModel *mod, MSA *msa, GFF_Feature *feat,
                                 double **scratch) {
  double retval = 0;
  int i;
  for (i = feat->start-1; i < feat->end; i++) /* offset of one */
    retval += col_compute_log_likelihood(mod, msa, msa->ss->tuple_idx[i], 
                                         scratch);
  return retval;
}

/* Compute the first and (optionally) second derivatives with respect
   to the scale parameter for the single-feature log likelihood
   function (feat_compute_log_likelihood).  This version assumes a
   single scale parameter; see below for the subtree case.  Return
   value is log likelihood, which is computed as a by-product.  Derivs
   will be stored in *first_deriv and *second_deriv.  If second_deriv
   == NULL, it will not be computed (saves some time).  */
double ff_scale_derivs(FeatFitData *d, double *first_deriv,
                       double *second_deriv, double ***scratch) {
  double d1, d2, retval = 0;
  double *pd2 = (second_deriv == NULL ? NULL : &d2);
  int i;
  (*first_deriv) = 0;
  if (second_deriv != NULL) (*second_deriv) = 0;
  for (i = d->feat->start-1; i < d->feat->end; i++) { /* offset of one */
    d->cdata->tupleidx = d->cdata->msa->ss->tuple_idx[i];
    retval += col_scale_derivs(d->cdata, &d1, pd2, scratch);
    (*first_deriv) += d1;
    if (second_deriv != NULL) (*second_deriv) += d2;
  }
  return retval;
}

/* Compute the first and (optionally) second derivatives with respect
   to the scale parameters for the single-feature log likelihood
   function (feat_compute_log_likelihood).  This version assumes scale
   parameters for the whole tree and for the subtree.  Return value is
   log likelihood, which is computed as a by-product.  Derivs will be
   stored in *gradient and *hessian.  If hessian == NULL,
   it will not be computed (saves some time).  */
double ff_scale_derivs_subtree(FeatFitData *d, Vector *gradient, 
                                 Matrix *hessian, double ***scratch) {
  double retval = 0;
  Vector *d1 = vec_new(2);
  Matrix *d2 = (hessian == NULL ? NULL : mat_new(2, 2));
  int i;
  vec_zero(gradient);
  if (hessian != NULL) mat_zero(hessian);
  for (i = d->feat->start-1; i < d->feat->end; i++) { /* offset of one */
    d->cdata->tupleidx = d->cdata->msa->ss->tuple_idx[i];
    retval += col_scale_derivs_subtree(d->cdata, d1, d2, scratch);
    vec_plus_eq(gradient, d1);
    if (hessian != NULL) mat_plus_eq(hessian, d2);
  }
  return retval;
}

/* Wrapper for likelihood function for use in parameter estimation */
double ff_likelihood_wrapper(Vector *params, void *data) {
  FeatFitData *d = (FeatFitData*)data;

  d->cdata->mod->scale = vec_get(params, 0);
  if (d->cdata->stype == SUBTREE) 
    d->cdata->mod->scale_sub = vec_get(params, 1);

  /* reestimate subst models on edges */
  tm_set_subst_matrices(d->cdata->mod); 

  return -1 * ff_compute_log_likelihood(d->cdata->mod, d->cdata->msa, 
                                        d->feat, d->cdata->fels_scratch[0]);
}

/* Wrapper for likelihood function for use in parameter estimation;
   version for use with opt_newton_1d */
double ff_likelihood_wrapper_1d(double x, void *data) {
  FeatFitData *d = (FeatFitData*)data;
  if (d->cdata->stype == SUBTREE)
    die("ERROR: ff_likelihood_wrapper_1d: d->cdata->stype is SUBTREE\n");

  d->cdata->mod->scale = x;

  /* reestimate subst models on edges */
  tm_set_subst_matrices(d->cdata->mod); 

  return -1 * ff_compute_log_likelihood(d->cdata->mod, d->cdata->msa, 
                                        d->feat, d->cdata->fels_scratch[0]);
}

/* Wrapper for gradient function for use in parameter estimation */
void ff_grad_wrapper(Vector *grad, Vector *params, void *data, 
                     Vector *lb, Vector *ub) {
  FeatFitData *d = (FeatFitData*)data;
  double deriv;

  if (d->cdata->stype == ALL) {
    ff_scale_derivs(d, &deriv, NULL, d->cdata->fels_scratch);
    vec_set(grad, 0, -deriv);   /* because working with neg lnl */
  }
  else {
    ff_scale_derivs_subtree(d, grad, NULL, d->cdata->fels_scratch);
    vec_scale(grad, -1);
  }
}

/* Perform a likelihood ratio test for each feature in a GFF,
   comparing the given null model with an alternative model that has a
   free scaling parameter for all branches.  Assumes a 0th order
   model, leaf-to-sequence mapping already available, prob matrices
   computed, sufficient stats available.  Computes p-values based
   using the chi-sq distribution and stores them in feat_pvals.  Will
   optionally store the individual scale factors in feat_scales and
   raw log likelihood ratios in feat_llrs if these variables are
   non-NULL.  Must define mode as CON (for 0 <= scale <= 1), ACC (for
   1 <= scale), NNEUT (0 <= scale), or CONACC (0 <= scale) */ 
void ff_lrts(TreeModel *mod, MSA *msa, GFF_Set *gff, mode_type mode, 
             double *feat_pvals, double *feat_scales, double *feat_llrs, 
             FILE *logf) {
  int i;
  FeatFitData *d;
  double null_lnl, alt_lnl, delta_lnl, this_scale = 1;

  /* init FeatFitData */
  d = ff_init_fit_data(mod, msa, ALL, mode, FALSE);

  /* iterate through features  */
  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *f = lst_get_ptr(gff->features, i);
    checkInterrupt();

    /* first check for actual substitution data in feature; if none,
       don't waste time computing likelihoods */
    if (!ff_has_data(mod, msa, f)) {
      delta_lnl = 0;
      this_scale = 1;
    }

    else {
      mod->scale = 1;
      tm_set_subst_matrices(mod);

      /* compute log likelihoods under null and alt hypotheses */
      null_lnl = ff_compute_log_likelihood(mod, msa, f, 
                                           d->cdata->fels_scratch[0]);

      vec_set(d->cdata->params, 0, d->cdata->init_scale);
      d->feat = f;

      opt_newton_1d(ff_likelihood_wrapper_1d, &d->cdata->params->data[0], d, 
                    &alt_lnl, SIGFIGS, d->cdata->lb->data[0], 
                    d->cdata->ub->data[0], logf, NULL, NULL);   
      /* turns out to be faster to use numerical rather than exact
         derivatives (judging by col case) */

      alt_lnl *= -1;
      this_scale = d->cdata->params->data[0];

      delta_lnl = alt_lnl - null_lnl;
      if (delta_lnl <= -0.01)
	die("ERROR ff_lrts: delta_lnl (%f) <= -0.01\n", delta_lnl);
      if (delta_lnl < 0) delta_lnl = 0;
    }

    /* compute p-vals via chi-sq */
    if (feat_pvals != NULL) {
      if (mode == NNEUT || mode == CONACC) 
        feat_pvals[i] = chisq_cdf(2*delta_lnl, 1, FALSE);
      else 
        feat_pvals[i] = half_chisq_cdf(2*delta_lnl, 1, FALSE);
      /* assumes 50:50 mix of chisq and point mass at zero, due to
         bounding of param */

      if (feat_pvals[i] < 1e-20)
        feat_pvals[i] = 1e-20;
      /* approx limit of eval of tail prob; pvals of 0 cause problems */

      if (mode == CONACC && this_scale > 1)
        feat_pvals[i] *= -1; /* mark as acceleration */
    }

    /* store scales and log likelihood ratios if necessary */
    if (feat_scales != NULL) feat_scales[i] = this_scale;
    if (feat_llrs != NULL) feat_llrs[i] = delta_lnl;
  }
  
  ff_free_fit_data(d);
  sfree(d);  //above function doesn't actually free object
}

/* Subtree version of LRT */
void ff_lrts_sub(TreeModel *mod, MSA *msa, GFF_Set *gff, mode_type mode, 
                 double *feat_pvals, double *feat_null_scales, 
                 double *feat_scales, double *feat_sub_scales, 
                 double *feat_llrs, FILE *logf) {
  int i;
  FeatFitData *d, *d2;
  double null_lnl, alt_lnl, delta_lnl;
  TreeModel *modcpy;
  List *inside=NULL, *outside=NULL;

  modcpy = tm_create_copy(mod);   /* need separate copy of tree model
                                     with different internal scaling
                                     data for supertree/subtree case */

  /* init ColFitData -- one for null model, one for alt */
  d = ff_init_fit_data(modcpy, msa, ALL, NNEUT, FALSE);
  d2 = ff_init_fit_data(mod, msa, SUBTREE, mode, FALSE); 
                                /* mod has the subtree info, modcpy
                                   does not */

  /* prepare lists of leaves inside and outside root, for use in
     checking for informative substitutions */
  if (mod->subtree_root != NULL) {
    inside = lst_new_ptr(mod->tree->nnodes);
    outside = lst_new_ptr(mod->tree->nnodes); 
    tr_partition_leaves(mod->tree, mod->subtree_root, inside, outside);
  }

  /* iterate through features  */
  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *f = lst_get_ptr(gff->features, i);
    checkInterrupt();

    /* first check for informative substitution data in feature; if none,
       don't waste time computing likelihoods */
    if (!ff_has_data_sub(mod, msa, f, inside, outside)) {
      delta_lnl = 0;
      d->cdata->params->data[0] = d2->cdata->params->data[0] = 
        d2->cdata->params->data[1] = 1;
    }

    else {
      /* compute log likelihoods under null and alt hypotheses */
      d->feat = f;
      vec_set(d->cdata->params, 0, d->cdata->init_scale);
      opt_newton_1d(ff_likelihood_wrapper_1d, &d->cdata->params->data[0], d, 
                    &null_lnl, SIGFIGS, d->cdata->lb->data[0], 
                    d->cdata->ub->data[0], logf, NULL, NULL);   
      null_lnl *= -1;

      d2->feat = f;
      vec_set(d2->cdata->params, 0, d->cdata->params->data[0]); 
                                /* init to previous estimate to save time */
      vec_set(d2->cdata->params, 1, d2->cdata->init_scale_sub);
      //      vec_set(d2->cdata->params, 1, 0.01);
      if (opt_bfgs(ff_likelihood_wrapper, d2->cdata->params, d2, &alt_lnl, 
                   d2->cdata->lb, d2->cdata->ub, logf, NULL, 
                   OPT_HIGH_PREC, NULL) != 0)
        ;                         /* do nothing; nonzero exit typically
                                     occurs when max iterations is
                                     reached; a warning is printed to
                                     the log */
      alt_lnl *= -1;

      delta_lnl = alt_lnl - null_lnl;

      /* This is a hack, it would be better to figure out why the
	 optimization sometimes fails here.
	 If we get a significantly negative lnL, re-initialize params
	 so that they are identical to null model params and re-start */
      if (delta_lnl <= -0.05) {
	d2->feat = f;
	vec_set(d2->cdata->params, 0, d->cdata->params->data[0]); 
	vec_set(d2->cdata->params, 1, 1.0);
	if (opt_bfgs(ff_likelihood_wrapper, d2->cdata->params, d2, &alt_lnl, 
		     d2->cdata->lb, d2->cdata->ub, logf, NULL, 
		     OPT_HIGH_PREC, NULL) != 0)
	  if (delta_lnl <= -0.1)
	    die("ERROR ff_lrts_sub: delta_lnl (%f) <= -0.1\n", delta_lnl);
      }
      if (delta_lnl < 0) delta_lnl = 0;
    }

    /* compute p-vals via chi-sq */
    if (feat_pvals != NULL) {
      if (mode == NNEUT || mode == CONACC) 
        feat_pvals[i] = chisq_cdf(2*delta_lnl, 1, FALSE);
      else 
        feat_pvals[i] = half_chisq_cdf(2*delta_lnl, 1, FALSE);
        /* assumes 50:50 mix of chisq and point mass at zero, due to
           bounding of param */

      if (feat_pvals[i] < 1e-20)
        feat_pvals[i] = 1e-20;
      /* approx limit of eval of tail prob; pvals of 0 cause problems */

      if (mode == CONACC && d2->cdata->params->data[1] > 1)
        feat_pvals[i] *= -1;    /* mark as acceleration */        
    }

    /* store scales and log likelihood ratios if necessary */
    if (feat_null_scales != NULL) 
      feat_null_scales[i] = d->cdata->params->data[0];
    if (feat_scales != NULL) 
      feat_scales[i] = d2->cdata->params->data[0];
    if (feat_sub_scales != NULL) 
      feat_sub_scales[i] = d2->cdata->params->data[1];
    if (feat_llrs != NULL) 
      feat_llrs[i] = delta_lnl;
  }
  
  ff_free_fit_data(d);
  ff_free_fit_data(d2);
  modcpy->estimate_branchlens = TM_BRANCHLENS_ALL; 
                                /* have to revert for tm_free to work
                                   correctly */
  tm_free(modcpy);
  if (inside != NULL) lst_free(inside);
  if (outside != NULL) lst_free(outside);
}

/* Score test */
void ff_score_tests(TreeModel *mod, MSA *msa, GFF_Set *gff, mode_type mode, 
                    double *feat_pvals, double *feat_derivs, 
                    double *feat_teststats) {
  int i;
  FeatFitData *d;
  double first_deriv, teststat, fim;

  /* init FeatFitData */
  d = ff_init_fit_data(mod, msa, ALL, NNEUT, FALSE);

  /* precompute FIM */
  fim = col_estimate_fim(mod);

  if (fim < 0) 
    die("ERROR: negative fisher information in col_score_tests\n");

  /* iterate through features  */
  for (i = 0; i < lst_size(gff->features); i++) {
    checkInterrupt();
    d->feat = lst_get_ptr(gff->features, i);

    /* first check for actual substitution data in feature; if none,
       don't waste time computing likelihoods */
    if (!ff_has_data(mod, msa, d->feat)) {
      teststat = 0;
      first_deriv = 1;
    }

    else {

      ff_scale_derivs(d, &first_deriv, NULL, d->cdata->fels_scratch);

      teststat = first_deriv*first_deriv / 
        ((d->feat->end - d->feat->start + 1) * fim);
      /* scale column-by-column FIM by length of feature (expected
         values are additive) */

      if ((mode == ACC && first_deriv < 0) ||
          (mode == CON && first_deriv > 0))
        teststat = 0;             /* derivative points toward boundary;
                                     truncate at 0 */
    }

    if (feat_pvals != NULL) {
      if (mode == NNEUT || mode == CONACC)
        feat_pvals[i] = chisq_cdf(teststat, 1, FALSE);
      else 
        feat_pvals[i] = half_chisq_cdf(teststat, 1, FALSE);
      /* assumes 50:50 mix of chisq and point mass at zero */
      
      if (feat_pvals[i] < 1e-20)
        feat_pvals[i] = 1e-20;
      /* approx limit of eval of tail prob; pvals of 0 cause problems */
        
      if (mode == CONACC && first_deriv > 0)
        feat_pvals[i] *= -1; /* mark as acceleration */
    }

    /* store scales and log likelihood ratios if necessary */
    if (feat_derivs != NULL) feat_derivs[i] = first_deriv;
    if (feat_teststats != NULL) feat_teststats[i] = teststat;
  }

  ff_free_fit_data(d);
}

/* Subtree version of score test */
void ff_score_tests_sub(TreeModel *mod, MSA *msa, GFF_Set *gff, mode_type mode,
                        double *feat_pvals, double *feat_null_scales, 
                        double *feat_derivs, double *feat_sub_derivs, 
                        double *feat_teststats, FILE *logf) {
  int i;
  FeatFitData *d, *d2;
  Vector *grad = vec_new(2);
  Matrix *fim = mat_new(2, 2);
  double lnl, teststat;
  FimGrid *grid;
  List *inside=NULL, *outside=NULL;
  TreeModel *modcpy = tm_create_copy(mod); /* need separate copy of tree model
                                              with different internal scaling
                                              data for supertree/subtree case */

  /* init FeatFitData -- one for null model, one for alt */
  d = ff_init_fit_data(modcpy, msa, ALL, NNEUT, FALSE);
  d2 = ff_init_fit_data(mod, msa, SUBTREE, NNEUT, FALSE); 
                                /* mod has the subtree info, modcpy
                                   does not */

  /* precompute Fisher information matrices for a grid of scale values */
  grid = col_fim_grid_sub(mod); 

  /* prepare lists of leaves inside and outside root, for use in
     checking for informative substitutions */
  if (mod->subtree_root != NULL) {
    inside = lst_new_ptr(mod->tree->nnodes);
    outside = lst_new_ptr(mod->tree->nnodes);
    tr_partition_leaves(mod->tree, mod->subtree_root, inside, outside);
  }

  /* iterate through features  */
  for (i = 0; i < lst_size(gff->features); i++) {
    checkInterrupt();
    d->feat = lst_get_ptr(gff->features, i);

    /* first check for informative substitution data in feature; if none,
       don't waste time computing likelihoods */
    if (!ff_has_data_sub(mod, msa, d->feat, inside, outside)) { 
      teststat = 0;
      vec_zero(grad);
    }

    else {
      vec_set(d->cdata->params, 0, d->cdata->init_scale);
      opt_newton_1d(ff_likelihood_wrapper_1d, &d->cdata->params->data[0], d, 
                    &lnl, SIGFIGS, d->cdata->lb->data[0], d->cdata->ub->data[0], 
                    logf, NULL, NULL);   
      /* turns out to be faster to use numerical rather than exact
         derivatives (judging by col case) */

      d2->feat = d->feat;
      d2->cdata->mod->scale = d->cdata->params->data[0];
      d2->cdata->mod->scale_sub = 1;
      tm_set_subst_matrices(d2->cdata->mod);
      ff_scale_derivs_subtree(d2, grad, NULL, d2->cdata->fels_scratch);

      fim = col_get_fim_sub(grid, d2->cdata->mod->scale); 
      mat_scale(fim, d->feat->end - d->feat->start + 1);
      /* scale column-by-column FIM by length of feature (expected
         values are additive) */
    
      teststat = grad->data[1]*grad->data[1] / 
        (fim->data[1][1] - fim->data[0][1]*fim->data[1][0]/fim->data[0][0]);

      if (teststat < 0) {
        fprintf(stderr, "WARNING: teststat < 0 (%f)\n", teststat);
        teststat = 0; 
      }

      if ((mode == ACC && grad->data[1] < 0) ||
          (mode == CON && grad->data[1] > 0))
        teststat = 0;             /* derivative points toward boundary;
                                     truncate at 0 */

      mat_free(fim);
    }

    if (feat_pvals != NULL) {
      if (mode == NNEUT || mode == CONACC)
        feat_pvals[i] = chisq_cdf(teststat, 1, FALSE);
      else 
        feat_pvals[i] = half_chisq_cdf(teststat, 1, FALSE);
        /* assumes 50:50 mix of chisq and point mass at zero */

      if (feat_pvals[i] < 1e-20)
        feat_pvals[i] = 1e-20;
      /* approx limit of eval of tail prob; pvals of 0 cause problems */

      if (mode == CONACC && grad->data[1] > 0)
        feat_pvals[i] *= -1; /* mark as acceleration */
    }

    /* store scales and log likelihood ratios if necessary */
    if (feat_null_scales != NULL) feat_null_scales[i] = d->cdata->params->data[0];
    if (feat_derivs != NULL) feat_derivs[i] = grad->data[0];
    if (feat_sub_derivs != NULL) feat_sub_derivs[i] = grad->data[1];
    if (feat_teststats != NULL) feat_teststats[i] = teststat;
  }

  ff_free_fit_data(d);
  ff_free_fit_data(d2);
  vec_free(grad); 
  modcpy->estimate_branchlens = TM_BRANCHLENS_ALL; 
                                /* have to revert for tm_free to work
                                   correctly */
  tm_free(modcpy);
  col_free_fim_grid(grid); 
  if (inside != NULL) lst_free(inside);
  if (outside != NULL) lst_free(outside);
}

/* Perform a GERP-like computation for each feature.  Computes expected
   number of subst. under neutrality (feat_nneut), expected number
   after rescaling by ML (feat_nobs), expected number of rejected
   substitutions (feat_nrejected), and number of species with data
   (feat_nspecies).  If any arrays are NULL, values will not be
   retained.  Gaps and missing data are handled by working with
   induced subtree.  */
void ff_gerp(TreeModel *mod, MSA *msa, GFF_Set *gff, mode_type mode, 
             double *feat_nneut, double *feat_nobs, double *feat_nrejected, 
             double *feat_nspec, FILE *logf) { 
  int i, j, nspec = 0;
  double nneut, scale, lnl;
  int *has_data = smalloc(mod->tree->nnodes * sizeof(int));
  FeatFitData *d;

  /* init FeatFitData */
  d = ff_init_fit_data(mod, msa, ALL, NNEUT, FALSE);

  /* iterate through features  */
  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *f = lst_get_ptr(gff->features, i);
    checkInterrupt();
    ff_find_missing_branches(mod, msa, f, has_data, &nspec);

    if (nspec < 3) 
      nneut = scale = 0;
    else {
      vec_set(d->cdata->params, 0, d->cdata->init_scale);
      d->feat = f;

      opt_newton_1d(ff_likelihood_wrapper_1d, &d->cdata->params->data[0], d, 
                    &lnl, SIGFIGS, d->cdata->lb->data[0], d->cdata->ub->data[0], 
                    logf, NULL, NULL);   
      /* turns out to be faster to use numerical rather than exact
         derivatives (judging by col case) */
      
      scale = d->cdata->params->data[0];
      for (j = 1, nneut = 0; j < mod->tree->nnodes; j++)  /* node 0 is root */
        if (has_data[j]) 
          nneut += ((TreeNode*)lst_get_ptr(mod->tree->nodes, j))->dparent;
    }

    if (feat_nspec != NULL) feat_nspec[i] = (double)nspec;
    if (feat_nneut != NULL) feat_nneut[i] = nneut;
    if (feat_nobs != NULL) feat_nobs[i] = scale * nneut;
    if (feat_nrejected != NULL) {
      feat_nrejected[i] = nneut * (1 - scale);
      if (mode == ACC) feat_nrejected[i] *= -1;
      else if (mode == NNEUT) feat_nrejected[i] = fabs(feat_nrejected[i]);
    }
  }
  ff_free_fit_data(d);
  sfree(has_data);
}

/* Create object with metadata and scratch memory for fitting scale
   factors */
FeatFitData *ff_init_fit_data(TreeModel *mod,  MSA *msa, scale_type stype, 
                              mode_type mode, int second_derivs) {
  FeatFitData *retval = smalloc(sizeof(FeatFitData));
  retval->feat = NULL;
  retval->cdata = col_init_fit_data(mod, msa, stype, mode, second_derivs);
  return retval;
}

/* Free metadata and memory for fitting scale factors */
void ff_free_fit_data(FeatFitData *d) {
  col_free_fit_data(d->cdata);
}

/* Identify branches wrt which a given feature is uninformative,
   in the sense that all leaves beneath these branches having only missing
   data.  Will set (preallocated) array has_data[i] = I(branch above
   node i is informative).  Will also set *nspec equal to number of
   leaves that have data. */
void ff_find_missing_branches(TreeModel *mod, MSA *msa, GFF_Feature *feat, 
                              int *has_data, int *nspec) {
  int i, j;
  List *traversal = tr_postorder(mod->tree);
  *nspec = 0;
  for (i = 0; i < lst_size(traversal); i++) {
    TreeNode *n = lst_get_ptr(traversal, i);
    if (!((n->lchild == NULL && n->rchild == NULL) || 
	  (n->lchild != NULL && n->rchild != NULL)))
      die("ERROR ff_find_missing_branches: lchild and rchild should both be NULL or not NULL\n");
    if (n->parent == NULL)      /* root */
      has_data[n->id] = FALSE;
    else if (n->lchild == NULL) {    /* leaf */
      has_data[n->id] = FALSE;       /* initialize to F, set to T if
                                        base in any col in feature */
      for (j = feat->start-1; j < feat->end; j++) {
        if (mod->rate_matrix->
            inv_states[(int)ss_get_char_tuple(msa, msa->ss->tuple_idx[j], 
                                              mod->msa_seq_idx[n->id], 0)] >= 0) {
          has_data[n->id] = TRUE;
          (*nspec)++;        
          break;
        }
      }
    }
    else {                      /* non-root ancestral node */
      if (has_data[n->lchild->id] || has_data[n->rchild->id])
        has_data[n->id] = TRUE;
      else 
        has_data[n->id] = FALSE;
    }
  }
}

/* returns TRUE if at least one column in feature has two or more
   actual bases (not gaps or missing data), otherwise returns FALSE */
int ff_has_data(TreeModel *mod, MSA *msa, GFF_Feature *f) {
  int i, j;
  for (j = f->start-1; j < f->end; j++) {
    int tupleidx = msa->ss->tuple_idx[j];
    int nbases = 0;
    for (i = 0; i < msa->nseqs && nbases < 2; i++) {
      int state = mod->rate_matrix->
        inv_states[(int)ss_get_char_tuple(msa, tupleidx, i, 0)];
      if (state >= 0) 
        nbases++;
    }
    if (nbases >= 2) return TRUE;
  }
  return FALSE;
}

/* returns TRUE if at least one column in feature meets the criteria
   of col_has_data_sub, otherwise returns FALSE */
int ff_has_data_sub(TreeModel *mod, MSA *msa, GFF_Feature *f, List *inside,
                    List *outside) {
  int j, retval = FALSE;
  for (j = f->start-1; j < f->end && !retval; j++) {
    if (col_has_data_sub(mod, msa, msa->ss->tuple_idx[j], inside, outside))
      retval = TRUE;
  }
  return retval;
}
