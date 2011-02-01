/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* Functions that output data computed by phyloP */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <getopt.h>
#include <misc.h>
#include <msa.h>
#include <sufficient_stats.h>
#include <subst_distrib.h>
#include <prob_vector.h>
#include <prob_matrix.h>
#include <phylo_p_print.h>
#include <list_of_lists.h>

void print_quantiles(FILE *outfile, Vector *distrib, ListOfLists *result) {
  int *quantiles = pv_quantiles(distrib);
  int i;
  if (outfile != NULL) {
    for (i = 0; i <= 100; i++) 
      fprintf(outfile, "%.2f\t%d\n", 1.0*i/100, quantiles[i]);
  }
  if (result != NULL) {
    double *tmpd = smalloc(101*sizeof(double));
    ListOfLists *group = lol_new(2);
    for (i=0; i<=100; i++) tmpd[i]=1.0*i/100;
    lol_push_dbl(group, tmpd, 101, "quantile");
    sfree(tmpd);
    lol_push_int(group, quantiles, 101, "nsub");
    lol_set_class(group, "data.frame");
    lol_push_lol(result, group, "nsub.quantile");
  }
  sfree(quantiles);
}

void print_prior_only(FILE *outfile, int nsites, char *mod_fname, 
		      Vector *prior_distrib,
		      ListOfLists *result) {
  int i, prior_min, prior_max;
  double prior_mean, prior_var;
  pv_stats(prior_distrib, &prior_mean, &prior_var);
  pv_confidence_interval(prior_distrib, 0.95, &prior_min, &prior_max);
  if (outfile != NULL) {
    fprintf(outfile, "#Let n be no. substitutions in %d sites given ", nsites);
    if (mod_fname != NULL)
      fprintf(outfile, "'%s'.\n", 
	      mod_fname);
    else fprintf(outfile, "the model\n");
    fprintf(outfile, "#E[n] = %.3f; Var[n] = %.3f; 95%% c.i. = [%d, %d]\n", 
	    prior_mean, prior_var, prior_min, prior_max);
    fprintf(outfile, "#n p(n)\n");
    for (i = 0; i < prior_distrib->size; i++)
      fprintf(outfile, "%d\t%f\n", i, prior_distrib->data[i]);
  }
  if (result != NULL) {
    ListOfLists *group = lol_new(2);
    int *tmpint = smalloc(prior_distrib->size*sizeof(int));
    for (i=0; i<prior_distrib->size; i++) tmpint[i]=i;
    lol_push_int(group, tmpint, prior_distrib->size, "nsub");
    sfree(tmpint);
    lol_push_dbl(group, prior_distrib->data, prior_distrib->size, "prior.distrib");
    lol_set_class(group, "data.frame");
    lol_push_lol(result, group, "prior");
  }
}


void print_post_only(FILE *outfile, char *mod_fname, char *msa_fname, 
		     Vector *post_distrib, double ci, double scale,
		     ListOfLists *result) {
  int i, post_min, post_max;
  double post_mean, post_var;
  if (ci == -1) ci = 0.95;      /* for purposes of stats */
  pv_stats(post_distrib, &post_mean, &post_var);
  pv_confidence_interval(post_distrib, ci, &post_min, &post_max);
  if (outfile != NULL) {
    fprintf(outfile, "#Let n be no. substitutions given ");
    if (mod_fname != NULL && msa_fname != NULL) 
      fprintf(outfile, "'%s' and '%s'.\n", 
	      mod_fname, msa_fname);
    else fprintf(outfile, "the model and alignment\n");
    fprintf(outfile, "#E[n] = %.3f; Var[n] = %.3f; %.1f%% c.i. = [%d, %d]\n", 
	    post_mean, post_var, ci*100, post_min, post_max);
    if (scale != -1)
      fprintf(outfile, "#estimated scale factor: %f\n", scale);
    fprintf(outfile, "#n p(n)\n");
    for (i = 0; i < post_distrib->size; i++)
      fprintf(outfile, "%d\t%f\n", i, post_distrib->data[i]);
  }
  if (result != NULL) {
    ListOfLists *group = lol_new(2);
    int *tmpint = smalloc(post_distrib->size*sizeof(int));
    for (i=0; i<post_distrib->size; i++)
      tmpint[i] = i;
    lol_push_int(group, tmpint, post_distrib->size, "nsub");
    lol_push_dbl(group, post_distrib->data, post_distrib->size, "post.distrib");
    lol_set_class(group, "data.frame");
    if (scale != -1) {
      ListOfLists *tmp = lol_new(2);
      lol_push_dbl(tmp, &scale, 1, "scale");
      lol_push_lol(tmp, group, "post.distrib");
      group = tmp;
    }
    lol_push_lol(result, group, "post.distrib");
  }
}

void print_p(FILE *outfile, char *mod_fname, char *msa_fname, 
	     Vector *prior_distrib, double post_mean, 
	     double post_var, double ci, double scale,
	     ListOfLists *result) {
  double prior_mean, prior_var, post_min, post_max, pcons, pacc;
  int prior_min, prior_max;

  pv_stats(prior_distrib, &prior_mean, &prior_var);
  pv_confidence_interval(prior_distrib, 0.95, &prior_min, &prior_max);

  if (ci != -1)
    /* avoid computing convolution; assume normality based on central
       limit theorem */
    norm_confidence_interval(post_mean, sqrt(post_var), ci, &post_min, &post_max);
  else 
    post_min = post_max = post_mean;

  post_min = floor(post_min); post_max = ceil(post_max);
  pcons = pv_p_value(prior_distrib, post_max, LOWER);
  pacc = pv_p_value(prior_distrib, post_min, UPPER);

  if (outfile != NULL) {
    fprintf(outfile, "\n*****\nP-values for number of substitutions observed in ");
    if (msa_fname != NULL && mod_fname != NULL)
      fprintf(outfile, "'%s' given '%s'\n*****\n\n", 
	    msa_fname, mod_fname);
    else fprintf(outfile, "the alignment given the model\n*****\n\n");
    fprintf(outfile, "p-value of conservation: %e\n", pcons);
    fprintf(outfile, "p-value of acceleration: %e\n\n", pacc);
    fprintf(outfile, 
	    "null distrib: mean = %f, var = %f, 95%% c.i. = [%d, %d]\n", 
	    prior_mean, prior_var, prior_min, prior_max);
    fprintf(outfile, "posterior distrib: mean = %f, var = %f", 
	    post_mean, post_var);
    if (ci != -1)
      fprintf(outfile, ", %.1f%% c.i. = [%.0f, %.0f]", ci*100, post_min, post_max);
    fprintf(outfile, "\n\n");
    if (scale != -1)
      fprintf(outfile, "estimated scale factor: %f\n\n", scale);
  }
  if (result != NULL) {
    ListOfLists *group = lol_new(12);
    lol_push_dbl(group, &pcons, 1, "pval.cons");
    lol_push_dbl(group, &pacc, 1, "pval.acc");
    lol_push_dbl(group, &prior_mean, 1, "prior.mean");
    lol_push_dbl(group, &prior_var, 1, "prior.var");
    lol_push_int(group, &prior_min, 1, "prior.min");
    lol_push_int(group, &prior_max, 1, "prior.max");
    lol_push_dbl(group, &post_mean, 1, "posterior.mean");
    lol_push_dbl(group, &post_var, 1, "posterior.var");
    if (ci != -1) {
      lol_push_dbl(group, &ci, 1, "confidence.int");
      lol_push_dbl(group, &post_min, 1, "posterior.ci.min");
      lol_push_dbl(group, &post_max, 1, "posterior.ci.max");
    }
    if (scale != -1) 
      lol_push_dbl(group, &scale, 1, "scale");
    lol_push_lol(result, group, "distribution");
  }
}

void print_prior_only_joint(FILE *outfile, char *node_name, int nsites, 
			    char *mod_fname, Matrix *prior_distrib,
			    ListOfLists *result) {
  int i, j, min_sup, max_sup, min_sub, max_sub;
  double mean_sup, var_sup, mean_sub, var_sub;
  Vector *marg_sub = vec_new(prior_distrib->nrows),
    *marg_sup = vec_new(prior_distrib->ncols);

  /* get marginal distributions for stats */
  vec_zero(marg_sup); vec_zero(marg_sub);
  for (i = 0; i < prior_distrib->nrows; i++) {
    for (j = 0; j < prior_distrib->ncols; j++) {
      marg_sub->data[i] += prior_distrib->data[i][j];
      marg_sup->data[j] += prior_distrib->data[i][j];
    }
  }

  /* get marginal stats */
  pv_stats(marg_sup, &mean_sup, &var_sup);
  pv_stats(marg_sub, &mean_sub, &var_sub);
  pv_confidence_interval(marg_sup, 0.95, &min_sup, &max_sup);
  pv_confidence_interval(marg_sub, 0.95, &min_sub, &max_sub);

  if (outfile != NULL) {
    fprintf(outfile, "#Let n1 be no. substitutions in supertree above '%s' (excluding leading branch) over %d site(s) given '%s'.\n", 
	    node_name, nsites, mod_fname);
    fprintf(outfile, "#Let n2 be no. substitutions in subtree beneath '%s' (including leading branch) over %d site(s) given '%s'.\n", 
	    node_name, nsites, mod_fname);
    fprintf(outfile, "#E[n1] = %.3f; Var[n1] = %.3f; 95%% c.i. = [%d, %d]\n", 
	    mean_sup, var_sup, min_sup, max_sup);
    fprintf(outfile, "#E[n2] = %.3f; Var[n2] = %.3f; 95%% c.i. = [%d, %d]\n", 
	    mean_sub, var_sub, min_sub, max_sub);
    fprintf(outfile, "\n#element at row n1 and col n2 in table below is p(n1, n2)\n");
    for (i = 0; i < prior_distrib->ncols; i++) 
      for (j = 0; j < prior_distrib->nrows; j++) 
	fprintf(outfile, "%f%c", prior_distrib->data[j][i], 
		j == prior_distrib->nrows - 1 ? '\n' : '\t');
  }
  if (result != NULL) {
    ListOfLists *group = lol_new(3);
    
    lol_push_int(group, &nsites, 1, "nsite");
    lol_push_charvec(group, &node_name, 1, "subtree.node");
    lol_push_matrix(group, prior_distrib, "joint.distrib");
    lol_push_lol(result, group, "joint.distrib");
  }

  vec_free(marg_sup);
  vec_free(marg_sub);
}

void print_post_only_joint(FILE *outfile, char *node_name, char *mod_fname, 
                           char *msa_fname, Matrix *post_distrib, 
                           double ci, double scale, double sub_scale,
			   ListOfLists *result) {
  int i, j, min_sup, max_sup, min_sub, max_sub;
  double mean_sup, var_sup, mean_sub, var_sub;
  Vector *marg_sup = vec_new(post_distrib->ncols),
    *marg_sub = vec_new(post_distrib->nrows);

  if (ci == -1) ci= 0.95;       /* for purposes of stats */

  /* get marginal distributions for stats */
  vec_zero(marg_sup); vec_zero(marg_sub);
  for (i = 0; i < post_distrib->nrows; i++) {
    checkInterruptN(i, 100);
    for (j = 0; j < post_distrib->ncols; j++) {
      marg_sub->data[i] += post_distrib->data[i][j];
      marg_sup->data[j] += post_distrib->data[i][j];
    }
  }

  /* get marginal stats */
  pv_stats(marg_sup, &mean_sup, &var_sup);
  pv_stats(marg_sub, &mean_sub, &var_sub);
  pv_confidence_interval(marg_sup, ci, &min_sup, &max_sup);
  pv_confidence_interval(marg_sub, ci, &min_sub, &max_sub);

  if (outfile != NULL) {
    fprintf(outfile, "#Let n1 be no. substitutions in supertree above '%s' (excluding leading branch) given ", node_name);
    if (mod_fname != NULL && msa_fname != NULL)
      fprintf(outfile, "'%s' and '%s'.\n",
	     mod_fname, msa_fname);
    else fprintf(outfile, "the model and alignment.\n");
    fprintf(outfile, "#Let n2 be no. substitutions in subtree beneath '%s' (including leading branch) given ", mod_fname);
    if (mod_fname != NULL && msa_fname != NULL) 
      fprintf(outfile, "'%s' and '%s'.\n", 
	      mod_fname, msa_fname);
    else fprintf(outfile, "the model and alignment\n");
    fprintf(outfile, "#E[n1] = %.3f; Var[n1] = %.3f; %.1f%% c.i. = [%d, %d]\n", 
	    mean_sup, var_sup, ci*100, min_sup, max_sup);
    fprintf(outfile, "#E[n2] = %.3f; Var[n2] = %.3f; %.1f%% c.i. = [%d, %d]\n", 
	    mean_sub, var_sub, ci*100, min_sub, max_sub);
    if (scale != -1)
      fprintf(outfile, "#estimated scale factors: %f [tree], %f [subtree]\n", scale, sub_scale);
    fprintf(outfile, "\n#element at row n1 and col n2 in table below is p(n1, n2)\n");
    for (i = 0; i < post_distrib->ncols; i++) {
      checkInterruptN(i, 100);
      for (j = 0; j < post_distrib->nrows; j++) 
	fprintf(outfile, "%f%c", post_distrib->data[j][i], 
		j == post_distrib->nrows - 1 ? '\n' : '\t');
    }
  }
  if (result != NULL) {
    ListOfLists *group = lol_new(2);
    
    lol_push_charvec(group, &node_name, 1, "subtree.node");
    lol_push_matrix(group, post_distrib, "joint.distrib");
    lol_push_lol(result, group, "joint.distrib");
  }

  vec_free(marg_sup);
  vec_free(marg_sub);
}

void print_p_joint(FILE *outfile, char *node_name, char *mod_fname, 
		   char *msa_fname, 
                   double ci, Matrix *prior_joint, 
                   double post_mean, double post_var, 
                   double post_mean_sup, double post_var_sup, 
                   double post_mean_sub, double post_var_sub,
                   double scale, double sub_scale,
		   ListOfLists *result) {

  double post_min_tot, post_max_tot, post_min_sup, post_max_sup, 
    post_min_sub, post_max_sub, cond_cons_p_sub, cond_anti_cons_p_sub, 
    prior_mean_sup, prior_var_sup, prior_mean_sub, prior_var_sub, 
    cons_p_sup, anti_cons_p_sup, cons_p_sub, anti_cons_p_sub;
  int prior_min_sup, prior_max_sup, prior_min_sub, prior_max_sub;
  Vector *prior_marg_sup, *prior_marg_sub, *cond;

  if (ci != -1) {
    norm_confidence_interval(post_mean, sqrt(post_var), ci, &post_min_tot, 
                             &post_max_tot);
    norm_confidence_interval(post_mean_sup, sqrt(post_var_sup), ci, 
                             &post_min_sup, &post_max_sup);
    norm_confidence_interval(post_mean_sub, sqrt(post_var_sub), ci, 
                             &post_min_sub, &post_max_sub);
  }
  else {
    post_min_tot = post_max_tot = post_mean;
    post_min_sup = post_max_sup = post_mean_sup;
    post_min_sub = post_max_sub = post_mean_sub;
  }

  post_min_tot = floor(post_min_tot); post_max_tot = ceil(post_max_tot);
  post_min_sup = floor(post_min_sup); post_max_sup = ceil(post_max_sup);
  post_min_sub = floor(post_min_sub); post_max_sub = ceil(post_max_sub);

  /* Conditional p-values of conservation.  To be conservative, base
     these on the largest reasonable estimate of the number of
     subst. in the subtree and the smallest reasonable estimate of the
     number of subst. in the whole tree*/

  cond = pm_x_given_tot(prior_joint, post_min_tot);
  cond_cons_p_sub = pv_p_value(cond, post_max_sub, LOWER);
  vec_free(cond);

  cond = pm_x_given_tot(prior_joint, post_max_tot);
  cond_anti_cons_p_sub = pv_p_value(cond, post_min_sub, UPPER);
  vec_free(cond);

  /* marginals of prior and stats */
  prior_marg_sup = pm_marg_y(prior_joint);
  prior_marg_sub = pm_marg_x(prior_joint);
  pv_stats(prior_marg_sup, &prior_mean_sup, &prior_var_sup);
  pv_stats(prior_marg_sub, &prior_mean_sub, &prior_var_sub);
  pv_confidence_interval(prior_marg_sup, 0.95, &prior_min_sup, &prior_max_sup);
  pv_confidence_interval(prior_marg_sub, 0.95, &prior_min_sub, &prior_max_sub);

  /* marginal p-values */
  cons_p_sup = pv_p_value(prior_marg_sup, post_max_sup, LOWER);
  anti_cons_p_sup = pv_p_value(prior_marg_sup, post_min_sup, UPPER);
  cons_p_sub = pv_p_value(prior_marg_sub, post_max_sub, LOWER);
  anti_cons_p_sub = pv_p_value(prior_marg_sub, post_min_sub, UPPER);

  vec_free(prior_marg_sup);
  vec_free(prior_marg_sub);

  if (outfile != NULL) {
    fprintf(outfile, "\n*****\nP-values for number of substitutions observed in ");
    if (msa_fname != NULL && mod_fname != NULL)
      fprintf(outfile, "'%s' given '%s',\n", 
	      msa_fname, mod_fname);
    else fprintf(outfile, "the alignment given the model\n");
    fprintf (outfile, "considering subtree/supertree beneath/above node '%s'\n*****\n\n", node_name);
    
    fprintf(outfile, "p-value of conservation in subtree: %e\n", cons_p_sub);
    fprintf(outfile, "p-value of acceleration in subtree: %e\n\n", anti_cons_p_sub);
    
    fprintf(outfile, "p-value of conservation in supertree: %e\n", cons_p_sup);
    fprintf(outfile, "p-value of acceleration in supertree: %e\n\n", anti_cons_p_sup);
    
    fprintf(outfile, "p-value of conservation in subtree given total: %e\n", cond_cons_p_sub);
    fprintf(outfile, "p-value of acceleration in subtree given total: %e\n\n", cond_anti_cons_p_sub);
    
    fprintf(outfile, "null distrib in subtree: mean = %f, var = %f, 95%% c.i. = [%d, %d]\n", 
	    prior_mean_sub, prior_var_sub, prior_min_sub, prior_max_sub);
    fprintf(outfile, "posterior distrib in subtree: mean = %f, var = %f", 
	    post_mean_sub, post_var_sub);
    if (ci != -1)
      fprintf(outfile, ", %.1f%% c.i. = [%.0f, %.0f]", ci*100, post_min_sub, post_max_sub);
    fprintf(outfile, "\n\nnull distrib in supertree: mean = %f, var = %f, 95%% c.i. = [%d, %d]\n", 
	    prior_mean_sup, prior_var_sup, prior_min_sup, prior_max_sup);
    fprintf(outfile, "posterior distrib in supertree: mean = %f, var = %f",
	    post_mean_sup, post_var_sup);
    if (ci != -1)
      fprintf(outfile, ", %.1f%% c.i. = [%.0f, %.0f]", ci*100, post_min_sup, post_max_sup);
    fprintf(outfile, "\n\n");
    if (scale != -1)
      fprintf(outfile, "estimated scale factors: %f [tree], %f [subtree]\n\n", scale,
	      sub_scale);
  }
  if (result != NULL) {
    ListOfLists *group = lol_new(30);
    lol_push_charvec(group, &node_name, 1, "subtree.node");
    lol_push_dbl(group, &cons_p_sub, 1, "pval.cons.subtree");
    lol_push_dbl(group, &anti_cons_p_sub, 1, "pval.acc.subtree");
    lol_push_dbl(group, &cons_p_sup, 1, "pval.cons.supertree");
    lol_push_dbl(group, &anti_cons_p_sup, 1, "pval.acc.supertree");
    lol_push_dbl(group, &cond_cons_p_sub, 1, "pval.cons.subtree.given.total");
    lol_push_dbl(group, &cond_anti_cons_p_sub, 1, "pval.acc.subtree.given.total");
    lol_push_dbl(group, &prior_mean_sub, 1, "prior.subtree.mean");
    lol_push_dbl(group, &prior_var_sub, 1, "prior.subtree.var");
    lol_push_int(group, &prior_min_sub, 1, "prior.subtree.ci95.min");
    lol_push_int(group, &prior_max_sub, 1, "prior.subtree.ci95.max");
    lol_push_dbl(group, &post_mean_sub, 1, "post.subtree.mean");
    lol_push_dbl(group, &post_var_sub, 1, "post.subtree.var");
    if (ci != -1) {
      lol_push_dbl(group, &ci, 1, "post.conf.int");
      lol_push_dbl(group, &post_min_sub, 1, "post.subtree.conf.min");
      lol_push_dbl(group, &post_max_sub, 1, "post.subtree.conf.max");
      lol_push_dbl(group, &post_min_sup, 1, "post.supertree.conf.min");
      lol_push_dbl(group, &post_max_sup, 1, "post.supertree.conf.max");
    }
    lol_push_dbl(group, &prior_mean_sup, 1, "prior.supertree.mean");
    lol_push_dbl(group, &prior_var_sup, 1, "prior.supertree.var");
    lol_push_int(group, &prior_min_sup, 1, "prior.supertree.ci95.min");
    lol_push_int(group, &prior_max_sup, 1, "prior.supertree.ci95.max");
    lol_push_dbl(group, &post_mean_sup, 1, "post.supertree.mean");
    lol_push_dbl(group, &post_var_sup, 1, "post.supertree.var");
    if (scale != -1) {
      lol_push_dbl(group, &scale, 1, "scale");
      lol_push_dbl(group, &sub_scale, 1, "subtree.scale");
    }
    lol_push_lol(result, group, "distrib.stats");
  }
}

/* Features output for SPH without subtree */
void print_feats_sph(FILE *outfile, p_value_stats *stats, GFF_Set *feats, 
                     mode_type mode, double epsilon, int output_gff,
		     ListOfLists *result) {
  int i;
  double *pvals = smalloc(lst_size(feats->features) * sizeof(double)),
    *post_means = NULL, *post_vars = NULL, *prior_means = NULL, 
    *prior_vars = NULL;

  if (result != NULL || !output_gff) {
    post_means = smalloc(lst_size(feats->features) * sizeof(double));
    post_vars = smalloc(lst_size(feats->features) * sizeof(double));
    prior_means = smalloc(lst_size(feats->features) * sizeof(double));
    prior_vars = smalloc(lst_size(feats->features) * sizeof(double));
  }
  for (i = 0; i < lst_size(feats->features); i++) {
    checkInterruptN(i, 100);
    if (result != NULL || !output_gff) {
      post_means[i] = stats[i].post_mean;
      post_vars[i] = stats[i].post_var;
      prior_means[i] = stats[i].prior_mean;
      prior_vars[i] = stats[i].prior_var;
    }

    if (mode == CON)
      pvals[i] = stats[i].p_cons;
    else if (mode == ACC)
      pvals[i] = stats[i].p_anti_cons;
    else if (mode == NNEUT)
      pvals[i] = 2 * (min(stats[i].p_cons, stats[i].p_anti_cons));
    else {
      if (mode != CONACC)
	die("ERROR print_feats_sph: unknown mode\n");
      if (stats[i].p_cons < stats[i].p_anti_cons)
        pvals[i] = stats[i].p_cons;
      else 
        pvals[i] = -stats[i].p_anti_cons;
    }

    if (pvals[i] == 0) {
      if (mode == ACC)
        pvals[i] = epsilon;
      else if (mode == CONACC)
        pvals[i] = -epsilon;
      else if (mode == NNEUT)
        pvals[i] = 2*epsilon;
      /* in these cases, reset pvals of zero to epsilon (or
         2*epsilon), because off scale of finite representation of
         distrib */
    }
  }
  if (output_gff && outfile != NULL) 
    print_gff_scores(outfile, feats, pvals, TRUE);
  if (result != NULL || !output_gff) 
    print_feats_generic(output_gff ? NULL : outfile, 
			"prior_mean\tprior_var\tpost_mean\tpost_var\tpval",
                        feats, NULL, result, FALSE, TRUE, 5, 
			"prior.mean", prior_means, "prior.var", prior_vars, 
			"post.mean", post_means, 
                        "post.var", post_vars, "pval", pvals);
  if (result != NULL || !output_gff) {
    sfree(post_means);
    sfree(post_vars);
    sfree(prior_means);
    sfree(prior_vars);
  }
  sfree(pvals);
}

/* Features output for SPH with subtree */
void print_feats_sph_subtree(FILE *outfile, p_value_joint_stats *stats, 
			     GFF_Set *feats, 
                             mode_type mode, double epsilon, int output_gff,
			     ListOfLists *result) {
  int i;
  double *pvals = smalloc(lst_size(feats->features) * sizeof(double)),
    *post_means_sub = NULL, *post_vars_sub = NULL, 
    *post_means_sup = NULL, *post_vars_sup = NULL, 
    *prior_means_sub = NULL, *prior_vars_sub = NULL,
    *prior_means_sup = NULL, *prior_vars_sup = NULL;

  if (!output_gff) {
    post_means_sub = smalloc(lst_size(feats->features) * sizeof(double));
    post_vars_sub = smalloc(lst_size(feats->features) * sizeof(double));
    post_means_sup = smalloc(lst_size(feats->features) * sizeof(double));
    post_vars_sup = smalloc(lst_size(feats->features) * sizeof(double));
    prior_means_sub = smalloc(lst_size(feats->features) * sizeof(double));
    prior_vars_sub = smalloc(lst_size(feats->features) * sizeof(double));
    prior_means_sup = smalloc(lst_size(feats->features) * sizeof(double));
    prior_vars_sup = smalloc(lst_size(feats->features) * sizeof(double));
  }
  for (i = 0; i < lst_size(feats->features); i++) {
    checkInterruptN(i, 100);
    if (!output_gff) {
      post_means_sub[i] = stats[i].post_mean_left;
      post_vars_sub[i] = stats[i].post_var_left;
      post_means_sup[i] = stats[i].post_mean_right;
      post_vars_sup[i] = stats[i].post_var_right;
      prior_means_sub[i] = stats[i].prior_mean_left;
      prior_vars_sub[i] = stats[i].prior_var_left;
      prior_means_sup[i] = stats[i].prior_mean_right;
      prior_vars_sup[i] = stats[i].prior_var_right;
    }

    if (mode == CON)
      pvals[i] = stats[i].p_cons_left;
    else if (mode == ACC)
      pvals[i] = stats[i].p_anti_cons_left;
    else if (mode == NNEUT)
      pvals[i] = 2 * (min(stats[i].p_cons_left, stats[i].p_anti_cons_left));
    else {
      if (mode != CONACC)
	die("ERROR print_feats_sph_subtree: unknown mode\n");
      if (stats[i].p_cons_left < stats[i].p_anti_cons_left)
        pvals[i] = stats[i].p_cons_left;
      else 
        pvals[i] = -stats[i].p_anti_cons_left;
    }

    if (pvals[i] == 0) {
      if (mode == ACC)
        pvals[i] = epsilon;
      else if (mode == CONACC)
        pvals[i] = -epsilon;
      else if (mode == NNEUT)
        pvals[i] = 2*epsilon;
      /* in these cases, reset pvals of zero to epsilon (or
         2*epsilon), because off scale of finite representation of
         distrib */
    }
  }
  if (output_gff) 
    print_gff_scores(outfile, feats, pvals, TRUE);
  else 
    print_feats_generic(outfile, "prior_mean_sub\tprior_var_sub\tprior_mean_sup\tprior_var_sup\tpost_mean_sub\tpost_var_sub\tpost_mean_sup\tpost_var_sup\t\tpval",
                        feats, NULL, result, FALSE, TRUE, 9, 
			"prior.mean.sub", prior_means_sub, 
			"prior.var.sub", prior_vars_sub, 
                        "prior.mean.sup", prior_means_sup, 
			"prior.var.sup", prior_vars_sup, 
			"post.mean.sub", post_means_sub, 
                        "post.var.sub", post_vars_sub, 
			"post.mean.sup", post_means_sup, 
			"post.var.sup", post_vars_sup, 
			"pval", pvals);
  if (!output_gff) {
    sfree(post_means_sub);
    sfree(post_vars_sub);
    sfree(post_means_sup);
    sfree(post_vars_sup);
    sfree(prior_means_sub);
    sfree(prior_vars_sub);
    sfree(prior_means_sup);
    sfree(prior_vars_sup);
  }
  sfree(pvals);
}


void print_wig(FILE *outfile, MSA *msa, double *vals, char *chrom, 
	       int refidx, int log_trans, ListOfLists *result) {
  int last, j, k;
  double val;
  List *posList=NULL, *scoreList=NULL;

  if (result != NULL) {
    posList = lst_new_int(msa->length);
    scoreList = lst_new_dbl(msa->length);
  }

  last = -INFTY;
  if (!(refidx >= 0 && refidx <= msa->nseqs))
    die("ERROR print_wig: bad refidx (%i)\n", refidx);
  for (j = 0, k = 0; j < msa->length; j++) {
    checkInterruptN(j, 1000);
    if (refidx == 0 || msa_get_char(msa, refidx-1, j) != GAP_CHAR) {
      if (refidx == 0 || !msa_missing_col(msa, refidx, j)) {
        if (k > last + 1 && outfile != NULL) 
          fprintf(outfile, "fixedStep chrom=%s start=%d step=1\n", chrom, 
		  k + msa->idx_offset + 1);
        val = vals[msa->ss->tuple_idx[j]];
        if (log_trans) {
          int sign = 1;
          if (val < 0) {
            val = -val;
            sign = -1;          /* propagate negative sign through */
          }
          val = fabs(-log10(val)) * sign; /* fabs prevents -0 for val == 1 */
        }
        if (outfile != NULL) fprintf(outfile, "%.3f\n", val);
	if (result != NULL) {
	  lst_push_int(posList, k + msa->idx_offset + 1);
	  lst_push_dbl(scoreList, val);
	}
        last = k;
      }
      k++;
    }
  }
  if (result != NULL) {
    ListOfLists *group = lol_new(2);
    lol_push(group, posList, "coord", INT_LIST);
    lol_push(group, scoreList, "score", DBL_LIST);
    lol_set_class(group, "data.frame");
    lol_push_lol(result, group, "wig");
  }
}


double *log10_pval(double *pval, int len) {
  double *scores = smalloc(len*sizeof(double)), sign;
  int i;

  for (i=0; i < len; i++) {
    if (pval[i] < 0) sign = -1.0;
    else sign=1.0;
    scores[i] = fabs(-log10(sign*pval[i])) * sign;
  }
  return scores;

}
	  


/* Print arbitrary columns of tuple-specific data in wig-like format */
void print_base_by_base(FILE *outfile, char *header, char *chrom, MSA *msa, 
                        char **formatstr, int refidx, ListOfLists *result,
			int log_trans_outfile, int log_trans_results,
			int ncols, ...) {
  int last, j, k, tup, col;
  va_list ap;
  double *data[ncols+1];
  List **resultList=NULL;
  char **colname;
  int get_log = (log_trans_outfile && outfile != NULL) || 
    (log_trans_results && result != NULL);
  
  if (!(refidx >= 0 && refidx <= msa->nseqs))
    die("ERROR print_base_by_base: bad refidx (%i)\n", refidx);

  if (result != NULL) {
    int resultLen = ncols + 1 + log_trans_results;
    resultList = smalloc(resultLen*sizeof(List*));
    resultList[0] = lst_new_int(msa->length);
    for (j=1; j<resultLen; j++) {
      resultList[j] = lst_new_dbl(msa->length);
    }
  }

  last = -INFTY;
  if (header != NULL && outfile != NULL) {
    fprintf(outfile, "%s", header);
    if (log_trans_outfile) fprintf(outfile, " score\n");
    else fprintf(outfile, "\n");
  }

  va_start(ap, ncols);
  colname = smalloc((ncols+1)*sizeof(char*));
  for (col = 0; col < ncols; col++) {
    colname[col] = va_arg(ap, char*);
    data[col] = va_arg(ap, double*);
  }
  if (get_log) {
    if (!(strcmp(colname[col-1], "pval")==0))
      die("ERROR print_base_by_base expected last col to be pval, got %s\n",
	  colname[col-1]);
    colname[col] = "score";
    data[col] = log10_pval(data[col-1], msa->ss->ntuples);
  }

  for (j = 0, k = 0; j < msa->length; j++) {
    checkInterruptN(j, 1000);
    if (refidx == 0 || msa_get_char(msa, refidx-1, j) != GAP_CHAR) {
      if (refidx == 0 || !msa_missing_col(msa, refidx, j)) {
        if (k > last + 1 && outfile != NULL) 
          fprintf(outfile, "fixedStep chrom=%s start=%d step=1\n", chrom, 
                 k + msa->idx_offset + 1);
        tup = msa->ss->tuple_idx[j];
	if (outfile != NULL) {
	  for (col = 0; col < ncols; col++) {
	    fprintf(outfile, (formatstr == NULL ? "%.5f" : formatstr[col]), data[col][tup]);
	    if (col <  ncols-1) fputc('\t', outfile);
	  }
	  if (log_trans_outfile) fprintf(outfile, "\t%.5f", data[col][tup]);
	  fputc('\n', outfile);
	}
	if (result != NULL) {
	  lst_push_int(resultList[0], k + msa->idx_offset + 1);
	  for (col=0; col < ncols; col++) 
	    lst_push_dbl(resultList[col+1], data[col][tup]);
	  if (log_trans_results)
	    lst_push_dbl(resultList[col+1], data[col][tup]);
	}
        last = k;
      }
      k++;
    }
  }
  va_end(ap);

  if (result != NULL) {
    ListOfLists *group = lol_new(ncols+1+log_trans_results);
    lol_push(group, resultList[0], "coord", INT_LIST);
    for (col=1; col<=ncols; col++) { 
      lol_push(group, resultList[col], colname[col-1], DBL_LIST);
    }
    if (log_trans_results) 
      lol_push(group, resultList[col], colname[col-1], DBL_LIST);
    lol_set_class(group, "data.frame");
    lol_push_lol(result, group, "baseByBase");
    sfree(resultList);
  }
  if (get_log) sfree(data[ncols]);
  sfree(colname);
}


/* Print a list of features and artibrary associated statistics */
void print_feats_generic(FILE *outfile, char *header, GFF_Set *gff, 
			 char **formatstr, ListOfLists *result, 
			 int log_trans_outfile, int log_trans_results, 
			 int ncols, ...) {
  int i, col;
  String *name;
  va_list ap;
  double *data[ncols+1];
  Regex *tag_val_re = str_re_new("[[:alnum:]_.]+[[:space:]]+(\"[^\"]*\"|[^[:space:]]+)");
  List *l = lst_new_ptr(2);
  char **colname;
  List **resultList=NULL;
  int get_log = (log_trans_outfile && outfile != NULL) || 
    (log_trans_results && result != NULL);

  if (outfile != NULL) {
    fprintf(outfile, "#chr\tstart\tend\tname");
    if (header != NULL) {
      fprintf(outfile, "\t%s", header);
      if (log_trans_outfile) fprintf(outfile, "\tscore\n");
      else fprintf(outfile, "\n");
    }
    else 
      fprintf(outfile, "\n");
  }
  if (result != NULL) {
    resultList = smalloc(5*sizeof(List*));
    resultList[0] = lst_new_ptr(lst_size(gff->features));
    resultList[1] = lst_new_int(lst_size(gff->features));
    resultList[2] = lst_new_int(lst_size(gff->features));
    resultList[3] = lst_new_ptr(lst_size(gff->features));
    resultList[4] = lst_new_ptr(lst_size(gff->features));
  }

  colname = smalloc((ncols+1)*sizeof(char*));
  va_start(ap, ncols);
  for (col = 0; col < ncols; col++) {
    colname[col] = va_arg(ap, char*);
    data[col] = va_arg(ap, double*);
  }
  if (get_log) {
    if (!(strcmp(colname[col-1], "pval")==0))
      die("ERROR print_feats_generic: expected last col to be pval, got %s\n",
	  colname[col-1]);
    colname[col] = "score";
    data[col] = log10_pval(data[col-1], lst_size(gff->features));
  }

  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *f = lst_get_ptr(gff->features, i);
    checkInterruptN(i, 100);

    /* try to extract feature name from attribute field */
    lst_clear(l);
    if (f->attribute->length > 0 && 
        str_re_match(f->attribute, tag_val_re, l, 1) >= 0) {
      name = lst_get_ptr(l, 1);
      str_remove_quotes(name);
    } else name=NULL;

    if (outfile != NULL) {
      fprintf(outfile, "%s\t%d\t%d\t%s\t", f->seqname->chars, f->start-1, f->end, 
	      name == NULL ? "." : name->chars);
      
      for (col = 0; col < ncols; col++) {
	fprintf(outfile, (formatstr == NULL ? "%.5f" : formatstr[col]), data[col][i]);
	if (col < ncols-1) fputc('\t', outfile);
      }
      if (log_trans_outfile) fprintf(outfile, "\t%.5f\n", data[col][i]);
      else fputc('\n', outfile);
	  
    }
    if (result != NULL) {
      char *tempstr;
      tempstr = copy_charstr(f->seqname->chars);
      lst_push_ptr(resultList[0], tempstr);
      lst_push_int(resultList[1], f->start);
      lst_push_int(resultList[2], f->end);
      lst_push_ptr(resultList[3], copy_charstr(f->feature == NULL ? "." : 
					       f->feature->chars));
      tempstr = copy_charstr(name == NULL ? "." : name->chars);
      lst_push_ptr(resultList[4], tempstr);
    }
    lst_free_strings(l);
  }
  
  if (result != NULL) {
    ListOfLists *group = lol_new(5+ncols+log_trans_results);
    lol_push(group, resultList[0], "chr", CHAR_LIST);
    lol_push(group, resultList[1], "start", INT_LIST);
    lol_push(group, resultList[2], "end", INT_LIST);
    lol_push(group, resultList[3], "feature", CHAR_LIST);
    lol_push(group, resultList[4], "name", CHAR_LIST);
    sfree(resultList);
    for (col=0; col < ncols; col++) 
      lol_push_dbl(group, data[col], lst_size(gff->features), 
			   colname[col]);
    if (log_trans_results)
      lol_push_dbl(group, data[col], lst_size(gff->features),
		   colname[col]);
    lol_set_class(group, "data.frame");
    lol_push_lol(result, group, "feature.stats");
  }

  va_end(ap);
  lst_free(l);
  str_re_free(tag_val_re);
  if (get_log) sfree(data[ncols]);
  sfree(colname);
}

/* Print GFF to stdout with feature scores defined by vals.  If
   log_trans == TRUE, take log transform (propagating negative
   signs) */
void print_gff_scores(FILE *outfile, GFF_Set *gff, double *vals, int log_trans) {
  int i;
  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *f = lst_get_ptr(gff->features, i);
    checkInterruptN(i, 100);
    f->score = vals[i];
    f->score_is_null = FALSE;
    if (log_trans) {
      int sign = 1;
      if (f->score < 0) {
        f->score = -f->score;
        sign = -1;          /* propagate negative sign through */
      }
      f->score = fabs(-log10(f->score)) * sign; /* fabs prevents -0
                                                   for val == 1 */
    }
  }
  if (outfile != NULL) gff_print_set(outfile, gff);
}
