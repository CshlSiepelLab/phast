/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
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
#include <phyloP.h>

void print_quantiles(Vector *distrib) {
  int *quantiles = pv_quantiles(distrib);
  int i;
  for (i = 0; i <= 100; i++)
    printf("%.2f\t%d\n", 1.0*i/100, quantiles[i]);
  free(quantiles);
}

void print_prior_only(int nsites, char *mod_fname, Vector *prior_distrib) {
  int i, prior_min, prior_max;
  double prior_mean, prior_var;
  pv_stats(prior_distrib, &prior_mean, &prior_var);
  pv_confidence_interval(prior_distrib, 0.95, &prior_min, &prior_max);
  printf("#Let n be no. substitutions in %d sites given '%s'.\n", 
         nsites, mod_fname);
  printf("#E[n] = %.3f; Var[n] = %.3f; 95%% c.i. = [%d, %d]\n", 
         prior_mean, prior_var, prior_min, prior_max);
  printf("#n p(n)\n");
  for (i = 0; i < prior_distrib->size; i++)
    printf("%d\t%f\n", i, prior_distrib->data[i]);
}

void print_post_only(char *mod_fname, char *msa_fname, Vector *post_distrib,
                     double ci, double scale) {
  int i, post_min, post_max;
  double post_mean, post_var;
  if (ci == -1) ci = 0.95;      /* for purposes of stats */
  pv_stats(post_distrib, &post_mean, &post_var);
  pv_confidence_interval(post_distrib, ci, &post_min, &post_max);
  printf("#Let n be no. substitutions given '%s' and '%s'.\n", 
         mod_fname, msa_fname);
  printf("#E[n] = %.3f; Var[n] = %.3f; %.1f%% c.i. = [%d, %d]\n", 
         post_mean, post_var, ci*100, post_min, post_max);
  if (scale != -1)
    printf("#estimated scale factor: %f\n", scale);
  printf("#n p(n)\n");
  for (i = 0; i < post_distrib->size; i++)
    printf("%d\t%f\n", i, post_distrib->data[i]);
}

void print_p(char *mod_fname, char *msa_fname, Vector *prior_distrib,
             double post_mean, double post_var, double ci, double scale) {
  double prior_mean, prior_var, post_min, post_max;
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

  printf("\n*****\nP-values for number of substitutions observed in '%s' given '%s'\n*****\n\n", 
         msa_fname, mod_fname);
  printf("p-value of conservation: %e\n", 
         pv_p_value(prior_distrib, post_max, LOWER));
  printf("p-value of acceleration: %e\n\n",
         pv_p_value(prior_distrib, post_min, UPPER));
  printf("null distrib: mean = %f, var = %f, 95%% c.i. = [%d, %d]\n", 
         prior_mean, prior_var, prior_min, prior_max);
  printf("posterior distrib: mean = %f, var = %f", post_mean, post_var);
  if (ci != -1)
    printf(", %.1f%% c.i. = [%.0f, %.0f]", ci*100, post_min, post_max);
  printf("\n\n");
  if (scale != -1)
    printf("estimated scale factor: %f\n\n", scale);
}

void print_prior_only_joint(char *node_name, int nsites, char *mod_fname, 
                            Matrix *prior_distrib) {
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

  printf("#Let n1 be no. substitutions in supertree above '%s' (excluding leading branch) over %d site(s) given '%s'.\n", 
         node_name, nsites, mod_fname);
  printf("#Let n2 be no. substitutions in subtree beneath '%s' (including leading branch) over %d site(s) given '%s'.\n", 
         node_name, nsites, mod_fname);
  printf("#E[n1] = %.3f; Var[n1] = %.3f; 95%% c.i. = [%d, %d]\n", 
         mean_sup, var_sup, min_sup, max_sup);
  printf("#E[n2] = %.3f; Var[n2] = %.3f; 95%% c.i. = [%d, %d]\n", 
         mean_sub, var_sub, min_sub, max_sub);
  printf("\n#element at row n1 and col n2 in table below is p(n1, n2)\n");
  for (i = 0; i < prior_distrib->ncols; i++) 
    for (j = 0; j < prior_distrib->nrows; j++) 
      printf("%f%c", prior_distrib->data[j][i], 
             j == prior_distrib->nrows - 1 ? '\n' : '\t');

  vec_free(marg_sup);
  vec_free(marg_sub);
}

void print_post_only_joint(char *node_name, char *mod_fname, 
                           char *msa_fname, Matrix *post_distrib, 
                           double ci, double scale, double sub_scale) {
  int i, j, min_sup, max_sup, min_sub, max_sub;
  double mean_sup, var_sup, mean_sub, var_sub;
  Vector *marg_sup = vec_new(post_distrib->ncols),
    *marg_sub = vec_new(post_distrib->nrows);

  if (ci == -1) ci= 0.95;       /* for purposes of stats */

  /* get marginal distributions for stats */
  vec_zero(marg_sup); vec_zero(marg_sub);
  for (i = 0; i < post_distrib->nrows; i++) {
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

  printf("#Let n1 be no. substitutions in supertree above '%s' (excluding leading branch) given '%s' and '%s'.\n", 
         node_name, mod_fname, msa_fname);
  printf("#Let n2 be no. substitutions in subtree beneath '%s' (including leading branch) given '%s' and '%s'.\n", 
         node_name, mod_fname, msa_fname);
  printf("#E[n1] = %.3f; Var[n1] = %.3f; %.1f%% c.i. = [%d, %d]\n", 
         mean_sup, var_sup, ci*100, min_sup, max_sup);
  printf("#E[n2] = %.3f; Var[n2] = %.3f; %.1f%% c.i. = [%d, %d]\n", 
         mean_sub, var_sub, ci*100, min_sub, max_sub);
  if (scale != -1)
    printf("#estimated scale factors: %f [tree], %f [subtree]\n", scale, sub_scale);
  printf("\n#element at row n1 and col n2 in table below is p(n1, n2)\n");
  for (i = 0; i < post_distrib->ncols; i++) 
    for (j = 0; j < post_distrib->nrows; j++) 
      printf("%f%c", post_distrib->data[j][i], 
             j == post_distrib->nrows - 1 ? '\n' : '\t');

  vec_free(marg_sup);
  vec_free(marg_sub);
}

void print_p_joint(char *node_name, char *mod_fname, char *msa_fname, 
                   double ci, Matrix *prior_joint, 
                   double post_mean, double post_var, 
                   double post_mean_sup, double post_var_sup, 
                   double post_mean_sub, double post_var_sub,
                   double scale, double sub_scale) {

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


  printf("\n*****\nP-values for number of substitutions observed in '%s' given '%s',\n", 
         msa_fname, mod_fname);
  printf ("considering subtree/supertree beneath/above node '%s'\n*****\n\n", node_name);

  printf("p-value of conservation in subtree: %e\n", cons_p_sub);
  printf("p-value of acceleration in subtree: %e\n\n", anti_cons_p_sub);

  printf("p-value of conservation in supertree: %e\n", cons_p_sup);
  printf("p-value of acceleration in supertree: %e\n\n", anti_cons_p_sup);

  printf("p-value of conservation in subtree given total: %e\n", cond_cons_p_sub);
  printf("p-value of acceleration in subtree given total: %e\n\n", cond_anti_cons_p_sub);

  printf("null distrib in subtree: mean = %f, var = %f, 95%% c.i. = [%d, %d]\n", 
         prior_mean_sub, prior_var_sub, prior_min_sub, prior_max_sub);
  printf("posterior distrib in subtree: mean = %f, var = %f", 
         post_mean_sub, post_var_sub);
  if (ci != -1)
    printf(", %.1f%% c.i. = [%.0f, %.0f]", ci*100, post_min_sub, post_max_sub);
  printf("\n\nnull distrib in supertree: mean = %f, var = %f, 95%% c.i. = [%d, %d]\n", 
         prior_mean_sup, prior_var_sup, prior_min_sup, prior_max_sup);
  printf("posterior distrib in supertree: mean = %f, var = %f",
         post_mean_sup, post_var_sup);
  if (ci != -1)
    printf(", %.1f%% c.i. = [%.0f, %.0f]", ci*100, post_min_sup, post_max_sup);
  printf("\n\n");
  if (scale != -1)
    printf("estimated scale factors: %f [tree], %f [subtree]\n\n", scale,
           sub_scale);
}

/* Features output for SPH without subtree */
void print_feats_sph(p_value_stats *stats, GFF_Set *feats, 
                     mode_type mode, double epsilon, int output_gff) {
  int i;
  double *pvals = smalloc(lst_size(feats->features) * sizeof(double)),
    *post_means = NULL, *post_vars = NULL, *prior_means = NULL, 
    *prior_vars = NULL;

  if (!output_gff) {
    post_means = smalloc(lst_size(feats->features) * sizeof(double));
    post_vars = smalloc(lst_size(feats->features) * sizeof(double));
    prior_means = smalloc(lst_size(feats->features) * sizeof(double));
    prior_vars = smalloc(lst_size(feats->features) * sizeof(double));
  }
  for (i = 0; i < lst_size(feats->features); i++) {
    if (!output_gff) {
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
      assert(mode == CONACC);
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
  if (output_gff) 
    print_gff_scores(feats, pvals, TRUE);
  else 
    print_feats_generic("prior_mean\tprior_var\tpost_mean\tpost_var\tpval",
                        feats, NULL, 5, prior_means, prior_vars, post_means, 
                        post_vars, pvals);
  if (!output_gff) {
    free(post_means);
    free(post_vars);
    free(prior_means);
    free(prior_vars);
  }
  free(pvals);
}

/* Features output for SPH with subtree */
void print_feats_sph_subtree(p_value_joint_stats *stats, GFF_Set *feats, 
                             mode_type mode, double epsilon, int output_gff) {
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
      assert(mode == CONACC);
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
    print_gff_scores(feats, pvals, TRUE);
  else 
    print_feats_generic("prior_mean_sub\tprior_var_sub\tprior_mean_sup\tprior_var_sup\tpost_mean_sub\tpost_var_sub\tpost_mean_sup\tpost_var_sup\t\tpval",
                        feats, NULL, 9, prior_means_sub, prior_vars_sub, 
                        prior_means_sup, prior_vars_sup, post_means_sub, 
                        post_vars_sub, post_means_sup, post_vars_sup, pvals);
  if (!output_gff) {
    free(post_means_sub);
    free(post_vars_sub);
    free(post_means_sup);
    free(post_vars_sup);
    free(prior_means_sub);
    free(prior_vars_sub);
    free(prior_means_sup);
    free(prior_vars_sup);
  }
  free(pvals);
}

void print_wig(MSA *msa, double *vals, char *chrom, int refidx, 
               int log_trans) {
  int last, j, k;
  double val;
  last = -INFTY;
  assert(refidx >= 0 && refidx <= msa->nseqs);
  for (j = 0, k = 0; j < msa->length; j++) {
    if (refidx == 0 || msa_get_char(msa, refidx-1, j) != GAP_CHAR) {
      if (refidx == 0 || !msa_missing_col(msa, refidx, j)) {
        if (k > last + 1) 
          printf("fixedStep chrom=%s start=%d step=1\n", chrom, 
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
        printf("%.3f\n", val);
        last = k;
      }
      k++;
    }
  }
}

/* Print arbitrary columns of tuple-specific data in wig-like format */
void print_base_by_base(char *header, char *chrom, MSA *msa, 
                        char **formatstr, int refidx, int ncols, ...) {
  int last, j, k, tup, col;
  va_list ap;
  double *data[ncols];
  assert(refidx >= 0 && refidx <= msa->nseqs);

  last = -INFTY;
  if (header != NULL)
    printf("%s\n", header);

  va_start(ap, ncols);
  for (col = 0; col < ncols; col++)
    data[col] = va_arg(ap, double*);

  for (j = 0, k = 0; j < msa->length; j++) {
    if (refidx == 0 || msa_get_char(msa, refidx-1, j) != GAP_CHAR) {
      if (refidx == 0 || !msa_missing_col(msa, refidx, j)) {
        if (k > last + 1) 
          printf("fixedStep chrom=%s start=%d step=1\n", chrom, 
                 k + msa->idx_offset + 1);
        tup = msa->ss->tuple_idx[j];
        for (col = 0; col < ncols; col++) {
          printf((formatstr == NULL ? "%.5f" : formatstr[col]), data[col][tup]);
          printf(col < ncols-1 ? "\t" : "\n");
        }
        last = k;
      }
      k++;
    }
  }
  va_end(ap);
}

/* Print a list of features and artibrary associated statistics */
void print_feats_generic(char *header, GFF_Set *gff, char **formatstr, 
                         int ncols, ...) {
  int i, col;
  String *name;
  va_list ap;
  double *data[ncols];
  Regex *tag_val_re = str_re_new("[[:alnum:]_.]+[[:space:]]+(\"[^\"]*\"|[^[:space:]]+)");
  List *l = lst_new_ptr(2);

  printf("#chr\tstart\tend\tname");
  if (header != NULL) 
    printf("\t%s\n", header);
  else 
    printf("\n");

  va_start(ap, ncols);
  for (col = 0; col < ncols; col++)
    data[col] = va_arg(ap, double*);

  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *f = lst_get_ptr(gff->features, i);

    /* try to extract feature name from attribute field */
    lst_clear(l);
    if (f->attribute->length > 0 && 
        str_re_match(f->attribute, tag_val_re, l, 1) >= 0) {
      name = lst_get_ptr(l, 1);
      str_remove_quotes(name);
    } else name=NULL;

    printf("%s\t%d\t%d\t%s\t", f->seqname->chars, f->start-1, f->end, 
           name == NULL ? "." : name->chars);

    for (col = 0; col < ncols; col++) {
      printf((formatstr == NULL ? "%.5f" : formatstr[col]), data[col][i]);
      printf(col < ncols-1 ? "\t" : "\n");
    }

    lst_free_strings(l);
  }
  va_end(ap);
  lst_free(l);
  str_re_free(tag_val_re);
}

/* Print GFF to stdout with feature scores defined by vals.  If
   log_trans == TRUE, take log transform (propagating negative
   signs) */
void print_gff_scores(GFF_Set *gff, double *vals, int log_trans) {
  int i;
  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *f = lst_get_ptr(gff->features, i);
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
  gff_print_set(stdout, gff);
}
