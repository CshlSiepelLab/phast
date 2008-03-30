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

/* maximum size of matrix for which to do explicit convolution of
   joint prior; beyond this size an approximation is used.
   Computational complexity is proportional to square of this number.
   This only comes into play when --features and --subtree are used
   together */
#define MAX_CONVOLVE_SIZE 22500


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

void print_prior_metadata(char *mod_fname, Vector *prior_distrib) {
  int i;
  double prior_mean, prior_var;
  Vector *cdf = pv_cdf(prior_distrib, UPPER);
  pv_stats(prior_distrib, &prior_mean, &prior_var);
  printf("#metadata for null model '%s'\n", mod_fname);
  printf("mean = %.3f\ns.d. = %.3f\nmax = %d\n-logP =", 
         prior_mean, prior_var, prior_distrib->size);
  for (i = 0; i < cdf->size; i++)
    printf("% .3f", fabs(-log10(cdf->data[i])));  
  printf("\n");
  vec_free(cdf);
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
  printf("p-value of anti-conservation: %e\n\n",
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
  printf("p-value of anti-conservation in subtree: %e\n\n", anti_cons_p_sub);

  printf("p-value of conservation in supertree: %e\n", cons_p_sup);
  printf("p-value of anti-conservation in supertree: %e\n\n", anti_cons_p_sup);

  printf("p-value of conservation in subtree given total: %e\n", cond_cons_p_sub);
  printf("p-value of anti-conservation in subtree given total: %e\n\n", cond_anti_cons_p_sub);

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

void print_p_feats(JumpProcess *jp, MSA *msa, GFF_Set *feats, double ci) {
  int i;
  Regex *tag_val_re = str_re_new("[[:alnum:]_.]+[[:space:]]+(\"[^\"]*\"|[^[:space:]]+)");
  p_value_stats *stats = sub_p_value_many(jp, msa, feats->features, ci);
  List *l = lst_new_ptr(2);

  msa_map_gff_coords(msa, feats, 0, 1, 0, NULL);

  printf("#chr\tstart\tend\tname\tp_cons\tp_anti_cons\tprior_mean\tprior_var\tprior_min\tprior_max\tpost_mean\tpost_var\tpost_min\tpost_max\n");
  for (i = 0; i < lst_size(feats->features); i++) {
    GFF_Feature *f = lst_get_ptr(feats->features, i);
    String *name = NULL;

    /* try to extract feature name from attribute field */
    lst_clear(l);
    if (f->attribute->length > 0 && 
        str_re_match(f->attribute, tag_val_re, l, 1) >= 0) {
      name = lst_get_ptr(l, 1);
      str_remove_quotes(name);
    }

    printf("%s\t%d\t%d\t%s\t%e\t%e\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t%d\t%d\n", 
           f->seqname->chars, f->start-1, f->end, 
           name == NULL ? "." : name->chars,
           stats[i].p_cons, stats[i].p_anti_cons, 
           stats[i].prior_mean, stats[i].prior_var, 
           stats[i].prior_min, stats[i].prior_max, 
           stats[i].post_mean, stats[i].post_var, 
           stats[i].post_min, stats[i].post_max);

    lst_free_strings(l);
  }
  lst_free(l);
  str_re_free(tag_val_re);
}

void print_p_joint_feats(JumpProcess *jp, MSA *msa, GFF_Set *feats, double ci) {
  int i;
  Regex *tag_val_re = str_re_new("[[:alnum:]_.]+[[:space:]]+(\"[^\"]*\"|[^[:space:]]+)");
  p_value_joint_stats *stats = 
    sub_p_value_joint_many(jp, msa, feats->features, 
                           ci, MAX_CONVOLVE_SIZE, NULL);
  List *l = lst_new_ptr(2);

  msa_map_gff_coords(msa, feats, 0, 1, 0, NULL);

  printf("#chr\tstart\tend\tname\tp_cons_sup\tp_anti_cons_sup\tp_cons_sub\tp_anti_cons_sub\tcond_p_cons_sub\tcond_p_anti_cons_sub\tcond_approx\tprior_mean_sup\tprior_var_sup\tprior_min_sup\tprior_max_sup\tprior_mean_sub\tprior_var_sub\tprior_min_sub\tprior_max_sub\tpost_mean_sup\tpost_var_sup\tpost_min_sup\tpost_max_sup\tpost_mean_sub\tpost_var_sub\tpost_min_sub\tpost_max_sub\n");
  for (i = 0; i < lst_size(feats->features); i++) {
    GFF_Feature *f = lst_get_ptr(feats->features, i);
    String *name = NULL;

    /* try to extract feature name from attribute field */
    lst_clear(l);
    if (f->attribute->length > 0 && 
        str_re_match(f->attribute, tag_val_re, l, 1) >= 0) {
      name = lst_get_ptr(l, 1);
      str_remove_quotes(name);
    }

    printf("%s\t%d\t%d\t%s\t%e\t%e\t%e\t%e\t%e\t%e\t%s\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t%d\t%d\n", 
           f->seqname->chars, f->start-1, f->end, 
           name == NULL ? "." : name->chars,
           stats[i].p_cons_right, stats[i].p_anti_cons_right, 
           stats[i].p_cons_left, stats[i].p_anti_cons_left, 
           stats[i].cond_p_cons_left, stats[i].cond_p_anti_cons_left, 
           stats[i].cond_p_approx ? "approx" : "exact",
           stats[i].prior_mean_right, stats[i].prior_var_right, 
           stats[i].prior_min_right, stats[i].prior_max_right, 
           stats[i].prior_mean_left, stats[i].prior_var_left, 
           stats[i].prior_min_left, stats[i].prior_max_left, 
           stats[i].post_mean_right, stats[i].post_var_right, 
           stats[i].post_min_right, stats[i].post_max_right,
           stats[i].post_mean_left, stats[i].post_var_left, 
           stats[i].post_min_left, stats[i].post_max_left);

    lst_free_strings(l);
  }
  lst_free(l);
  str_re_free(tag_val_re);
}

void print_wig(MSA *msa, double *vals, char *chrom, int log_trans) {
  int last, j, k;
  double val;
  last = -INFTY;
  for (j = 0, k = 0; j < msa->length; j++) {
    if (msa_get_char(msa, 0, j) != GAP_CHAR) {
      if (!msa_missing_col(msa, 1, j)) {
        if (k > last + 1) 
          printf("fixedStep chrom=%s start=%d step=1\n", chrom, 
                 k + msa->idx_offset + 1);
        val = vals[msa->ss->tuple_idx[j]];
        if (log_trans) val = fabs(-log10(val));
        printf("%.3f\n", val);
        last = k;
      }
      k++;
    }
  }
}

/* print arbitrary columns of tuple-specific data in wig-like format */
void print_base_by_base(char *header, char *chrom, MSA *msa, int ncols, 
                        char **formatstr, ...) {
  int last, j, k, tup, col;
  va_list ap;
  double *data[ncols];

  last = -INFTY;
  if (header != NULL)
    printf("%s\n", header);

  va_start(ap, ncols);
  for (col = 0; col < ncols; col++)
    data[col] = va_arg(ap, double*);

  for (j = 0, k = 0; j < msa->length; j++) {
    if (msa_get_char(msa, 0, j) != GAP_CHAR) {
      if (!msa_missing_col(msa, 1, j)) {
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
