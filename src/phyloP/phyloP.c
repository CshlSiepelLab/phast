#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <msa.h>
#include <maf.h>
#include <tree_model.h>
#include <sufficient_stats.h>
#include <subst_distrib.h>
#include <prob_vector.h>
#include <prob_matrix.h>
#include "phyloP.help"

/* maximum size of matrix for which to do explicit convolution of
   joint prior; beyond this size an approximation is used.
   Computational complexity is proportional to square of this number.
   This only comes into play when --features and --subtree are used
   together */
#define MAX_CONVOLVE_SIZE 22500

void print_prior_only(int nsites, char *mod_fname, Vector *prior_distrib);
void print_post_only(char *mod_fname, char *msa_fname, Vector *post_distrib,
                     double ci);
void print_p(char *mod_fname, char *msa_fname, Vector *prior_distrib,
             double post_mean, double post_var, double ci);
void print_prior_only_joint(char *node_name, int nsites, char *mod_fname, 
                            Matrix *prior_distrib);
void print_post_only_joint(char *node_name, char *mod_fname, 
                           char *msa_fname, Matrix *post_distrib, 
                           double ci);
void print_p_joint(char *node_name, char *mod_fname, char *msa_fname, 
                   double ci, Matrix *prior_joint, 
                   double post_mean, double post_var, 
                   double post_mean_sup, double post_var_sup, 
                   double post_mean_sub, double post_var_sub);
void print_p_feats(JumpProcess *jp, MSA *msa, GFF_Set *feats, double ci);
void print_p_joint_feats(JumpProcess *jp, MSA *msa, GFF_Set *feats, double ci);
void print_quantiles(Vector *distrib);

int main(int argc, char *argv[]) {
  /* variables for options with defaults */
  msa_format_type msa_format = FASTA;
  int nsites = -1, prior_only = FALSE, post_only = FALSE, quantiles = FALSE;
  double ci = -1;
  char *subtree_name = NULL;
  GFF_Set *feats = NULL;

  /* other variables */
  FILE *msa_f = NULL;
  TreeModel *mod;
  MSA *msa;
  TreeNode *subtree_root = NULL;
  Vector *prior_distrib, *post_distrib;
  Matrix *prior_joint_distrib, *post_joint_distrib;
  JumpProcess *jp;
  char c;
  int opt_idx;

  struct option long_opts[] = {
    {"msa-format", 1, 0, 'i'},
    {"null", 1, 0, 'n'},
    {"posterior", 0, 0, 'p'},
    {"confidence-interval", 1, 0, 'c'},
    {"subtree", 1, 0, 's'},
    {"features", 1, 0, 'f'},
    {"quantiles", 0, 0, 'q'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "i:n:pc:s:f:qh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == -1)
        die("ERROR: unrecognized alignment format.\n");
      break;
    case 'n':
      nsites = get_arg_int_bounds(optarg, 1, INFTY);
      prior_only = TRUE;
      break;
    case 'p':
      post_only = TRUE;
      break;
    case 'c':
      ci = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 's':
      subtree_name = optarg;
      break;
    case 'f':
      feats = gff_read_set(fopen_fname(optarg, "r"));
      break;
    case 'q':
      quantiles = TRUE;
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'phyloP -h'.\n");
    }
  }

  if ((prior_only && optind > argc - 1) || 
      (!prior_only && optind != argc - 2))
    die("ERROR: bad arguments.  Try 'phyloP -h'.\n");

  if (quantiles && !prior_only && !post_only)
    die("ERROR: --quantiles can only be used with --null or --posterior.\n");
  if (quantiles && subtree_name != NULL)
    die("ERROR: --quantiles cannot be used with --subtree.\n");
  if (feats != NULL && (prior_only || post_only))
    die("ERROR: --features cannot be used with --null or --posterior.\n");

  mod = tm_new_from_file(fopen_fname(argv[optind], "r"));
  jp = sub_define_jump_process(mod);

  if (!prior_only) {
    msa_f = fopen_fname(argv[optind+1], "r");
    if (msa_format == MAF) 
      msa = maf_read(msa_f, NULL, 1, NULL, NULL, NULL, -1, 
                     feats == NULL ? FALSE : TRUE, /* --features requires order */
                     NULL, NO_STRIP, FALSE); 
    else 
      msa = msa_new_from_file(msa_f, msa_format, NULL);

    if (msa->ss == NULL)
      ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1);

    if (msa_alph_has_lowercase(msa)) msa_toupper(msa);     
    msa_remove_N_from_alph(msa);

    if (feats != NULL && msa->ss->tuple_idx == NULL)
      die("ERROR: ordered alignment required.\n");
  }

  if (nsites == -1) nsites = msa->length;

  if (feats != NULL) 
    msa_map_gff_coords(msa, feats, 1, 0, 0, NULL);

  if (subtree_name == NULL) {
    double post_mean, post_var;

    if (feats == NULL) {
      /* compute distributions and stats*/
      if (!post_only) 
        prior_distrib = sub_prior_distrib_alignment(jp, nsites);

      if (post_only)	 /* don't need explicit distrib for p-value */
        post_distrib = sub_posterior_distrib_alignment(jp, msa);
      else if (!prior_only)
        sub_posterior_stats_alignment(jp, msa, &post_mean, &post_var);

      /* print output */
      if (quantiles)
        print_quantiles(prior_only ? prior_distrib : post_distrib);
      else if (prior_only)
        print_prior_only(nsites, argv[optind], prior_distrib);
      else if (post_only)
        print_post_only(argv[optind], argv[optind+1], post_distrib, ci);
      else
        print_p(argv[optind], argv[optind+1], prior_distrib, 
                post_mean, post_var, ci);
    }
    else                        /* --features case */
      print_p_feats(jp, msa, feats, ci);
  }
  else {			/* supertree/subtree mode */
    double post_mean, post_var, post_mean_sup, post_var_sup, 
      post_mean_sub, post_var_sub;
    TreeNode *tmp;

    if (!tm_is_reversible(mod->subst_mod))
      die("ERROR: reversible model required with --subtree.\n");

    /* reroot tree */
    tr_name_ancestors(mod->tree);
    subtree_root = tr_get_node(mod->tree, subtree_name);
    if (subtree_root == NULL) 
      die("ERROR: no node named '%s'.\n", subtree_name);

    tr_reroot(mod->tree, subtree_root, TRUE);
    mod->tree = subtree_root->parent; /* take parent because including
                                         branch */

    /* swap left and right children.  This is necessary because
       routines for computing joint distrib assume branch to right has
       length zero, but because branch is included, tr_reroot will put
       zero length branch on left */
    tmp = mod->tree->lchild;
    mod->tree->lchild = mod->tree->rchild;
    mod->tree->rchild = tmp;

    if (feats == NULL) {
      /* compute distributions and stats */
      if (!post_only)
        prior_joint_distrib = sub_prior_joint_distrib_alignment(jp, nsites);

      if (post_only)
        post_joint_distrib = sub_posterior_joint_distrib_alignment(jp, msa);
      else if (!prior_only)
        sub_posterior_joint_stats_alignment(jp, msa, &post_mean, &post_var,
                                            &post_mean_sub, &post_var_sub, 
                                            &post_mean_sup, &post_var_sup);

      /* print output */
      if (prior_only) 
        print_prior_only_joint(subtree_root->name, nsites, argv[optind],
                               prior_joint_distrib);
      else if (post_only) 
        print_post_only_joint(subtree_root->name, argv[optind], 
                              argv[optind+1], post_joint_distrib, ci);
      else 
        print_p_joint(subtree_root->name, argv[optind], argv[optind+1],
                      ci, prior_joint_distrib, post_mean, post_var, 
                      post_mean_sup, post_var_sup, post_mean_sub, 
                      post_var_sub);
    }
    else                        /* --features case */
      print_p_joint_feats(jp, msa, feats, ci);
  }
    
  return 0;
}

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
                     double ci) {
  int i, post_min, post_max;
  double post_mean, post_var;
  if (ci == -1) ci = 0.95;      /* for purposes of stats */
  pv_stats(post_distrib, &post_mean, &post_var);
  pv_confidence_interval(post_distrib, ci, &post_min, &post_max);
  printf("#Let n be no. substitutions given '%s' and '%s'.\n", 
         mod_fname, msa_fname);
  printf("#E[n] = %.3f; Var[n] = %.3f; %.1f%% c.i. = [%d, %d]\n", 
         post_mean, post_var, ci*100, post_min, post_max);
  printf("#n p(n)\n");
  for (i = 0; i < post_distrib->size; i++)
    printf("%d\t%f\n", i, post_distrib->data[i]);
}

void print_p(char *mod_fname, char *msa_fname, Vector *prior_distrib,
             double post_mean, double post_var, double ci) {
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
                           double ci) {
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
                   double post_mean_sub, double post_var_sub) {

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
}

void print_p_feats(JumpProcess *jp, MSA *msa, GFF_Set *feats, double ci) {
  int i;
  Regex *tag_val_re = str_re_new("[[:alnum:]_.]+[[:space:]]+(\"[^\"]*\"|[^[:space:]]+)");
  p_value_stats *stats = sub_p_value_many(jp, msa, feats->features, ci);
  List *l = lst_new_ptr(2);

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
