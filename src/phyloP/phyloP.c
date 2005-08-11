#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <msa.h>
#include <maf.h>
#include <tree_model.h>
#include <subst_distrib.h>
#include <prob_vector.h>
#include <prob_matrix.h>
#include "phyloP.help"

void print_quantiles(Vector *distrib);
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
  jp = sub_define_jump_process(mod, max(20, 10 * tr_total_len(mod->tree)));

  if (!prior_only) {
    msa_f = fopen_fname(argv[optind+1], "r");
    if (msa_format == MAF) 
      msa = maf_read(msa_f, NULL, 1, NULL, NULL, -1, FALSE, NULL, 
                     NO_STRIP, FALSE); 
    else 
      msa = msa_new_from_file(msa_f, msa_format, NULL);
  }

  if (nsites == -1) nsites = msa->length;

  if (subtree_name == NULL) {
    double post_mean, post_var;

    /* compute distributions and stats*/
    if (!post_only) 
      prior_distrib = sub_prior_distrib_alignment(jp, nsites);

    if (post_only)	 /* don't need explicit distrib for p-value */
      post_distrib = sub_posterior_distrib_alignment(jp, msa);
    else if (!prior_only)
      sub_posterior_stats_alignment(jp, msa, &post_mean, &post_var);

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
  else {			/* supertree/subtree mode */
    double post_mean, post_var, post_mean_sup, post_var_sup, 
      post_mean_sub, post_var_sub;

    if (!tm_is_reversible(mod->subst_mod))
      die("ERROR: reversible model required with --subtree.\n");

    subtree_root = tr_get_node(mod->tree, subtree_name);
    if (subtree_root == NULL) 
      die("ERROR: no node named '%s'.\n", subtree_name);

    tr_reroot(mod->tree, subtree_root);
    mod->tree = subtree_root;

    if (!post_only)
      prior_joint_distrib = sub_prior_joint_distrib_alignment(jp, nsites);

    if (post_only)
      post_joint_distrib = sub_posterior_joint_distrib_alignment(jp, msa);
    else if (!prior_only)
      sub_posterior_joint_stats_alignment(jp, msa, &post_mean, &post_var,
                                          &post_mean_sup, &post_var_sup, 
                                          &post_mean_sub, &post_var_sub);

    if (prior_only) 
      print_prior_only_joint(subtree_root->name, nsites, argv[optind],
                             prior_joint_distrib);
    else if (post_only) 
      print_post_only_joint(subtree_root->name, argv[optind], 
                            argv[optind+1], post_joint_distrib, ci);
    else 
      print_p_joint(subtree_root->name, argv[optind], argv[optind+1]
                    , ci, prior_joint_distrib, post_mean, post_var, 
                    post_mean_sup, post_var_sup, post_mean_sub, post_var_sub);
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
  pv_confidence_interval(prior_distrib, ci, &prior_min, &prior_max);

  /* for posterior, avoid computing convolution and assume normality
     based on central limit theorem */
  if (ci != -1)
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
  printf("null distrib: mean = %f, var = %f, %.1f%% c.i. = [%d, %d]\n", 
         prior_mean, prior_var, ci*100, prior_min, prior_max);
  printf("posterior distrib: mean = %f, var = %f, %.1f%% c.i. = [%.0f, %.0f]\n\n",
         post_mean, post_var, ci*100, post_min, post_max);
}

void print_prior_only_joint(char *node_name, int nsites, char *mod_fname, 
                            Matrix *prior_distrib) {
  int i, j, min_sup, max_sup, min_sub, max_sub;
  double mean_sup, var_sup, mean_sub, var_sub;
  Vector *marg_sup = vec_new(prior_distrib->nrows),
    *marg_sub = vec_new(prior_distrib->ncols);

  /* get marginal distributions for stats */
  vec_zero(marg_sup); vec_zero(marg_sub);
  for (i = 0; i < prior_distrib->nrows; i++) {
    for (j = 0; j < prior_distrib->ncols; j++) {
      marg_sup->data[i] += prior_distrib->data[i][j];
      marg_sub->data[j] += prior_distrib->data[i][j];
    }
  }

  /* get marginal stats */
  pv_stats(marg_sup, &mean_sup, &var_sup);
  pv_stats(marg_sub, &mean_sub, &var_sub);
  pv_confidence_interval(marg_sup, 0.95, &min_sup, &max_sup);
  pv_confidence_interval(marg_sub, 0.95, &min_sub, &max_sub);

  printf("#Let n1 be no. substitutions in supertree above '%s' over %d site(s) given '%s'.\n", 
         node_name, nsites, mod_fname);
  printf("#Let n2 be no. substitutions in subtree beneath '%s' over %d site(s) given '%s'.\n", 
         node_name, nsites, mod_fname);
  printf("#E[n1] = %.3f; Var[n1] = %.3f; 95%% c.i. = [%d, %d]\n", 
         mean_sup, var_sup, min_sup, max_sup);
  printf("#E[n2] = %.3f; Var[n2] = %.3f; 95%% c.i. = [%d, %d]\n", 
         mean_sub, var_sub, min_sub, max_sub);
  printf("\n#element at row n1 and col n2 in table below is p(n1, n2)\n");
  for (i = 0; i < prior_distrib->nrows; i++) 
    for (j = 0; j < prior_distrib->ncols; j++) 
      printf("%f%c", prior_distrib->data[i][j], 
             j == prior_distrib->ncols - 1 ? '\n' : '\t');

  vec_free(marg_sup);
  vec_free(marg_sub);
}

void print_post_only_joint(char *node_name, char *mod_fname, 
                           char *msa_fname, Matrix *post_distrib, 
                           double ci) {
  int i, j, min_sup, max_sup, min_sub, max_sub;
  double mean_sup, var_sup, mean_sub, var_sub;
  Vector *marg_sup = vec_new(post_distrib->nrows),
    *marg_sub = vec_new(post_distrib->ncols);

  if (ci == -1) ci= 0.95;       /* for purposes of stats */

  /* get marginal distributions for stats */
  vec_zero(marg_sup); vec_zero(marg_sub);
  for (i = 0; i < post_distrib->nrows; i++) {
    for (j = 0; j < post_distrib->ncols; j++) {
      marg_sup->data[i] += post_distrib->data[i][j];
      marg_sub->data[j] += post_distrib->data[i][j];
    }
  }

  /* get marginal stats */
  pv_stats(marg_sup, &mean_sup, &var_sup);
  pv_stats(marg_sub, &mean_sub, &var_sub);
  pv_confidence_interval(marg_sup, ci, &min_sup, &max_sup);
  pv_confidence_interval(marg_sub, ci, &min_sub, &max_sub);

  printf("#Let n1 be no. substitutions in supertree above '%s' given '%s' and '%s'.\n", 
         node_name, mod_fname, msa_fname);
  printf("#Let n2 be no. substitutions in subtree beneath '%s' given '%s' and '%s'.\n", 
         node_name, mod_fname, msa_fname);
  printf("#E[n1] = %.3f; Var[n1] = %.3f; %.1f%% c.i. = [%d, %d]\n", 
         mean_sup, var_sup, ci*100, min_sup, max_sup);
  printf("#E[n2] = %.3f; Var[n2] = %.3f; %.1f%% c.i. = [%d, %d]\n", 
         mean_sub, var_sub, ci*100, min_sub, max_sub);
  printf("\n#element at row n1 and col n2 in table below is p(n1, n2)\n");
  for (i = 0; i < post_distrib->nrows; i++) 
    for (j = 0; j < post_distrib->ncols; j++) 
      printf("%f%c", post_distrib->data[i][j], 
             j == post_distrib->ncols - 1 ? '\n' : '\t');

  vec_free(marg_sup);
  vec_free(marg_sub);
}

void print_p_joint(char *node_name, char *mod_fname, char *msa_fname, 
                   double ci, Matrix *prior_joint, 
                   double post_mean, double post_var, 
                   double post_mean_sup, double post_var_sup, 
                   double post_mean_sub, double post_var_sub) {

  double post_min_tot, post_max_tot, post_min_sup, post_max_sup, 
    post_min_sub, post_max_sub, cons_p_sup, anti_cons_p_sup, 
    cons_p_sub, anti_cons_p_sub, prior_mean_sup, prior_var_sup, 
    prior_mean_sub, prior_var_sub;
  int prior_min_sup, prior_max_sup, prior_min_sub, prior_max_sub;
  Vector *prior_marg_sup, *prior_marg_sub, *cond;

  /* To be conservative, base the conservation p value on the largest
     reasonable estimate of the number of subst. in the subtree and
     the smallest reasonable estimate of the number of subst. in the
     whole tree, and base the anti-conservation p-value on the
     smallest reasonable estimate of the number of subst. in the
     supertree and the largest reasonable estimate of the number of
     subst. in the whole tree */

  if (ci == -1) {
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

  cond = pm_x_given_tot(prior_joint, post_min_tot);
  cons_p_sup = pv_p_value(cond, post_max_sup, LOWER);
  vec_free(cond);

  cond = pm_x_given_tot(prior_joint, post_max_tot);
  anti_cons_p_sup = pv_p_value(cond, post_min_sup, UPPER);
  vec_free(cond);

  cond = pm_y_given_tot(prior_joint, post_min_tot);
  cons_p_sub = pv_p_value(cond, post_max_sub, LOWER);
  vec_free(cond);

  cond = pm_y_given_tot(prior_joint, post_max_tot);
  anti_cons_p_sub = pv_p_value(cond, post_min_sub, UPPER);
  vec_free(cond);

  prior_marg_sup = pm_marg_x(prior_joint);
  prior_marg_sub = pm_marg_y(prior_joint);
  pv_stats(prior_marg_sup, &prior_mean_sup, &prior_var_sup);
  pv_stats(prior_marg_sub, &prior_mean_sub, &prior_var_sub);
  pv_confidence_interval(prior_marg_sup, ci, &prior_min_sup, &prior_max_sup);
  pv_confidence_interval(prior_marg_sub, ci, &prior_min_sub, &prior_max_sub);
  vec_free(prior_marg_sup);
  vec_free(prior_marg_sub);


  printf("\n*****\nP-values for number of substitutions observed in '%s' given '%s',\n", 
         msa_fname, mod_fname);
  printf ("considering subtree beneath node '%s', supertree above node '%s'\n*****\n\n", node_name, node_name);
  printf("p-value of conservation in subtree: %e\n", cons_p_sub);
  printf("p-value of anti-conservation in subtree: %e\n\n", anti_cons_p_sub);
  printf("p-value of conservation in supertree: %e\n", cons_p_sup);
  printf("p-value of anti-conservation in supertree: %e\n\n", anti_cons_p_sup);

  printf("null distrib in subtree: mean = %f, var = %f, %.1f%% c.i. = [%d, %d]\n", 
         prior_mean_sub, prior_var_sub, ci*100, prior_min_sub, prior_max_sub);
  printf("posterior distrib in subtree: mean = %f, var = %f, %.1f%% c.i. = [%.0f, %.0f]\n\n",
         post_mean_sub, post_var_sub, ci*100, post_min_sub, post_max_sub);
  printf("null distrib in supertree: mean = %f, var = %f, %.1f%% c.i. = [%d, %d]\n", 
         prior_mean_sup, prior_var_sup, ci*100, prior_min_sup, prior_max_sup);
  printf("posterior distrib in supertree: mean = %f, var = %f, %.1f%% c.i. = [%.0f, %.0f]\n\n",
         post_mean_sup, post_var_sup, ci*100, post_min_sup, post_max_sup);
}
