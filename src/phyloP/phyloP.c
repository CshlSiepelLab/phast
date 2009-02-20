/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

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
#include "phyloP.h"
#include "fit_column.h"
#include "fit_feature.h"
#include "phyloP.help"

/* default values for epsilon; can tolerate larger value with --wig-scores or
   --base-by-base */
#define DEFAULT_EPSILON 1e-10
#define DEFAULT_EPSILON_BASE_BY_BASE 1e-6

typedef enum{SPH, LRT, SCORE, GERP} method_type;

/* maximum size of matrix for which to do explicit convolution of
   joint prior; beyond this size an approximation is used.
   Computational complexity is proportional to square of this number.
   This only comes into play when --features and --subtree are used
   together */
#define MAX_CONVOLVE_SIZE 22500

TreeModel* fit_tree_model(TreeModel *source_mod, MSA *msa, 
                          char *subtree_name, double *scale,
                          double *sub_scale);

int main(int argc, char *argv[]) {
  /* variables for options with defaults */
  msa_format_type msa_format = FASTA;
  int nsites = -1, prior_only = FALSE, post_only = FALSE, quantiles = FALSE,
    fit_model = FALSE, base_by_base = FALSE, output_wig = FALSE, 
    default_epsilon = TRUE, output_gff = FALSE, refidx = 1;
  double ci = -1, epsilon = DEFAULT_EPSILON;
  char *subtree_name = NULL, *chrom = NULL;
  GFF_Set *feats = NULL;
  method_type method = SPH;
  mode_type mode = CON;
  FILE *logf = NULL;

  /* other variables */
  FILE *msa_f = NULL;
  TreeModel *mod, *mod_fitted = NULL;
  MSA *msa;
  Vector *prior_distrib, *post_distrib;
  Matrix *prior_joint_distrib, *post_joint_distrib;
  JumpProcess *jp, *jp_post;
  List *pruned_names;
  char c;
  int j, old_nleaves, opt_idx;
  double scale = -1, sub_scale = -1;
  double prior_mean, prior_var;
  double *pvals = NULL, *post_means = NULL, *post_vars = NULL,
    *llrs = NULL, *scales = NULL, *sub_scales = NULL, 
    *null_scales = NULL, *derivs = NULL, *sub_derivs = NULL,
    *teststats = NULL, *nneut = NULL, *nobs = NULL,
    *nspec = NULL, *nrejected = NULL;

  struct option long_opts[] = {
    {"method", 1, 0, 'm'},
    {"mode", 1, 0, 'o'},
    {"msa-format", 1, 0, 'i'},
    {"null", 1, 0, 'n'},
    {"posterior", 0, 0, 'p'},
    {"confidence-interval", 1, 0, 'c'},
    {"subtree", 1, 0, 's'},
    {"features", 1, 0, 'f'},
    {"fit-model", 0, 0, 'F'},
    {"epsilon", 1, 0, 'e'},
    {"quantiles", 0, 0, 'q'},
    {"wig-scores", 0, 0, 'w'},
    {"base-by-base", 0, 0, 'b'},
    {"refidx", 1, 0, 'r'},
    {"chrom", 1, 0, 'N'},
    {"log", 1, 0, 'l'},
    {"gff-scores", 0, 0, 'g'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "m:o:i:n:pc:s:f:Fe:l:r:qwgbN:h", 
                          long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'm':
      if (!strcmp(optarg, "SPH"))
        method = SPH;
      else if (!strcmp(optarg, "LRT"))
        method = LRT;
      else if (!strcmp(optarg, "SCORE"))
        method = SCORE;
      else if (!strcmp(optarg, "GERP"))
        method = GERP;
      else die("ERROR: bad argument to --method (-m).\n");
      break;
    case 'o':
      if (!strcmp(optarg, "CON"))
        mode = CON;
      else if (!strcmp(optarg, "ACC"))
        mode = ACC;
      else if (!strcmp(optarg, "NNEUT"))
        mode = NNEUT;
      else if (!strcmp(optarg, "CONACC"))
        mode = CONACC;
      else die("ERROR: bad argument to --mode (-o).\n");
      break;
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == -1)
        die("ERROR: unrecognized alignment format.\n");
      break;
    case 'g':
      output_gff = TRUE;
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
    case 'F':
      fit_model = TRUE;
      break;
    case 'e':
      epsilon = get_arg_dbl_bounds(optarg, 0, 0.1);
      default_epsilon = FALSE;
      break;
    case 'q':
      quantiles = TRUE;
      break;
    case 'w':
      base_by_base = TRUE;
      output_wig = TRUE;
      break;
    case 'b':
      base_by_base = TRUE;
      break;
    case 'N':
      chrom = optarg;
      break;
    case 'r':
      refidx = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'l':
      if (!strcmp(optarg, "-"))
        logf = stderr;
      else
        logf = fopen_fname(optarg, "w+");
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'phyloP -h'.\n");
    }
  }

  if (method != SPH && (prior_only || post_only || fit_model || 
                        !default_epsilon || quantiles || ci != -1))
    die("ERROR: bad arguments.  Try 'phyloP -h'.\n");

  if ((prior_only && optind > argc - 1) || 
      (!prior_only && optind != argc - 2))
    die("ERROR: bad arguments.  Try 'phyloP -h'.\n");

  if (quantiles && !prior_only && !post_only)
    die("ERROR: --quantiles can only be used with --null or --posterior.\n");
  if (quantiles && subtree_name != NULL)
    die("ERROR: --quantiles cannot be used with --subtree.\n");
  if (feats != NULL && (prior_only || post_only || fit_model))
    die("ERROR: --features cannot be used with --null, --posterior, or --fit-model.\n");
  if (base_by_base && (prior_only || post_only || ci != -1 || feats != NULL))
    die("ERROR: --wig-scores and --base-by-base cannot be used with --null, --posterior, --features, --quantiles, or --confidence-interval.\n");
  if (method == GERP && subtree_name != NULL)
    die("ERROR: --subtree not supported with --method GERP.\n");

  mod = tm_new_from_file(fopen_fname(argv[optind], "r"));

  if (!prior_only) {
    msa_f = fopen_fname(argv[optind+1], "r");
    if (msa_format == MAF) 
      msa = maf_read(msa_f, NULL, 1, NULL, NULL, NULL, -1, 
                     (feats == NULL && base_by_base==0) ? FALSE : TRUE, /* --features requires order */
                     NULL, NO_STRIP, FALSE); 
    else 
      msa = msa_new_from_file(msa_f, msa_format, NULL);

    /* if base_by_base and undefined chrom, use filename root as chrom */
    if (base_by_base && chrom == NULL) {
      String *tmpstr = str_new_charstr(argv[optind+1]);
      if (str_equals_charstr(tmpstr, "-")) chrom = "NA";
      else {
        str_remove_path(tmpstr);
        str_shortest_root(tmpstr, '.');
        chrom = tmpstr->chars;    
      }
    }

    if (msa->ss == NULL)
      ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1);

    if (msa_alph_has_lowercase(msa)) msa_toupper(msa);     
    msa_remove_N_from_alph(msa);

    if ((feats != NULL || base_by_base) && msa->ss->tuple_idx == NULL)
      die("ERROR: ordered alignment required.\n");

    /* prune tree, if necessary */
    pruned_names = lst_new_ptr(msa->nseqs);
    old_nleaves = (mod->tree->nnodes + 1) / 2;
    tm_prune(mod, msa, pruned_names);
    if (lst_size(pruned_names) >= old_nleaves)
      die("ERROR: no match for leaves of tree in alignment.\n");
    else if (lst_size(pruned_names) > 0) {
      fprintf(stderr, "WARNING: pruned away leaves with no match in alignment (");
      for (j = 0; j < lst_size(pruned_names); j++)
        fprintf(stderr, "%s%s", ((String*)lst_get_ptr(pruned_names, j))->chars, 
                j < lst_size(pruned_names) - 1 ? ", " : ").\n");
    }
  }

  /* set subtree if necessary */
  if (subtree_name != NULL && method != SPH) {
    /* (SPH is a special case -- requires rerooting) */
    mod->subtree_root = tr_get_node(mod->tree, subtree_name);
    if (mod->subtree_root == NULL)
      die("ERROR: no node named '%s'.\n", subtree_name);
  }

  if (feats != NULL) 
    msa_map_gff_coords(msa, feats, 1, 0, 0, NULL);
  /* NOTE: msa offset not currently handled */

  /* SPH method */
  if (method == SPH) {
    /* fit model to whole data set if necessary */
    if (fit_model && (!base_by_base && feats == NULL)) 
      mod_fitted = fit_tree_model(mod, msa, subtree_name, &scale, &sub_scale);

    /* set up for subtree mode */
    if (subtree_name != NULL) {
      if (!tm_is_reversible(mod->subst_mod))
        die("ERROR: reversible model required with --subtree.\n");
      tr_name_ancestors(mod->tree);
      sub_reroot(mod, subtree_name);
      if (mod_fitted != NULL) sub_reroot(mod_fitted, subtree_name);
      /* note: rerooting has to be done before creating jump process */
      if (fit_model && base_by_base) 
        mod->subtree_root = mod->tree->lchild; /* for rescaling */
    }

    /* if base-by-base, use larger default epsilon */
    if (base_by_base && default_epsilon)
      epsilon = DEFAULT_EPSILON_BASE_BY_BASE;

    /* jump process for prior */
    jp = sub_define_jump_process(mod, epsilon, tr_total_len(mod->tree));

    /* jump process for posterior -- use fitted model if necessary */
    if (mod_fitted != NULL)
      jp_post = sub_define_jump_process(mod_fitted, epsilon, 
                                        tr_max_branchlen(mod->tree));
    else if (fit_model && base_by_base) 
      jp_post = sub_define_jump_process(mod, epsilon, 
                                        10 * tr_max_branchlen(mod->tree));
    else
      jp_post = jp;
    
    if (nsites == -1) nsites = msa->length;

    /* now actually compute and print output */
    if (subtree_name == NULL) {   /* full-tree mode */
      if (base_by_base) {
        /* compute p-vals and (optionally) posterior means/variances per
           tuple, then print to stdout in wig or wig-like format */  
        pvals = smalloc(msa->ss->ntuples * sizeof(double));
        if (!output_wig) {
          post_means = smalloc(msa->ss->ntuples * sizeof(double));
          post_vars = smalloc(msa->ss->ntuples * sizeof(double));
        }
        sub_pval_per_site(jp, msa, mode, fit_model, &prior_mean, &prior_var, 
                          pvals, post_means, post_vars, logf);

        if (output_wig) 
          print_wig(msa, pvals, chrom, refidx, TRUE);
        else {
          char str[1000];
          sprintf(str, "#neutral mean = %.3f var = %.3f\n#post_mean post_var pval", 
                  prior_mean, prior_var);
          print_base_by_base(str, chrom, msa, NULL, refidx, 3, post_means, 
                             post_vars, pvals);
        }
      }
      else if (feats == NULL) {
        double post_mean, post_var;

        /* compute distributions and stats*/
        if (!post_only) 
          prior_distrib = sub_prior_distrib_alignment(jp, nsites);

        if (post_only)
          post_distrib = sub_posterior_distrib_alignment(jp_post, msa);
        else if (!prior_only) /* don't need explicit distrib for p-value */
          sub_posterior_stats_alignment(jp_post, msa, &post_mean, &post_var);

        /* print output */
        if (quantiles)
          print_quantiles(prior_only ? prior_distrib : post_distrib);
        else if (prior_only) 
          print_prior_only(nsites, argv[optind], prior_distrib);
        else if (post_only)
          print_post_only(argv[optind], argv[optind+1], post_distrib, ci, scale);
        else
          print_p(argv[optind], argv[optind+1], prior_distrib, 
                  post_mean, post_var, ci, scale);
      }
      else {                        /* --features case */
        p_value_stats *stats = sub_p_value_many(jp, msa, feats->features, ci);
        msa_map_gff_coords(msa, feats, 0, 1, 0, NULL);
        print_feats_sph(stats, feats, mode, epsilon, output_gff);
      }
    }
    else {			/* SPH and supertree/subtree */
      if (base_by_base) {
        /* compute p-vals and (optionally) posterior means/variances per
           tuple, then print to stdout in wig or wig-like format */  
        double *post_means_sub = NULL, *post_vars_sub = NULL, 
          *post_means_sup = NULL, *post_vars_sup = NULL;
        double prior_mean_sub, prior_var_sub, prior_mean_sup, prior_var_sup;
        pvals = smalloc(msa->ss->ntuples * sizeof(double));
        if (!output_wig) {
          post_means_sub = smalloc(msa->ss->ntuples * sizeof(double)); 
          post_means_sup = smalloc(msa->ss->ntuples * sizeof(double)); 
          post_vars_sub = smalloc(msa->ss->ntuples * sizeof(double)); 
          post_vars_sup = smalloc(msa->ss->ntuples * sizeof(double)); 
        }
        sub_pval_per_site_subtree(jp, msa, mode, fit_model, &prior_mean_sub, 
                                  &prior_var_sub, &prior_mean_sup, 
                                  &prior_var_sup, pvals, post_means_sub, 
                                  post_vars_sub, post_means_sup, post_vars_sup, 
                                  logf);

        if (output_wig) 
          print_wig(msa, pvals, chrom, refidx, TRUE);
        else {
          char str[1000];
          sprintf(str, "#neutral mean_sub = %.3f var_sub = %.3f mean_sup = %.3f  var_sup = %.3f\n#post_mean_sub post_var_sub post_mean_sup post_var_sup pval", 
                  prior_mean_sub, prior_var_sub, prior_mean_sup, prior_var_sup);
          print_base_by_base(str, chrom, msa, NULL, refidx, 5, post_means_sub, 
                             post_vars_sub, post_means_sup, post_vars_sup, 
                             pvals);
        }
      }

      else if (feats == NULL) {
        double post_mean, post_var, post_mean_sup, post_var_sup, 
          post_mean_sub, post_var_sub;

        /* compute distributions and stats */
        if (!post_only)
          prior_joint_distrib = sub_prior_joint_distrib_alignment(jp, nsites);

        if (post_only)
          post_joint_distrib = sub_posterior_joint_distrib_alignment(jp_post, msa);
        else if (!prior_only)
          sub_posterior_joint_stats_alignment(jp_post, msa, &post_mean, &post_var,
                                              &post_mean_sub, &post_var_sub, 
                                              &post_mean_sup, &post_var_sup);

        /* print output */
        if (prior_only) 
          print_prior_only_joint(subtree_name, nsites, argv[optind],
                                 prior_joint_distrib);
        else if (post_only) 
          print_post_only_joint(subtree_name, argv[optind], 
                                argv[optind+1], post_joint_distrib, ci,
                                scale, sub_scale);
        else 
          print_p_joint(subtree_name, argv[optind], argv[optind+1],
                        ci, prior_joint_distrib, post_mean, post_var, 
                        post_mean_sup, post_var_sup, post_mean_sub, 
                        post_var_sub, scale, sub_scale);
      }
      else {                      /* --features case */
        p_value_joint_stats *jstats = 
          sub_p_value_joint_many(jp, msa, feats->features, 
                                 ci, MAX_CONVOLVE_SIZE, NULL);
        msa_map_gff_coords(msa, feats, 0, 1, 0, NULL);
        print_feats_sph_subtree(jstats, feats, mode, epsilon, output_gff);
      }
    }
  } /* end SPH */

  /* LRT method */
  else if (method == LRT) {
    if (base_by_base) { 
      pvals = smalloc(msa->ss->ntuples * sizeof(double));
      if (!output_wig) {
        llrs = smalloc(msa->ss->ntuples * sizeof(double));
        scales = smalloc(msa->ss->ntuples * sizeof(double));
      }
      if (subtree_name == NULL) { /* no subtree case */
        col_lrts(mod, msa, mode, pvals, scales, llrs, logf);
        if (output_wig) 
          print_wig(msa, pvals, chrom, refidx, TRUE);
        else 
          print_base_by_base("#scale lnlratio pval", chrom, msa, NULL, refidx, 
                             3, scales, llrs, pvals);
      }
      else {                    /* subtree case */
        if (!output_wig) {
          sub_scales = smalloc(msa->ss->ntuples * sizeof(double));
          null_scales = smalloc(msa->ss->ntuples * sizeof(double));
        }

        col_lrts_sub(mod, msa, mode, pvals, null_scales, scales, sub_scales, 
                     llrs, logf);

        if (output_wig) 
          print_wig(msa, pvals, chrom, refidx, TRUE);
        else 
          print_base_by_base("#null_scale alt_scale alt_subscale lnlratio pval", 
                             chrom, msa, NULL, refidx, 5, null_scales, scales, 
                             sub_scales, llrs, pvals);
      }
    }
    else if (feats != NULL) {   /* feature-by-feature evaluation */
      pvals = smalloc(lst_size(feats->features) * sizeof(double));
      if (!output_gff) {
        scales = smalloc(lst_size(feats->features) * sizeof(double));
        llrs = smalloc(lst_size(feats->features) * sizeof(double));
      }
      if (subtree_name == NULL) {  /* no subtree case */
        ff_lrts(mod, msa, feats, mode, pvals, scales, llrs, logf);
        msa_map_gff_coords(msa, feats, 0, 1, 0, NULL);
        if (output_gff) 
          print_gff_scores(feats, pvals, TRUE);
        else
          print_feats_generic("scale\tlnlratio\tpval", feats, NULL, 3,
                              scales, llrs, pvals);
      }
      else {                    /* subtree case */
        if (!output_gff) {
          null_scales = smalloc(lst_size(feats->features) * sizeof(double));
          sub_scales = smalloc(lst_size(feats->features) * sizeof(double));
        }
        ff_lrts_sub(mod, msa, feats, mode, pvals, null_scales, scales, 
                    sub_scales, llrs, logf);
        msa_map_gff_coords(msa, feats, 0, 1, 0, NULL);
        if (output_gff) 
          print_gff_scores(feats, pvals, TRUE);
        else
          print_feats_generic("null_scale\talt_scale\talt_subscale\tlnlratio\tpval",
                              feats, NULL, 5, null_scales, scales, sub_scales,
                              llrs, pvals);
      }
    }
  } /* end LRT method */

  /* SCORE method */
  else if (method == SCORE) {
    if (base_by_base) {

      pvals = smalloc(msa->ss->ntuples * sizeof(double));
      if (!output_wig) {
        teststats = smalloc(msa->ss->ntuples * sizeof(double));
        derivs = smalloc(msa->ss->ntuples * sizeof(double));
      }

      if (subtree_name == NULL) { /* no subtree case */
        col_score_tests(mod, msa, mode, pvals, derivs, 
                        teststats);
        if (output_wig) 
          print_wig(msa, pvals, chrom, refidx, TRUE);
        else 
          print_base_by_base("#deriv teststat pval", chrom, msa, NULL, refidx, 
                             3, derivs, teststats, pvals);
      }
      else {                    /* subtree case */
        if (!output_wig) {
          null_scales = smalloc(msa->ss->ntuples * sizeof(double));
          sub_derivs = smalloc(msa->ss->ntuples * sizeof(double));
        }

        col_score_tests_sub(mod, msa, mode, pvals, null_scales, derivs, 
                            sub_derivs, teststats, logf);

        if (output_wig) 
          print_wig(msa, pvals, chrom, refidx, TRUE);
        else 
          print_base_by_base("#scale deriv subderiv teststat pval", chrom, 
                             msa, NULL, refidx, 5, null_scales, derivs, 
                             sub_derivs, teststats, pvals);
      }
    }
    else if (feats != NULL) {   /* feature by feature evaluation */
      pvals = smalloc(lst_size(feats->features) * sizeof(double));
      if (!output_gff) {
        teststats = smalloc(lst_size(feats->features) * sizeof(double));
        derivs = smalloc(lst_size(feats->features) * sizeof(double));
      }
      if (subtree_name == NULL) { /* no subtree case */
        ff_score_tests(mod, msa, feats, mode, pvals, derivs, teststats);
        msa_map_gff_coords(msa, feats, 0, 1, 0, NULL);
        if (output_gff) 
          print_gff_scores(feats, pvals, TRUE);
        else
          print_feats_generic("deriv\tteststat\tpval", feats, NULL, 3,
                              derivs, teststats, pvals);
      }
      else {                     /* subtree case */
        if (!output_gff) {
          null_scales = smalloc(lst_size(feats->features) * sizeof(double));
          sub_derivs = smalloc(lst_size(feats->features) * sizeof(double));
        }
        ff_score_tests_sub(mod, msa, feats, mode, pvals, null_scales, derivs, 
                           sub_derivs, teststats, logf);
        msa_map_gff_coords(msa, feats, 0, 1, 0, NULL);
        if (output_gff) 
          print_gff_scores(feats, pvals, TRUE);
        else
          print_feats_generic("scale\tderiv\tsubderiv\tteststat\tpval",
                              feats, NULL, 5, null_scales, derivs, sub_derivs,
                              teststats, pvals);
      }
    }
  } /* end SCORE */

  /* GERP method */
  else if (method == GERP) { 
    char *formatstr[4] = {"%.3f", "%.3f", "%.3f", "%.0f"};
    if (base_by_base) {
      nrejected = smalloc(msa->ss->ntuples * sizeof(double));
      if (!output_wig) {
        nneut = smalloc(msa->ss->ntuples * sizeof(double));
        nobs = smalloc(msa->ss->ntuples * sizeof(double));
        nspec = smalloc(msa->ss->ntuples * sizeof(double));
      }
      col_gerp(mod, msa, mode, nneut, nobs, nrejected, nspec, logf);
      if (output_wig) 
        print_wig(msa, nrejected, chrom, refidx, FALSE);
      else {
        print_base_by_base("#nneut nobs nrej nspec", chrom, msa, formatstr, 
                           refidx, 4, nneut, nobs, nrejected, nspec);
      }
    }
    else if (feats != NULL) {   /* feature by feature evaluation */
      nrejected = smalloc(lst_size(feats->features) * sizeof(double));
      if (!output_gff) {
        nneut = smalloc(lst_size(feats->features) * sizeof(double));
        nobs = smalloc(lst_size(feats->features) * sizeof(double));
        nspec = smalloc(lst_size(feats->features) * sizeof(double));
      }
      ff_gerp(mod, msa, feats, mode, nneut, nobs, nrejected, nspec, logf);
      msa_map_gff_coords(msa, feats, 0, 1, 0, NULL);
      if (output_gff) 
        print_gff_scores(feats, nrejected, FALSE);
      else 
        print_feats_generic("nneut\tnobs\tnrej\tnspec", feats, formatstr, 4,
                            nneut, nobs, nrejected, nspec);
    }
  } /* end GERP */
    
  return 0;
}

/* estimate scale parameters for model from data */
TreeModel* fit_tree_model(TreeModel *source_mod, MSA *msa, 
                          char *subtree_name, double *scale, 
                          double *sub_scale) {
  Vector *params;
  TreeModel *retval = tm_create_copy(source_mod);
  double oldscale;

  tm_free_rmp(retval);
  retval->estimate_branchlens = TM_SCALE_ONLY;

  if (subtree_name != NULL) {
    retval->subtree_root = tr_get_node(retval->tree, subtree_name);
    if (retval->subtree_root == NULL) 
      die("ERROR: no node named '%s'.\n", subtree_name);
    /* also make sure the supertree has nonzero branch length in the
       unrooted tree */
    if ((retval->tree->lchild == retval->subtree_root && retval->tree->rchild->lchild == NULL) || 
        (retval->tree->rchild == retval->subtree_root && retval->tree->lchild->lchild == NULL))
      die("ERROR: supertree contains no branches (in unrooted tree).\n");
  }

  retval->estimate_ratemat = FALSE;
  tm_init_rmp(source_mod);           /* (no. params changes) */
  params = tm_params_new_init_from_model(retval);

  tm_fit(retval, msa, params, -1, OPT_HIGH_PREC, NULL);

  oldscale = vec_get(params, 0);

  if (subtree_name == NULL) {
    /* correction for variance in estimates.  Based on simulation
       experiments, the observed standard deviation in estimated scale
       factors is about 4/3 that expected from true variation in the
       number of substitutions, across various element lengths.  If we
       assume Gaussian estimation error with mean zero, and we assume
       Gaussian numbers of substitutions (approximately true by CLT),
       then we can correct by scaling the difference from the mean by
       3/4.  If no correction is made, you see an excess of very small
       and very large p-values due to a thickening of the tails in the
       distribution of estimated numbers of substitutions  */
    *scale = (oldscale - 1) * 0.75 + 1; 
    tm_scale(retval, *scale/oldscale, 0);  /* (tree has already been
                                              scaled by oldscale) */    
  }
  else {
    /* no correction in subtree case: simulation experiments indicate
       that the method is conservative enough to compensate for any
       estimation error */
    *scale = oldscale;
    *sub_scale = vec_get(params, 1) * oldscale; 
                                /* actual scale of sub is product of
                                   overall scale and of
                                   subtree-specific parameter  */
  }

  vec_free(params);
  return retval;
}


