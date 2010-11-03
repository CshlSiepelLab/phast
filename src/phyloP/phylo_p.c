/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
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
#include "list_of_lists.h"
#include "phyloP.h"
#include "phylo_p.h"
#include "fit_column.h"
#include "fit_feature.h"
#include "trees.h"


/* initialize phyloP options to default (may be different for rphast) */
struct phyloP_struct *phyloP_struct_new(int rphast) {
  struct phyloP_struct *p = smalloc(sizeof(struct phyloP_struct));
  p->msa = NULL;
  p->nsites = -1;

  p->prior_only = FALSE;
  p->post_only = FALSE;
  p->quantiles_only = FALSE;
  
  p->output_wig = FALSE;
  p->output_gff = FALSE;

  p->fit_model = FALSE;
  p->base_by_base = FALSE;
  p->refidx = 1;
  p->refidx_feat = 1;
  p->ci = -1;
  p->epsilon = -1;
  p->subtree_name = NULL;
  p->chrom = NULL;
  p->branch_name = NULL;
  p->feats = NULL;
  p->method = SPH;
  p->mode = CON;
  p->outfile = rphast ? NULL : stdout;
  p->logf = NULL;
  p->mod = NULL;
  p->cats_to_do = NULL;
  p->cm = NULL;
  p->help = rphast ? "?phyloP" : "phyloP -h";
  p->mod_fname = NULL;
  p->msa_fname = NULL;

  p->results = rphast ? lol_new(20) : NULL;
  return p;
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
    if (retval->subtree_root == NULL)  {
      tr_name_ancestors(retval->tree);
      retval->subtree_root = tr_get_node(retval->tree, subtree_name);
      if (retval->subtree_root == NULL)
	die("ERROR: no node named '%s'.\n", subtree_name);
    }
    /* also make sure the supertree has nonzero branch length in the
       unrooted tree */
    if ((retval->tree->lchild == retval->subtree_root && retval->tree->rchild->lchild == NULL) || 
        (retval->tree->rchild == retval->subtree_root && retval->tree->lchild->lchild == NULL))
      die("ERROR: supertree contains no branches (in unrooted tree).\n");
  }

  retval->estimate_ratemat = FALSE;
  tm_init_rmp(source_mod);           /* (no. params changes) */
  params = tm_params_new_init_from_model(retval);

  tm_fit(retval, msa, params, -1, OPT_HIGH_PREC, NULL, 1);

  oldscale = vec_get(params, retval->scale_idx);

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
    tm_scale_branchlens(retval, *scale/oldscale, 0);  /* (tree has already been
                                              scaled by oldscale) */    
  }
  else {
    /* no correction in subtree case: simulation experiments indicate
       that the method is conservative enough to compensate for any
       estimation error */
    *scale = oldscale;
    *sub_scale = vec_get(params, retval->scale_idx+1) * oldscale; 
                                /* actual scale of sub is product of
                                   overall scale and of
                                   subtree-specific parameter  */
  }

  vec_free(params);
  return retval;
}


void phyloP(struct phyloP_struct *p) {
  /* variables for options that are passed through p */
  int nsites, fit_model, base_by_base, refidx;
  int prior_only, post_only, quantiles_only,
    output_wig, output_gff;
  double ci, epsilon;
  char *subtree_name, *chrom;
  List *branch_name;
  GFF_Set *feats;
  method_type method;
  mode_type mode;
  FILE *logf;
  List *cats_to_do;
  CategoryMap *cm;
  char *mod_fname, *msa_fname, *help;
  FILE *outfile;
  ListOfLists *results;

  /* other variables */
  TreeModel *mod, *mod_fitted = NULL;
  MSA *msa;
  Vector *prior_distrib=NULL, *post_distrib=NULL;
  Matrix *prior_joint_distrib=NULL, *post_joint_distrib=NULL;
  JumpProcess *jp, *jp_post;
  List *pruned_names;
  int j, old_nleaves;
  double scale = -1, sub_scale = -1;
  double prior_mean, prior_var;
  double *pvals = NULL, *post_means = NULL, *post_vars = NULL,
    *llrs = NULL, *scales = NULL, *sub_scales = NULL, 
    *null_scales = NULL, *derivs = NULL, *sub_derivs = NULL,
    *teststats = NULL, *nneut = NULL, *nobs = NULL,
    *nspec = NULL, *nrejected = NULL;

  msa = p->msa;
  nsites = p->nsites;
  fit_model = p->fit_model;
  base_by_base = p->base_by_base;
  refidx = p->refidx;
  ci = p->ci;
  epsilon = p->epsilon;
  subtree_name = p->subtree_name;
  chrom = p->chrom;
  branch_name = p->branch_name;
  feats = p->feats;
  method = p->method;
  mode = p->mode;
  logf = p->logf;
  mod = p->mod;
  cats_to_do = p->cats_to_do;
  cm = p->cm;
  help = p->help;
  mod_fname = p->mod_fname;
  msa_fname = p->msa_fname;
  outfile = p->outfile;
  results = p->results;

  prior_only = p->prior_only;
  post_only = p->post_only;
  quantiles_only = p->quantiles_only;
  output_wig = p->output_wig;
  output_gff = p->output_gff;

  if (msa == NULL && !prior_only)
    die("Need either --prior-only or an alignment\n");
  if (msa != NULL && (refidx < 0 || refidx > msa->nseqs))
    die("refidx should be 0 (for alignment frame of refernce), or between 1 and msa->nseqs (%i)", msa->nseqs);
  if (chrom == NULL) {
    if (msa != NULL && refidx >= 1) chrom = msa->names[refidx-1];
    else chrom="chr1";
  }
  if (method != SPH && (fit_model || epsilon>=0 || ci!=-1 || 
			prior_only || post_only || quantiles_only))
    die("ERROR: given arguments only available in SPH mode.  Try '%s'.\n", help);
  if (quantiles_only && !prior_only && !post_only)
    die("ERROR: --quantiles can only be used with --null or --posterior.\n");
  if (subtree_name != NULL && quantiles_only)
    die("ERROR: --quantiles cannot be used with --subtree.\n");
  if (feats != NULL && (fit_model || prior_only || post_only))
    die("ERROR: --features cannot be used with --null, --posterior, or --fit-model.\n");
  if (base_by_base && (ci != -1 || feats!=NULL || prior_only || post_only))
    die("ERROR: --wig-scores and --base-by-base cannot be used with --null, --posterior, --features, --quantiles, or --confidence-interval.\n");
  if (method == GERP && subtree_name != NULL)
    die("ERROR: --subtree not supported with --method GERP.\n");
  if ((method == GERP || method == SPH) && branch_name != NULL) 
    die("ERROR --branch not supported with --method GERP or --method SPH\n");
  if (branch_name != NULL && subtree_name != NULL)
    die("ERROR: can use only one of --subtree or --branch options\n");
  if (method != SPH && feats==NULL && !base_by_base)
    die("ERROR: need base-by-base, wig-scores, or features unless method is SPH\n");
  if (prior_only && msa==NULL && nsites < 0)
    die("ERROR: need to specify nsites or msa to get prior");
  if (!prior_only) {
    if (msa->ss == NULL)
      ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1, 0);

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
      String *warnstr = str_new(1000);
      str_cpy_charstr(warnstr, "WARNING: pruned away leaves with no match in alignment (");
      for (j = 0; j < lst_size(pruned_names); j++) {
	str_append(warnstr, (String*)lst_get_ptr(pruned_names, j));
	str_append_charstr(warnstr, j < lst_size(pruned_names) ? ", " : ").\n");
      }
      phast_warning(warnstr->chars);
      str_free(warnstr);
    }
    lst_free(pruned_names);
  }

  /* set subtree if necessary */
  if (subtree_name != NULL && method != SPH) {
    /* (SPH is a special case -- requires rerooting) */
    mod->subtree_root = tr_get_node(mod->tree, subtree_name);
    if (mod->subtree_root == NULL) {
      tr_name_ancestors(mod->tree);
      mod->subtree_root = tr_get_node(mod->tree, subtree_name);
      if (mod->subtree_root == NULL)
	die("ERROR: no node named '%s'.\n", subtree_name);
    }
  }
  if (branch_name != NULL) {
    TreeNode *n;
    char *nodeName;
    tr_name_ancestors(mod->tree);
    mod->in_subtree = smalloc(mod->tree->nnodes * sizeof(int));
    for (j=0; j<mod->tree->nnodes; j++)
      mod->in_subtree[j] = 0;
    for (j=0; j<lst_size(branch_name); j++) {
      nodeName = ((String*)lst_get_ptr(branch_name, j))->chars;
      n = tr_get_node(mod->tree, nodeName);
      if (n == NULL) {
	tr_name_ancestors(mod->tree);
	n = tr_get_node(mod->tree, nodeName);
	if (n == NULL) {
	  die("ERROR: no node named %s\n", nodeName);
	}
      }
      mod->in_subtree[n->id] = 1;
    }
    for (j=0; j<mod->tree->nnodes; j++)
      if (mod->in_subtree[j] == 0) break;
    if (j == mod->tree->nnodes)
      die("ERROR: ERROR: cannot name all branches with --branch option\n");
  }

  if (feats != NULL) {
    if (msa->idx_offset > 0)
      gff_add_offset(feats, -(msa->idx_offset), msa_seqlen(msa, 0));
    msa_map_gff_coords(msa, feats, p->refidx_feat, 0, 0, NULL);
  }

  /* SPH method */
  if (method == SPH) {
    /* fit model to whole data set if necessary */
    if (fit_model && (!base_by_base && feats == NULL)) 
      mod_fitted = fit_tree_model(mod, msa, subtree_name, &scale, &sub_scale);

    /* set up for subtree mode */
    if (subtree_name != NULL) {
      if (!tm_is_reversible(mod))
        die("ERROR: reversible model required with --subtree.\n");
      tr_name_ancestors(mod->tree);
      sub_reroot(mod, subtree_name);
      if (mod_fitted != NULL) sub_reroot(mod_fitted, subtree_name);
      /* note: rerooting has to be done before creating jump process */
      if (fit_model && base_by_base) 
        mod->subtree_root = mod->tree->lchild; /* for rescaling */
    }
    /* if base-by-base, use larger default epsilon */
    if (base_by_base && epsilon < 0)
      epsilon = DEFAULT_EPSILON_BASE_BY_BASE;
    else if (epsilon < 0) epsilon = DEFAULT_EPSILON;
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
	if (results != NULL || !output_wig) {
          post_means = smalloc(msa->ss->ntuples * sizeof(double));
          post_vars = smalloc(msa->ss->ntuples * sizeof(double));
        }
        sub_pval_per_site(jp, msa, mode, fit_model, &prior_mean, &prior_var, 
                          pvals, post_means, post_vars, logf);

        if (outfile != NULL && output_wig)
          print_wig(outfile, msa, pvals, chrom, refidx, TRUE, NULL);
	if ((outfile != NULL && !output_wig) || results!=NULL) {
	  char str[1000];
	  sprintf(str, "#neutral mean = %.3f var = %.3f\n#post_mean post_var pval", 
                  prior_mean, prior_var);
          print_base_by_base(output_wig ? NULL : outfile, 
			     str, chrom, msa, NULL, 
			     refidx, results, FALSE, TRUE, 3,
			     "post.mean", post_means, 
                             "post.var", post_vars, "pval", pvals);
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
        if (quantiles_only)
          print_quantiles(outfile, prior_only ? prior_distrib : post_distrib,
			  results);
        else if (prior_only) 
          print_prior_only(outfile, nsites, mod_fname, prior_distrib, results);
        else if (post_only)
          print_post_only(outfile, mod_fname, msa_fname, post_distrib, 
			  ci, scale, results);
        else
          print_p(outfile, mod_fname, msa_fname, prior_distrib, 
                  post_mean, post_var, ci, scale, results);
	if (prior_distrib != NULL) vec_free(prior_distrib);
      }
      else {                        /* --features case */
        p_value_stats *stats = sub_p_value_many(jp, msa, feats->features, ci);
        msa_map_gff_coords(msa, feats, 0, p->refidx_feat, 0, NULL);
	if (msa->idx_offset > 0)
	  gff_add_offset(feats, msa->idx_offset, 0);
        print_feats_sph(outfile, stats, feats, mode, epsilon, output_gff, 
			results);
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
        if (results != NULL || !output_wig) {
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
          print_wig(outfile, msa, pvals, chrom, refidx, TRUE, results);
	if (results != NULL || !output_wig) {
          char str[1000];
          sprintf(str, "#neutral mean_sub = %.3f var_sub = %.3f mean_sup = %.3f  var_sup = %.3f\n#post_mean_sub post_var_sub post_mean_sup post_var_sup pval", 
                  prior_mean_sub, prior_var_sub, prior_mean_sup, prior_var_sup);
          print_base_by_base(output_wig ? NULL : outfile, 
			     str, chrom, msa, NULL, 
			     refidx, results, FALSE, TRUE, 5, 
			     "post.mean.sub", post_means_sub, 
                             "post.var.sub", post_vars_sub, 
			     "post.mean.sup", post_means_sup, 
			     "post.var.sup", post_vars_sup, 
                             "pval", pvals);
	  sfree(post_means_sub);
	  sfree(post_means_sup);
	  sfree(post_vars_sub);
	  sfree(post_vars_sup);
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
          print_prior_only_joint(outfile, subtree_name, nsites, mod_fname,
                                 prior_joint_distrib, results);
        else if (post_only) 
          print_post_only_joint(outfile, subtree_name, mod_fname, 
                                msa_fname, post_joint_distrib, ci,
                                scale, sub_scale, results);
        else 
          print_p_joint(outfile, subtree_name, mod_fname, msa_fname,
                        ci, prior_joint_distrib, post_mean, post_var, 
                        post_mean_sup, post_var_sup, post_mean_sub, 
                        post_var_sub, scale, sub_scale, results);
      }
      else {                      /* --features case */
        p_value_joint_stats *jstats = 
          sub_p_value_joint_many(jp, msa, feats->features, 
                                 ci, MAX_CONVOLVE_SIZE, NULL);
        msa_map_gff_coords(msa, feats, 0, p->refidx_feat, 0, NULL);
	if (msa->idx_offset > 0)
	  gff_add_offset(feats, msa->idx_offset, 0);
        print_feats_sph_subtree(outfile, jstats, feats, mode, epsilon, 
				output_gff, results);
      }
    }
  } /* end SPH */

  /* LRT method */
  else if (method == LRT) {
    if (base_by_base) { 
      pvals = smalloc(msa->ss->ntuples * sizeof(double));
      if (results != NULL || !output_wig) {
        llrs = smalloc(msa->ss->ntuples * sizeof(double));
        scales = smalloc(msa->ss->ntuples * sizeof(double));
      }
      if (subtree_name == NULL && branch_name == NULL) { /* no subtree case */
        col_lrts(mod, msa, mode, pvals, scales, llrs, logf);
        if (output_wig) 
          print_wig(outfile, msa, pvals, chrom, refidx, TRUE, NULL);
	if (results != NULL || !output_wig)
          print_base_by_base(output_wig ? NULL : outfile, 
			     "#scale lnlratio pval", 
			     chrom, msa, NULL, refidx, results, FALSE, TRUE, 3, 
                             "scale", scales, "lnlratio", llrs, "pval", pvals);
      }
      else {                    /* subtree case */
        if (results != NULL || !output_wig) {
          sub_scales = smalloc(msa->ss->ntuples * sizeof(double));
          null_scales = smalloc(msa->ss->ntuples * sizeof(double));
        }
        col_lrts_sub(mod, msa, mode, pvals, null_scales, scales, sub_scales, 
                     llrs, logf);

        if (output_wig) 
          print_wig(outfile, msa, pvals, chrom, refidx, TRUE, NULL);
	if (results != NULL || !output_wig)
          print_base_by_base(output_wig ? NULL : outfile, 
			     "#null_scale alt_scale alt_subscale lnlratio pval", 
                             chrom, msa, NULL, refidx, results, FALSE, TRUE, 5, 
			     "null.scale", null_scales, "alt.scale", scales,
                             "alt.subscale", sub_scales, "lnlratio", llrs,
			     "pval", pvals);
      }
    }
    else if (feats != NULL) {   /* feature-by-feature evaluation */
      pvals = smalloc(lst_size(feats->features) * sizeof(double));
      if (results != NULL || !output_gff) {
        scales = smalloc(lst_size(feats->features) * sizeof(double));
        llrs = smalloc(lst_size(feats->features) * sizeof(double));
      }
      if (subtree_name == NULL && branch_name == NULL) {  /* no subtree case */
        ff_lrts(mod, msa, feats, mode, pvals, scales, llrs, logf);
        msa_map_gff_coords(msa, feats, 0, p->refidx_feat, 0, NULL);
	if (msa->idx_offset > 0)
	  gff_add_offset(feats, msa->idx_offset, 0);
        if (output_gff) 
          print_gff_scores(outfile, feats, pvals, TRUE);
	if (results != NULL || !output_gff)
          print_feats_generic(output_gff ? NULL : outfile, 
			      "scale\tlnlratio\tpval",
			      feats, NULL, results, FALSE, TRUE, 3, 
                              "scale", scales, "lnlratio", llrs, "pval", pvals);
      }
      else {                    /* subtree case */
        if (results != NULL || !output_gff) {
          null_scales = smalloc(lst_size(feats->features) * sizeof(double));
          sub_scales = smalloc(lst_size(feats->features) * sizeof(double));
        }
        ff_lrts_sub(mod, msa, feats, mode, pvals, null_scales, scales, 
                    sub_scales, llrs, logf);
        msa_map_gff_coords(msa, feats, 0, p->refidx_feat, 0, NULL);
	if (msa->idx_offset > 0)
	  gff_add_offset(feats, msa->idx_offset, 0);
        if (output_gff) 
          print_gff_scores(outfile, feats, pvals, TRUE);
	if (results != NULL || !output_gff)
          print_feats_generic(output_gff ? NULL : outfile, 
			      "null_scale\talt_scale\talt_subscale\tlnlratio\tpval",
                              feats, NULL, results, FALSE, TRUE, 5, 
			      "null.scale", null_scales, "alt.scale", scales, 
			      "alt.subscale", sub_scales,"lnlratio", llrs,
                              "pval", pvals);
      }
    }
  } /* end LRT method */

  /* SCORE method */
  else if (method == SCORE) {
    if (base_by_base) {
      pvals = smalloc(msa->ss->ntuples * sizeof(double));
      if (results != NULL || !output_wig) {
        teststats = smalloc(msa->ss->ntuples * sizeof(double));
        derivs = smalloc(msa->ss->ntuples * sizeof(double));
      }

      if (subtree_name == NULL && branch_name == NULL) { /* no subtree case */
        col_score_tests(mod, msa, mode, pvals, derivs, 
                        teststats);
        if (output_wig) 
          print_wig(outfile, msa, pvals, chrom, refidx, TRUE, NULL);
	if (results != NULL || !output_wig)
          print_base_by_base(output_wig ? NULL : outfile, 
			     "#deriv teststat pval", 
			     chrom, msa, NULL, refidx, results, FALSE, TRUE, 3, 
                             "deriv", derivs, "teststat", teststats, 
			     "pval", pvals);
      }
      else {                    /* subtree case */
        if (results != NULL || !output_wig) {
          null_scales = smalloc(msa->ss->ntuples * sizeof(double));
          sub_derivs = smalloc(msa->ss->ntuples * sizeof(double));
        }

        col_score_tests_sub(mod, msa, mode, pvals, null_scales, derivs, 
                            sub_derivs, teststats, logf);

        if (output_wig) 
          print_wig(outfile, msa, pvals, chrom, refidx, TRUE, NULL);
	if (results != NULL || !output_wig)
          print_base_by_base(output_wig ? NULL : outfile, 
			     "#scale deriv subderiv teststat pval", 
			     chrom, msa, NULL, refidx, results, FALSE, TRUE, 5, 
			     "scale", null_scales, "deriv", derivs, 
                             "subderiv", sub_derivs, "teststat", teststats, 
			     "pval", pvals);
      }
    }
    else if (feats != NULL) {   /* feature by feature evaluation */
      pvals = smalloc(lst_size(feats->features) * sizeof(double));
      if (results != NULL || !output_gff) {
        teststats = smalloc(lst_size(feats->features) * sizeof(double));
        derivs = smalloc(lst_size(feats->features) * sizeof(double));
      }
      if (subtree_name == NULL && branch_name == NULL) { /* no subtree case */
        ff_score_tests(mod, msa, feats, mode, pvals, derivs, teststats);
        msa_map_gff_coords(msa, feats, 0, p->refidx_feat, 0, NULL);
	if (msa->idx_offset > 0)
	  gff_add_offset(feats, msa->idx_offset, 0);
        if (output_gff) 
          print_gff_scores(outfile, feats, pvals, TRUE);
	if (results != NULL || !output_gff)
          print_feats_generic(output_gff ? NULL : outfile, 
			      "deriv\tteststat\tpval", feats, NULL, 
			      results, FALSE, TRUE, 3,
                              "deriv", derivs, "teststat", teststats, 
			      "pval", pvals);
      }
      else {                     /* subtree case */
        if (results != NULL || !output_gff) {
          null_scales = smalloc(lst_size(feats->features) * sizeof(double));
          sub_derivs = smalloc(lst_size(feats->features) * sizeof(double));
        }
        ff_score_tests_sub(mod, msa, feats, mode, pvals, null_scales, derivs, 
                           sub_derivs, teststats, logf);
        msa_map_gff_coords(msa, feats, 0, p->refidx_feat, 0, NULL);
	if (msa->idx_offset > 0)
	  gff_add_offset(feats, msa->idx_offset, 0);
        if (output_gff) 
          print_gff_scores(outfile, feats, pvals, TRUE);
	if (results != NULL || !output_gff)
          print_feats_generic(output_gff ? NULL : outfile, 
			      "scale\tderiv\tsubderiv\tteststat\tpval",
                              feats, NULL, results, FALSE, TRUE, 5, 
			      "scale", null_scales, "deriv", derivs, 
			      "subderiv", sub_derivs, "teststat", teststats,
                              "pval", pvals);
      }
    }
  } /* end SCORE */

  /* GERP method */
  else if (method == GERP) { 
    char *formatstr[4] = {"%.3f", "%.3f", "%.3f", "%.0f"};
    if (base_by_base) {
      nrejected = smalloc(msa->ss->ntuples * sizeof(double));
      if (results != NULL || !output_wig) {
        nneut = smalloc(msa->ss->ntuples * sizeof(double));
        nobs = smalloc(msa->ss->ntuples * sizeof(double));
        nspec = smalloc(msa->ss->ntuples * sizeof(double));
      }
      col_gerp(mod, msa, mode, nneut, nobs, nrejected, nspec, logf);
      if (output_wig) 
        print_wig(outfile, msa, nrejected, chrom, refidx, FALSE, NULL);
      if (results != NULL || !output_wig) {
        print_base_by_base(output_wig ? NULL : outfile, 
			   "#nneut nobs nrej nspec", chrom, 
			   msa, formatstr, refidx, results, FALSE, FALSE, 4, 
			   "nneut", nneut, "nobs", nobs, "nrej", nrejected, 
			   "nspec", nspec);
      }
    }
    else if (feats != NULL) {   /* feature by feature evaluation */
      nrejected = smalloc(lst_size(feats->features) * sizeof(double));
      if (results != NULL || !output_gff) {
        nneut = smalloc(lst_size(feats->features) * sizeof(double));
        nobs = smalloc(lst_size(feats->features) * sizeof(double));
        nspec = smalloc(lst_size(feats->features) * sizeof(double));
      }
      ff_gerp(mod, msa, feats, mode, nneut, nobs, nrejected, nspec, logf);
      msa_map_gff_coords(msa, feats, 0, p->refidx_feat, 0, NULL);
      if (msa->idx_offset > 0)
	gff_add_offset(feats, msa->idx_offset, 0);
      if (output_gff) 
        print_gff_scores(outfile, feats, nrejected, FALSE);
      if (results != NULL || !output_gff)
        print_feats_generic(output_gff ? NULL : outfile, 
			    "nneut\tnobs\tnrej\tnspec", 
			    feats, formatstr, results, FALSE, TRUE, 4,
                            "nneut", nneut, "nobs", nobs, 
			    "nrej", nrejected, "nspec", nspec);
    }
  } /* end GERP */
  if (pvals != NULL) sfree(pvals);
  if (post_means != NULL) sfree(post_means);
  if (post_vars != NULL) sfree(post_vars);
  if (llrs != NULL) sfree(llrs);
  if (scales != NULL) sfree(scales);
  if (sub_scales != NULL) sfree(sub_scales);
  if (null_scales != NULL) sfree(null_scales);
  if (teststats != NULL) sfree(teststats);
  if (derivs != NULL) sfree(derivs);
  if (sub_derivs != NULL) sfree(sub_derivs);
  if (nrejected != NULL) sfree(nrejected);
  if (nneut != NULL) sfree(nneut);
  if (nobs != NULL) sfree(nobs);
  if (nspec != NULL) sfree(nspec);
} 


