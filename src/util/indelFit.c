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
#include <tree_model.h>
#include <msa.h>
#include <gff.h>
#include <indel_history.h>
#include <indel_mod.h>
#include "indelFit.help"

int *get_cats(IndelHistory *ih, GFF_Set *feats, CategoryMap *cm,
              char *reference);

int main(int argc, char *argv[]) {
  char c;
  int opt_idx, i, cat;
  TreeNode *tree;
  IndelHistory *ih;
  IndelModel *im;
  IndelSuffStats *ss;
  double lnl = INFTY;
  int ncats = 1;
  CategoryMap *cm = NULL;
  int *cats = NULL;

  /* variables for options with defaults */
  double alpha = 0.02, beta = 0.04, tau = 0.05;
  int lnl_only = FALSE, columns = FALSE;
  FILE *logf = NULL;
  GFF_Set *feats = NULL;
  char *reference = NULL;

  struct option long_opts[] = {
    {"alpha", 1, 0, 'a'},
    {"beta", 1, 0, 'b'},
    {"tau", 1, 0, 't'},
    {"lnl", 0, 0, 'L'},
    {"columns", 0, 0, 'c'},
    {"features", 1, 0, 'f'},
    {"reference", 1, 0, 'r'},
    {"log", 1, 0, 'l'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "a:b:t:Lcf:r:l:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'a':
      alpha = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'b':
      beta = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 't':
      tau = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'L':
      lnl_only = TRUE;
      break;
    case 'c':
      columns = TRUE;
      break;
    case 'f':
      feats = gff_read_set(phast_fopen(optarg, "r"));
      break;
    case 'r':
      reference = optarg;
      break;
    case 'l':
      logf = phast_fopen(optarg, "w+");
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'indelFit -h'.\n");
    }
  }

  if (optind != argc - 2) 
    die("ERROR: Two arguments required.  Try 'indelFit -h'.\n");
  else if (lnl_only && logf != NULL)
    die("WARNING: --log ignored.\n");
  if (feats != NULL && (columns || lnl_only))
    die("ERROR: can't use --features with --lnl or --columns.\n");

  set_seed(-1);

  ih = ih_new_from_file(phast_fopen(argv[optind], "r"));
  tree = tr_new_from_file(phast_fopen(argv[optind+1], "r"));

  /* ensure trees are compatible */
  if (tree->nnodes != ih->tree->nnodes)
    die("ERROR: trees for indel model and indel history don't match.\n");
  for (i = 0; i < tree->nnodes; i++) {
    TreeNode *n1 = lst_get_ptr(tree->nodes, i);
    TreeNode *n2 = lst_get_ptr(ih->tree->nodes, i);
    if (n1->name[0] != '\0' && n2->name[0] != '\0' && 
        strcmp(n1->name, n2->name) != 0)
      die("ERROR: trees for indel model and indel history don't match.\n");
  }

  if (feats != NULL) {
    cm = cm_new_from_features(feats);
    cats = get_cats(ih, feats, cm, reference);
    ncats = cm->ncats + 1;
  }

  for (cat = 0; cat < ncats; cat++) {
    im = im_new_all(alpha, beta, tau, tree);
    ss = ncats == 1 ? im_suff_stats(ih) : im_suff_stats_cat(ih, cats, cat);

    if (lnl_only) 
      lnl = im_likelihood(im, ss);

    else {
      im_estimate(im, ih, ss, logf);
      lnl = im->training_lnl;
    }

    if (columns) {
      double *col_logl = smalloc(ih->ncols * sizeof(double));
      lnl = log(2) * im_column_logl(ih, im, col_logl);
      printf("#alpha = %f, beta = %f, tau = %f\n", im->alpha, im->beta, 
             im->tau);
      printf("#total lnl = %f\n", lnl);
      printf("#pos lnl\n");
      for (i = 0; i < ih->ncols; i++)
        printf("%d\t%f\n", i, col_logl[i] * log(2));
    }
    else {
      if (ncats > 1)
        printf("Category %d (%s): ", cat, cm_get_feature(cm, cat)->chars);
      printf("alpha = %f, beta = %f, tau = %f, lnl = %f\n", im->alpha, 
             im->beta, im->tau, lnl);
    }
  }

  return 0;
}

/* get array of categories for indel history based on given
   feature set */
int *get_cats(IndelHistory *ih, GFF_Set *feats, CategoryMap *cm,
              char *reference) {
  int *retval;
  int i;
  TreeNode *node;
  char *seq = smalloc(ih->ncols * sizeof(char));
  MSA *dummy_msa;

  if (reference != NULL) {
    node = tr_get_node(ih->tree, reference);
    if (node == NULL)
      die("ERROR: node '%s' not found in tree.\n", reference);

    /* make a dummy MSA based on the indel history */  
    for (i = 0; i < ih->ncols; i++) {
      if (ih->indel_strings[node->id][i] == BASE) seq[i] = 'A';
      else seq[i] = GAP_CHAR;
    }
    dummy_msa = msa_new(&seq, NULL, 1, ih->ncols, NULL);

    /* map features to alignment space */
    msa_map_gff_coords(dummy_msa, feats, 1, 0, 0);
  }
  else 
    /* (no need to map features, no need for seq) */
    dummy_msa = msa_new(NULL, NULL, 1, ih->ncols, NULL);

  /* label categories based on features */
  msa_label_categories(dummy_msa, feats, cm);

  /* pull out categories and free alignment */
  retval = dummy_msa->categories;
  dummy_msa->categories = NULL;
  dummy_msa->seqs = NULL;
  msa_free(dummy_msa);
  sfree(seq);

  return retval;
}
