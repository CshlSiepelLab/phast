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
#include <phast/misc.h>
#include <phast/msa.h>
#include <phast/maf.h>
#include <phast/nj.h>
#include <phast/tree_model.h>
#include <phast/subst_mods.h>
#include <phast/sufficient_stats.h>
#include "nj_var.help"

#define DEFAULT_NSAMPLES 100
#define DEFAULT_DIM 3
#define DEFAULT_BATCHSIZE 20
#define DEFAULT_LEARNRATE 0.001
#define DEFAULT_NBATCHES_CONV 10
#define DEFAULT_MIN_NBATCHES 30
#define DEFAULT_KAPPA 4

int main(int argc, char *argv[]) {
  signed char c;
  int opt_idx, i, ntips = 0, nsamples = DEFAULT_NSAMPLES, dim = DEFAULT_DIM,
    batchsize = DEFAULT_BATCHSIZE, nbatches_conv = DEFAULT_NBATCHES_CONV,
    min_nbatches = DEFAULT_MIN_NBATCHES;
  unsigned int nj_only = FALSE, random_start = FALSE;
  MSA *msa = NULL;

  char *alphabet = "ACGT";
  char **names = NULL;
  msa_format_type format = UNKNOWN_FORMAT;
  FILE *infile = NULL, *indistfile = NULL, *outdistfile = NULL, *logfile = NULL,
    *postmeanfile = NULL;
  Matrix *D = NULL;
  TreeNode *tree;
  List *namestr, *trees;
  subst_mod_type subst_mod = JC69;
  TreeModel *mod = NULL;
  double kappa = DEFAULT_KAPPA, learnrate = DEFAULT_LEARNRATE;
  MarkovMatrix *rmat = NULL;
  Vector *mu = NULL;
  Matrix *sigma = NULL;
  TreeNode *init_tree = NULL;

  struct option long_opts[] = {
    {"format", 1, 0, 'i'},
    {"batchsize", 1, 0, 'b'},
    {"nbatches-conv", 1, 0, 'c'},
    {"distances", 1, 0, 'd'},
    {"dimensionality", 1, 0, 'D'},
    {"hky85", 0, 0, 'k'}, 
    {"logfile", 1, 0, 'l'},
    {"mean", 1, 0, 'm'},
    {"min-nbatches", 1, 0, 'M'},
    {"names", 1, 0, 'n'},
    {"nj-only", 0, 0, 'j'},
    {"out-dists", 1, 0, 'o'},
    {"nsamples", 1, 0, 's'},
    {"learnrate", 1, 0, 'r'},
    {"random-start", 0, 0, 'R'},
    {"tree", 1, 0, 't'},
    {"treemodel", 1, 0, 'T'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "b:c:d:D:hi:jkl:m:M:n:o:r:Rt:T:s:", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'b':
      batchsize = atoi(optarg);
      if (batchsize <= 0)
        die("ERROR: --batchsize must be positive\n");
      break;
    case 'c':
      nbatches_conv = atoi(optarg);
      if (nbatches_conv <= 0)
        die("ERROR: --nbatches-conv must be positive\n");
      break;
    case 'd':
      indistfile = phast_fopen(optarg, "r");
      break;
    case 'D':
      dim = atoi(optarg);
      if (dim <= 0)
        die("ERROR: --dimensionality must be positive\n");
      break;
    case 'i':
      format = msa_str_to_format(optarg);
      if (format == UNKNOWN_FORMAT) 
        die("ERROR: bad input format.\n");
      break;
    case 'j':
      nj_only = TRUE;
      break;
    case 'k':
      subst_mod = HKY85;
      break;
    case 'l':
      logfile = phast_fopen(optarg, "w");
      break;
    case 'm':
      postmeanfile = phast_fopen(optarg, "w");
      break;
    case 'M':
      min_nbatches = atoi(optarg);
      if (min_nbatches <= 0)
        die("ERROR: --min-nbatches must be positive\n");
      break;
    case 'n':
      namestr = get_arg_list(optarg);
      ntips = lst_size(namestr);
      names = smalloc(ntips * sizeof(char*));
      for (i = 0; i < ntips; i++)
        names[i] = ((String*)lst_get_ptr(namestr, i))->chars;
      break;
    case 'o':
      outdistfile = phast_fopen(optarg, "w");
      break;
    case 'r':
      learnrate = atof(optarg);
      if (learnrate <= 0)
        die("ERROR: --learnrate must be positive\n");
      break;
    case 'R':
      random_start = TRUE;
      break;
    case 's':
      nsamples = atoi(optarg);
      if (nsamples <= 0)
        die("ERROR: --nsamples must be > 0\n");
      break;
    case 't':
      init_tree = tr_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'T':
      mod = tm_new_from_file(phast_fopen(optarg, "r"), 1);
      init_tree = mod->tree;
      break;
    case 'h':
      printf("%s", HELP); 
      exit(0);
    case '?':
      die("Bad argument.  Try 'nj_var -h'.\n");
    }
  }

  if (init_tree != NULL && indistfile != NULL)
    die("Cannot specify both --tree/-treemod and --distances\n");

  if (nj_only && (indistfile != NULL || init_tree != NULL))
    fprintf(stderr, "No alignment needed....\n");
  else {
    if (optind != argc - 1)
      die("ERROR: alignment file required.\n");
    
    infile = phast_fopen(argv[optind], "r");
    if (format == UNKNOWN_FORMAT)
      format = msa_format_for_content(infile, 1);
    if (format == MAF) 
      msa = maf_read(phast_fopen(argv[optind], "r"), NULL, 1, alphabet,
		     NULL, NULL, -1, TRUE, NULL, NO_STRIP, FALSE);
    else
      msa = msa_new_from_file_define_format(phast_fopen(argv[optind], "r"), 
					    format, alphabet);

    if (msa->ss == NULL)
      ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1, 0);
  }

    
  /* get a distance matrix */
  if (init_tree != NULL)
    D = nj_tree_to_distances(init_tree);
  else if (indistfile != NULL) {
    D = mat_new_from_file(indistfile, ntips, ntips);
    if (names == NULL || ntips <= 0)
      die("ERROR: --names required with --distances.  Try 'nj_var -h'.\n");
  }
  else if (msa != NULL)
    D = nj_compute_JC_matr(msa);
  else
    die("ERROR: no distance matrix available\n");
  
  /* make sure valid */
  nj_test_D(D);

  /* we'll need a tree in any case */
  if (names == NULL && msa != NULL)
    names = msa->names;
  tree = nj_infer_tree(D, names);

  if (nj_only == TRUE) /* just print in this case */
    tr_print(stdout, tree, TRUE);

  else {  /* full variational inference */
    if (msa == NULL)
      die("ERROR: Alignment required for variational inference\n");

    /* set up a tree model if necessary */
    if (mod == NULL) {
      rmat = mm_new(strlen(msa->alphabet), msa->alphabet, CONTINUOUS);    
      mod = tm_new(tree, rmat, NULL, subst_mod, msa->alphabet, 1, 1, NULL, -1);
      tm_init_backgd(mod, msa, -1); 
    
      if (subst_mod == JC69)
	tm_set_JC69_matrix(mod);
      else 
	tm_set_HKY_matrix(mod, kappa, -1);   /* FIXME: estimate kappa from msa */
    }

    /* initialize parameters of multivariate normal */    
    mu = vec_new(msa->nseqs * dim);
    sigma = mat_new(msa->nseqs * dim, msa->nseqs * dim);
    if (random_start == TRUE) {
      nj_sample_std_mvn(mu); vec_scale(mu, 0.1);
      mat_set_identity(sigma); mat_scale(sigma, 0.1);
    }
    else
      nj_estimate_mvn_from_distances(D, dim, mu, sigma);

    nj_variational_inf(mod, msa, D, mu, sigma, dim, batchsize, learnrate,
		       nbatches_conv, min_nbatches, logfile);
    trees = nj_var_sample(nsamples, dim, mu, sigma, msa->names);
    for (i = 0; i < nsamples; i++)
      tr_print(stdout, (TreeNode*)lst_get_ptr(trees, i), TRUE);

    if (postmeanfile != NULL)
      tr_print(postmeanfile, nj_mean(mu, dim, msa->names), TRUE);
  }

  if (outdistfile != NULL) {
    if (nj_only == FALSE)
      nj_points_to_distances(mu, D);  /* reset D to posterior mean */
    mat_print(D, outdistfile);
  }
    
  return (0);
}
