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
#include <phast/mvn.h>
#include "nj_var.help"

#define DEFAULT_NSAMPLES 100
#define DEFAULT_DIM 3
#define DEFAULT_BATCHSIZE 20
#define DEFAULT_LEARNRATE 0.001
#define DEFAULT_NBATCHES_CONV 10
#define DEFAULT_MIN_NBATCHES 30
#define DEFAULT_KAPPA 4
#define DEFAULT_RANK 3

int main(int argc, char *argv[]) {
  signed char c;
  int opt_idx, i, ntips = 0, nsamples = DEFAULT_NSAMPLES, dim = DEFAULT_DIM,
    batchsize = DEFAULT_BATCHSIZE, nbatches_conv = DEFAULT_NBATCHES_CONV,
    min_nbatches = DEFAULT_MIN_NBATCHES, rank = DEFAULT_RANK;
  unsigned int nj_only = FALSE, random_start = FALSE,
    hyperbolic = FALSE, embedding_only = FALSE, importance_sampling = FALSE,
    mvn_dump = FALSE;
  MSA *msa = NULL;
  enum covar_type covar_param = DIAG;

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
  double kappa = DEFAULT_KAPPA, learnrate = DEFAULT_LEARNRATE, negcurvature = 1;
  MarkovMatrix *rmat = NULL;
  multi_MVN *mmvn = NULL;
  TreeNode *init_tree = NULL;
  CovarData *covar_data = NULL;

  struct option long_opts[] = {
    {"format", 1, 0, 'i'},
    {"batchsize", 1, 0, 'b'},
    {"nbatches-conv", 1, 0, 'c'},
    {"dimensionality", 1, 0, 'D'},
    {"distances", 1, 0, 'd'},
    {"embedding-only", 0, 0, 'e'},
    {"hky85", 0, 0, 'k'}, 
    {"hyperbolic", 0, 0, 'H'},
    {"importance-sampling", 0, 0, 'I'},
    {"logfile", 1, 0, 'l'},
    {"mean", 1, 0, 'm'},
    {"min-nbatches", 1, 0, 'M'},
    {"names", 1, 0, 'n'},
    {"negcurvature", 1, 0, 'K'},
    {"nj-only", 0, 0, 'j'},
    {"out-dists", 1, 0, 'o'},
    {"nsamples", 1, 0, 's'},
    {"learnrate", 1, 0, 'r'},
    {"random-start", 0, 0, 'R'},
    {"distance-covar", 1, 0, 'S'},
    {"tree", 1, 0, 't'},
    {"treemodel", 1, 0, 'T'},
    {"mvn-dump", 0, 0, 'V'},
    {"rank", 1, 0, 'W'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "b:c:d:D:ehHIi:jkK:l:m:M:n:o:r:Rt:T:VW:S:s:", long_opts, &opt_idx)) != -1) {
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
    case 'e':
      embedding_only = TRUE;
      break;
    case 'H':
      hyperbolic = TRUE;
      break;
    case 'I':
      importance_sampling = TRUE;
      break;
    case 'j':
      nj_only = TRUE;
      break;
    case 'k':
      subst_mod = HKY85;
      break;
    case 'K':
      negcurvature = atof(optarg);
      if (negcurvature < 0)
        die("ERROR: --negcurvature must be nonnegative\n");
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
    case 'S':
      if (!strcmp(optarg, "DIAG"))
        covar_param = DIAG;
      else if (!strcmp(optarg, "CONST"))
        covar_param = CONST;
      else if (!strcmp(optarg, "DIST"))
        covar_param = DIST;
      else if (!strcmp(optarg, "LOWR"))
        covar_param = LOWR;
      else die("ERROR: bad argument to --covar (-S).\n");
      break;
    case 't':
      init_tree = tr_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'T':
      mod = tm_new_from_file(phast_fopen(optarg, "r"), 1);
      init_tree = mod->tree;
      break;
    case 'V':
      mvn_dump = TRUE;
      break;
    case 'W':
      rank = atoi(optarg);
      if (rank <= 0)
        die("ERROR: --rank must be positive\n");
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

  if (hyperbolic == TRUE && negcurvature == 0) 
    hyperbolic = FALSE;
  /* for convenience in scripting; nonhyperbolic considered special case of hyperbolic */

  if (rank != DEFAULT_RANK && covar_parm != LOWR)
    fprintf("WARNING: --rank ignored when --covar is not LOWR\n");
  
  if ((nj_only || embedding_only) &&
      (indistfile != NULL || init_tree != NULL)) {
    if (optind != argc) 
      die("ERROR: No alignment needed in this case.  Too many arguments.  Try 'nj_var -h'.\n");
  }
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

    names = msa->names;
    ntips = msa->nseqs;
  }
  
  if (msa == NULL && names == NULL) {
    if (init_tree) {
      List *namelst = tr_leaf_names(init_tree); /* have to convert to char arrays */
      ntips = lst_size(namelst);
      names = smalloc(sizeof(char*)*ntips);
      for (i = 0; i < ntips; i++) {
        String *str = lst_get_ptr(namelst, i);
        names[i] = smalloc(sizeof(char) * (str->length+1));
        strcpy(names[i], str->chars);
      }
    }
    else 
      die("ERROR: must specify alignment, --tree/--treemod, or --names.\n");
  }

  /* at this point, names and ntips must be defined even if we don't have an alignment */

    
  /* get a distance matrix */
  if (init_tree != NULL)
    D = nj_tree_to_distances(init_tree, names, ntips);  
  else if (indistfile != NULL) 
    D = mat_new_from_file(indistfile, ntips, ntips);
  else if (msa != NULL)
    D = nj_compute_JC_matr(msa);
  else
    die("ERROR: no distance matrix available\n");
  
  /* we must have a distance matrix now; make sure valid */
  nj_test_D(D);

  covar_data = nj_new_covar_data(covar_param, D, dim, rank);
  
  if (embedding_only == TRUE) {
    /* in this case, embed the distances now */
    if (outdistfile == NULL)
      die("ERROR: must use --out-dists with -embedding-only\n");

    mmvn = mmvn_new(ntips, dim, covar_data->mvn_type);
    nj_estimate_mmvn_from_distances(D, dim, mmvn, negcurvature,
                                    covar_data, hyperbolic);
  }

  else {
    /* we'll need a starting NJ tree for either variational inference
       or NJ */
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
      mmvn = mmvn_new(ntips, dim, covar_data->mvn_type);
      if (random_start == TRUE) 
        mvn_sample_std(mmvn->mvn->mu); vec_scale(mmvn->mvn->mu, 0.1);
      else 
        nj_estimate_mmvn_from_distances(D, dim, mmvn, negcurvature,
                                        covar_data, hyperbolic);

      if (mvn_dump) {  /* in this case, just dump the MVN and associated data for inspection */
        mmvn_print(mmvn, stdout, FALSE, TRUE);
        nj_dump_covar_data(covar_data, stdout);
        exit(0);
      }
      
      nj_variational_inf(mod, msa, D, mmvn, dim, hyperbolic,
                         negcurvature, batchsize, learnrate,
                         nbatches_conv, min_nbatches, 
                         covar_data, logfile);

      if (importance_sampling == TRUE) {
        /* sample 100x as many then importance sample; make free param? */
        Vector *logdens = vec_new(100*nsamples);
        List *origtrees = nj_var_sample(100*nsamples, dim, mmvn,
                                        msa->names, hyperbolic,
                                        negcurvature, logdens);
        trees = nj_importance_sample(nsamples, origtrees, logdens, mod, msa, logfile);        
      }

      else /* otherwise just sample directly from posterior */
        trees = nj_var_sample(nsamples, dim, mmvn, msa->names,
                              hyperbolic, negcurvature, NULL);

      for (i = 0; i < nsamples; i++)
        tr_print(stdout, (TreeNode*)lst_get_ptr(trees, i), TRUE);

      if (postmeanfile != NULL) {
        Vector *mu_full = vec_new(mmvn->d * mmvn->n);
        mmvn_save_mu(mmvn, mu_full);
        tr_print(postmeanfile, nj_mean(mu_full, dim, msa->names,
                                       hyperbolic, negcurvature), TRUE);
        vec_free(mu_full);
      }
    }
  }

  if (outdistfile != NULL) {
    if (embedding_only == TRUE || nj_only == FALSE)
      /* in this case need to reset D */
      nj_mmvn_to_distances(mmvn, D, hyperbolic, negcurvature);

    mat_print(D, outdistfile);
  }
    
  return (0);
}
