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
#include <tree_likelihoods.h>
#include "consEntropy.help"

/* solve for new expected length given L_min*H using Newton's method */
double solve_newton(double expected_len, double target_coverage, double H, double LminH) {
  double L_min, odds, mu1, mu2, func, deriv;
  int niter = 0;
  H = H * log(2); /* switch to natural log scale -- makes derivatives simpler */
  LminH = LminH * log(2);
  L_min = LminH/H;
  odds = target_coverage / (1-target_coverage);
  mu1 = 1/expected_len;
  fprintf(stderr, "\n( Solving for new omega: %f ", 1/mu1);
  while (TRUE) {
    func = (L_min+1) * log(1-odds*mu1) - (L_min-1) * log(1-mu1) - log(odds*mu1) - log(mu1) - LminH;
    deriv = -(L_min+1) * odds / (1 - odds * mu1) + (L_min-1) / (1-mu1) - 2/mu1;
    mu2 = mu1 - func/deriv;
    if (mu2 < 0) mu2 = 1e-3;
    else if (mu2 > 1) mu2 = 1 - 1e-3;
    fprintf(stderr, "%f ", 1/mu2);
    if (fabs(mu2-mu1) < 1e-4) break;
    mu1 = mu2;
    niter++;
    if (niter >= 30) die("ERROR: too many iterations, not converging; try without --NH.");
  }
  fprintf(stderr, ")\n\n");
  return 1/mu2;
}

int main(int argc, char *argv[]) {
  char c;
  int i, j, opt_idx, nleaves, alph_size, nlabels;
  TreeModel *cons_mod, *noncons_mod;
  MSA *msa;
  double *cons_lprob, *noncons_lprob;
  char *leaf_labels;
  double H = -1, H_alt = -1, checksum1, checksum2, target_coverage, 
    expected_len, mu, nu, L_min, L_max, new_exp_len, LminH = -1;

  struct option long_opts[] = {
    {"LminH", 1, 0, 'L'},
    {"NH", 1, 0, 'N'},		/* backward compatibility */
    {"H", 1, 0, 'H'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "H:N::h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'H':
      H = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 'N':			/* backward compatibility */
    case 'L':
      LminH = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if ((H == -1 && optind != argc - 4) || (H != -1 && optind != argc - 2))
    die("Missing mandatory arguments.  Try '%s -h'.\n", argv[0]);

  set_seed(-1);
    
  target_coverage = get_arg_dbl_bounds(argv[optind], 0, 1);
  expected_len = get_arg_dbl_bounds(argv[optind+1], 0, INFTY);

  if (H == -1) {                /* compute relative entropy */
    cons_mod = tm_new_from_file(phast_fopen(argv[optind+2], "r"), 1);
    noncons_mod = tm_new_from_file(phast_fopen(argv[optind+3], "r"), 1);

    nleaves = (cons_mod->tree->nnodes + 1)/2;
    leaf_labels = smalloc((nleaves + 1) * sizeof(char));
    leaf_labels[nleaves] = '\0';
    alph_size = strlen(cons_mod->rate_matrix->states);
    nlabels = int_pow(alph_size, nleaves);

    /* define dummy MSA */
    msa = msa_new(NULL, NULL, nleaves, nlabels, cons_mod->rate_matrix->states);
    msa->seqs = smalloc(nleaves * sizeof(void*));
    msa->names = smalloc(nleaves * sizeof(void*));
    for (j = 0; j < nleaves; j++) 
      msa->seqs[j] = smalloc((nlabels+1) * sizeof(char));
    for (j = 0, i = 0; j < cons_mod->tree->nnodes; j++) {
      TreeNode *n = lst_get_ptr(cons_mod->tree->nodes, j);
      if (n->lchild == NULL && n->rchild == NULL)
        msa->names[i++] = copy_charstr(n->name);
    }

    /* enumerate all possible columns and put in MSA */
    for (i = 0; i < nlabels; i++) {
      get_tuple_str(leaf_labels, i, nleaves, cons_mod->rate_matrix->states);
      for (j = 0; j < nleaves; j++) msa->seqs[j][i] = leaf_labels[j];
    }

    /* compute log likelihoods of all columns */
    cons_lprob = smalloc(nlabels * sizeof(double));
    noncons_lprob = smalloc(nlabels * sizeof(double));
    tl_compute_log_likelihood(cons_mod, msa, cons_lprob, NULL, -1, NULL);
    tl_compute_log_likelihood(noncons_mod, msa, noncons_lprob, NULL, -1, NULL);

    H = H_alt = 0;		/* H is relative entropy of cons wrt
				   noncons; H_alt is relative entropy
				   of noncons wrt cons */
    checksum1 = checksum2 = 0;
    for (i = 0; i < nlabels; i++) {
      double tmp = exp2(cons_lprob[i]); /* tl_compute_log_likelihood uses base 2 */
      double tmp2 = exp2(noncons_lprob[i]);
      checksum1 += tmp;
      checksum2 += tmp2;
      H += tmp * (cons_lprob[i] - noncons_lprob[i]);
      H_alt += tmp2 * (noncons_lprob[i] - cons_lprob[i]);
    }

    if (fabs(checksum1 - 1) > 1e-4 || fabs(checksum1 - 1) > 1e-4)
      die("ERROR: checksum failed (%f or %f not 1 +/- 1.0e-4).\n", checksum1, checksum2);
  }

  mu = 1/expected_len;
  nu = mu * target_coverage / (1-target_coverage);
  L_min = (log2(nu) + log2(mu) - log2(1-nu) - log2(1-mu)) / (log2(1-nu) - log2(1-mu) - H);
  L_max = (log2(nu) + log2(mu) - log2(1-nu) - log2(1-mu)) / (log2(1-mu) - log2(1-nu) - H_alt);

  if (LminH > -1) 
    new_exp_len = solve_newton(expected_len, target_coverage, H, LminH);

  printf("Transition parameters: gamma=%f, omega=%f, mu=%f, nu=%f\n", 
         target_coverage, expected_len, mu, nu);
  printf("Relative entropy: H=%f bits/site\n", H);
  printf("Expected min. length: L_min=%f sites\n", L_min);
  printf("Expected max. length: L_max=%f sites\n", L_max);
  printf("Phylogenetic information threshold: PIT=L_min*H=%f bits\n", L_min*H);
  if (LminH > -1)
    printf("Recommended expected length: omega=%f sites (for L_min*H=%f)\n", new_exp_len, LminH);
  return 0;
}
