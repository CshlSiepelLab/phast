#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <tree_model.h>
#include <msa.h>
#include <tree_likelihoods.h>

void usage(char *prog) {
  printf("\n\
PROGRAM: %s\n\
\n\
DESCRIPTION:\n\
    For use with phastCons.  Given phylogenetic models for conserved\n\
    and non-conserved states, the target coverage, and the (prior)\n\
    expected length of a conserved element, compute the relative\n\
    entropy (H) of the phylogenetic models and the expected number (N)\n\
    of conserved sites required to predict conserved element.\n\
    Also will make a recommendation for a new prior expected length\n\
    based on a constant value of NH (see --NH).\n\
\n\
USAGE: %s [OPTIONS] <cons.mod> <noncons.mod> <target-coverage> \\\n\
                <expected-length>\n\
\n\
OPTIONS:\n\
    --NH, -n <value>\n\
        Report the expected length that would produce the specified\n\
        value of NH, assuming H remains constant (it generally won't).\n\
        Can be used iteratively to converge on a desired value of NH.\n\
\n\
    --help, -h\n\
        Print this help message.\n\
\n\
NOTE:\n\
    The relative entropy is currently computed by brute force, i.e.,\n\
    by enumerating all possible labelings of the leaves of the tree.\n\
    This approach won't be feasible with large trees.\n\n", prog, prog);
  exit(0);
}

/* solve for new expected length given NH using Newton's method */
double solve_newton(double expected_len, double target_coverage, double H, double NH) {
  double N, odds, mu1, mu2, func, deriv;
  H = H * log(2);               /* switch to natural log scale -- makes derivatives simpler */
  NH = NH * log(2);
  N = NH/H;
  odds = target_coverage / (1-target_coverage);
  mu1 = 1/expected_len;
  fprintf(stderr, "\n( Solving for new omega: %f ", 1/mu1);
  while (TRUE) {
    func = (N+1) * log(1-odds*mu1) - (N-1) * log(1-mu1) - log(odds*mu1) - log(mu1) - NH;
    deriv = -(N+1) * odds / (1 - odds * mu1) + (N-1) / (1-mu1) - 2/mu1;
    mu2 = mu1 - func/deriv;
    if (mu2 < 0) mu2 = 1e-3;
    else if (mu2 > 1) mu2 = 1 - 1e-3;
    fprintf(stderr, "%f ", 1/mu2);
    if (fabs(mu2-mu1) < 1e-4) break;
    mu1 = mu2;
  }
  fprintf(stderr, ")\n\n");
  return 1/mu2;
}

int main(int argc, char *argv[]) {
  char c;
  int i, j, opt_idx, nleaves, alph_size, nlabels;
  double H, checksum1, checksum2;
  TreeModel *cons_mod, *noncons_mod;
  MSA *msa;
  double *cons_lprob, *noncons_lprob;
  char *leaf_labels;
  double target_coverage, expected_len, mu, nu, N, new_exp_len, NH = -1;

  struct option long_opts[] = {
    {"NH", 1, 0, 'n'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "n:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'n':
      NH = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 'h':
      usage(argv[0]);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if (optind != argc - 4) 
    die("Four arguments required.  Try '%s -h'.\n", argv[0]);
    
  cons_mod = tm_new_from_file(fopen_fname(argv[optind], "r"));
  noncons_mod = tm_new_from_file(fopen_fname(argv[optind+1], "r"));
  target_coverage = get_arg_dbl_bounds(argv[optind+2], 0, 1);
  expected_len = get_arg_dbl_bounds(argv[optind+3], 0, INFTY);

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
      msa->names[i++] = strdup(n->name);
  }

  /* enumerate all possible columns and put in MSA */
  for (i = 0; i < nlabels; i++) {
    get_tuple_str(leaf_labels, i, nleaves, cons_mod->rate_matrix->states);
    for (j = 0; j < nleaves; j++) msa->seqs[j][i] = leaf_labels[j];
  }

  /* compute log likelihoods of all columns */
  cons_lprob = smalloc(nlabels * sizeof(double));
  noncons_lprob = smalloc(nlabels * sizeof(double));
  tl_compute_log_likelihood(cons_mod, msa, cons_lprob, -1, NULL);
  tl_compute_log_likelihood(noncons_mod, msa, noncons_lprob, -1, NULL);

  H = 0;
  checksum1 = checksum2 = 0;
  for (i = 0; i < nlabels; i++) {
    double tmp = exp2(cons_lprob[i]); /* tl_compute_log_likelihood uses base 2 */
    checksum1 += tmp;
    checksum2 += exp2(noncons_lprob[i]);
    H += tmp * (cons_lprob[i] - noncons_lprob[i]);
  }

  if (fabs(checksum1 - 1) > 1e-4 || fabs(checksum1 - 1) > 1e-4)
    die("ERROR: checksum failed (%f or %f not 1 +/- 1.0e-4).\n", checksum1, checksum2);

  mu = 1/expected_len;
  nu = mu * target_coverage / (1-target_coverage);
  N = (log2(nu) + log2(mu) - log2(1-nu) - log2(1-mu)) / (log2(1-nu) - log2(1-mu) - H);

  if (NH > -1) 
    new_exp_len = solve_newton(expected_len, target_coverage, H, NH);

  printf("Transition parameters: gamma=%f, omega=%f, mu=%f, nu=%f\n", 
         target_coverage, expected_len, mu, nu);
  printf("Relative entropy: H=%f bits/site\n", H);
  printf("Required length: N=%f sites\n", N);
  printf("Total entropy: NH=%f bits\n", N*H);
  if (NH > -1)
    printf("Recommended expected length: omega=%f sites (for NH=%f)\n", new_exp_len, NH);
  return 0;
}
