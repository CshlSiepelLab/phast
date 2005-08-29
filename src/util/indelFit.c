#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <tree_model.h>
#include <indel_history.h>
#include <indel_mod.h>
#include "indelFit.help"

int main(int argc, char *argv[]) {
  char c;
  int opt_idx;
  TreeNode *tree;
  IndelHistory *ih;
  IndelModel *im;
  double lnl = INFTY;

  /* variables for options with defaults */
  double alpha = 0.05, beta = 0.05, tau = 0.2;
  int lnl_only = FALSE, columns = FALSE;
  FILE *logf = NULL;

  struct option long_opts[] = {
    {"alpha", 1, 0, 'a'},
    {"beta", 1, 0, 'b'},
    {"tau", 1, 0, 't'},
    {"lnl", 0, 0, 'L'},
    {"columns", 0, 0, 'c'},
    {"log", 1, 0, 'l'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "a:b:o:Lcl:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'a':
      alpha = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'b':
      beta = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'o':
      tau = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'L':
      lnl_only = TRUE;
      break;
    case 'c':
      columns = TRUE;
      break;
    case 'l':
      logf = fopen_fname(optarg, "w+");
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'indelFit -h'.\n");
    }
  }

  if (optind != argc - 2) 
    die("ERROR: Two arguments required.  Try 'indelFit -h'.\n");
  else if (lnl_only && logf != NULL)
    die("WARNING: --log ignored.\n");

  ih = ih_new_from_file(fopen_fname(argv[optind], "r"));
  tree = tr_new_from_file(fopen_fname(argv[optind+1], "r"));

  im = im_new(alpha, beta, tau, tree);

  if (lnl_only) {
    IndelSuffStats *ss = im_suff_stats(ih);
    lnl = im_likelihood(im, ss);
  }
  else {
    im_estimate(im, ih, logf);
    lnl = im->training_lnl;
  }

  if (columns) {
    int i;
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
    printf("alpha = %f, beta = %f, tau = %f\n", im->alpha, im->beta, 
           im->tau);
    printf("total lnl = %f\n", lnl);
  }

  return 0;
}
