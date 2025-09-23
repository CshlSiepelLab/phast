#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <phast/misc.h>
#include <phast/nj.h>
#include <phast/sufficient_stats.h>
#include <phast/sparse_matrix.h>
#include <phast/lists.h>
#include <phast/msa.h>
#include <phast/sufficient_stats.h>
#include <phast/subst_mods.h>
#include <phast/tree_model.h>
#include "evalTrees.help"

int main(int argc, char *argv[]) {
  TreeNode *tree;
  TreeModel *mod = NULL;
  double kappa = -1, ll, lltot;
  String *line = str_new(STR_VERY_LONG_LEN);
  int opt_idx, lineno = 0;
  CovarData *data;
  Matrix *D;
  char c;
  FILE *treefile, *msafile;
  MSA *msa;
  MarkovMatrix *rmat;
  msa_format_type format;
  
  struct option long_opts[] = {
    {"hky-kappa", 1, 0, 'k'},
    {"tree-model", 1, 0, 'm'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "k:m:h", 
                          long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'k':
      kappa = atof(optarg);
      if (kappa < 0)
        die("ERROR: --hky-kappa must be > 0.\n");
      break;
    case 'm':
      mod = tm_new_from_file(phast_fopen(optarg, "r"), 1);
      break;
    case 'h':
      printf("%s", HELP); 
      exit(0);
    case '?':
      die("Bad argument.  Try 'evalTrees -h'.\n");
    }
  }

  if (optind != argc - 2)
    die("Missing required argument.  Try '%s -h'.\n", argv[0]);

  if (mod != NULL && kappa > 0)
    die("Select either --tree-model or --hky-kappa but not both.\n");

  /* open input files, read alignment */
  treefile = phast_fopen(argv[optind++], "r");
  msafile = phast_fopen(argv[optind], "r");
  format = msa_format_for_content(msafile, 1);
  msa = msa_new_from_file_define_format(msafile, format, DEFAULT_ALPHABET);

  /* set up alignment */
  if (msa->ss == NULL)
    ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1, 0);

  /* this is mostly a dummy; only the alignment is used */
  D = mat_new(5, 5);
  data = nj_new_covar_data(CONST, D, 1, msa, NULL, NULL, FALSE, 1.0, 3,
                           1.0, FALSE, -1, FALSE, FALSE, FALSE);

  lltot = 0;
  while (str_readline(line, treefile) != EOF) {
    lineno++;
    str_double_trim(line);

    if (line->chars[0] != '(')
      die("ERROR in line %d: Input does not look like a Newick-formatted tree.\n",
          lineno);
    if (line->chars[line->length-1] == ';')
      line->chars[--line->length] = '\0';

    tree = tr_new_from_string(line->chars);

    if (mod == NULL) { /* do this the first time through; need a tree to initialize */
      rmat = mm_new(strlen(DEFAULT_ALPHABET), DEFAULT_ALPHABET, CONTINUOUS);
      mod = tm_new(tree, rmat, NULL, kappa > 0 ? HKY85 : JC69, DEFAULT_ALPHABET,
                   1, 1, NULL, -1);
      tm_init_backgd(mod, msa, -1);
      if (kappa > 0) /* create HKY model */ {
        fprintf(stderr, "Using HKY85 with kappa = %f...\n", kappa);
        tm_set_HKY_matrix(mod, kappa, -1);
      }
      else {           /* create JC model */
        fprintf(stderr, "Using JC69...\n");
        tm_set_JC69_matrix(mod);
      }
    }
    else
      nj_reset_tree_model(mod, tree);

    /* have to force index rebuild because node ids can change */
    sfree(mod->msa_seq_idx);
    tm_build_seq_idx(mod, msa);
    
    ll = nj_compute_log_likelihood(mod, data, NULL);

    /* temporary */
    printf("Tree %d: ll = %f\n", lineno, ll);

    lltot += ll;    
  }

  printf("Successfully processed %d trees.\nAverage loglikelihood: %f\nPer site: %f\n",
         lineno, lltot / lineno, lltot / (lineno * msa->length));
  fprintf(stderr, "Done processing %d trees.\n", lineno);
}


  
