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
#include <phast/rf.h>
#include "evalTrees.help"

static inline
void print_stats(FILE *F, double mean, double stdev, double median,
                 double min, double max, double min_95CI, double max_95CI,
                 double q25, double q75) {
  fprintf(F, "Mean: %f\n", mean);
  fprintf(F, "Std: %f\n", stdev);
  fprintf(F, "Median: %f\n", median);
  fprintf(F, "Range: %f - %f\n", min, max);
  fprintf(F, "95%%_CI: %f - %f\n", min_95CI, max_95CI);
  fprintf(F, "50%%_CI: %f - %f\n", q25, q75);
}

int main(int argc, char *argv[]) {
  TreeNode *tree;
  TreeModel *mod = NULL;
  double kappa = -1, ll;
  String *line = str_new(STR_VERY_LONG_LEN);
  int opt_idx, lineno = 0, i, j, nleaves = 0, npairs = 0;
  CovarData *data;
  char c;
  FILE *treefile, *msafile = NULL;
  MarkovMatrix *rmat;
  msa_format_type format;
  TreeNode *topol_ref = NULL;
  MSA *evalaln = NULL;
  double mean, stdev, median, min, max, min_95CI, max_95CI, q25, q75;
  char *topolfname, *msafname, *treefname;
  List *rfdists, *lldists;
  char **names = NULL;
  Matrix *D = NULL;
  List **Dij_list = NULL;
  int is_crispr = FALSE;
  CrisprMutTable *crispr_muts = NULL;
  CrisprMutModel *crispr_mod = NULL;
  enum crispr_model_type crispr_modtype = SITEWISE;
  enum crispr_mutrates_type crispr_muttype = UNIF;
  
  struct option long_opts[] = {
    {"hky-kappa", 1, 0, 'k'},
    {"crispr", 0, 0, 'k'},
    {"tree-model", 1, 0, 'm'},
    {"model-fit", 1, 0, 'f'},
    {"topology", 1, 0, 't'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "f:k:m:t:ch", 
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
    case 'f':
      msafname = optarg;
      msafile = phast_fopen(optarg, "r");
      break;
    case 't':
      topolfname = optarg;
      topol_ref = tr_new_from_file(phast_fopen(topolfname, "r"));
      break;
    case 'c':
      is_crispr = TRUE;
      break;
    case 'h':
      printf("%s", HELP); 
      exit(0);
    case '?':
      die("Bad argument.  Try 'evalTrees -h'.\n");
    }
  }

  if (optind != argc - 1)
    die("Missing required argument.  Try '%s -h'.\n", argv[0]);

    /* open dna or crispr file */
  if (is_crispr == TRUE) {
    if (msafile == NULL)
      die("Option --crispr requires --model-fit.\n");
    crispr_muts = cpr_read_table(msafile);
  }
  else if (msafile != NULL) {
    format = msa_format_for_content(msafile, 1);
    evalaln = msa_new_from_file_define_format(msafile, format, DEFAULT_ALPHABET);
  }
  
  if (evalaln == NULL && (mod != NULL || kappa > 0)) 
    die("Options --tree-model and --hky-kappa require --model-fit.\n");
  
  if (evalaln != NULL && mod != NULL && kappa > 0)
    die("Options --tree-model and --hky-kappa are mutually exclusive.\n");

  if (evalaln != NULL && topol_ref != NULL)
    die("Options --model-fit and --topology are mutually exclusive.\n");

  if (is_crispr && (mod != NULL || kappa > 0))
    die("Options --tree-model and --hky-kappa are incompatible with --crispr.\n");
  
  /* open tree file */
  treefname = argv[optind];
  fprintf(stderr, "Reading trees from %s...\n", treefname);
  treefile = phast_fopen(treefname, "r");  
  
  /* set up for --model-fit */
  if (evalaln != NULL || is_crispr == TRUE) {
    fprintf(stderr, "Evaluating model fit on %s...\n", msafname);

    if (evalaln != NULL && evalaln->ss == NULL)
      ss_from_msas(evalaln, 1, TRUE, NULL, NULL, NULL, -1, 0);
    else if (is_crispr)
      crispr_mod = cpr_new_model(crispr_muts, NULL, crispr_modtype, crispr_muttype);
      
    /* this is mostly a dummy; only the msa or crispr mod is used */
    D = mat_new(5, 5);
    data = nj_new_covar_data(CONST, D, 1, evalaln, crispr_mod, NULL, FALSE,
                             1.0, 3, 1.0, FALSE, -1, FALSE, FALSE, FALSE);
    lldists = lst_new_dbl(1000);
  }
  else if (topol_ref != NULL) {
    rfdists = lst_new_dbl(1000);
    fprintf(stderr, "Evaluating RF distance to %s...\n", topolfname);
  }
  else
    fprintf(stderr, "Computing pairwise-distance stats...\n");
  
  while (str_readline(line, treefile) != EOF) {
    str_double_trim(line);

    if (line->length == 0)
      continue;

    lineno++;

    if (line->chars[0] != '(')
      die("ERROR in line %d: Input does not look like a Newick-formatted tree.\n",
          lineno);
    if (line->chars[line->length-1] == ';')
      line->chars[--line->length] = '\0';

    tree = tr_new_from_string(line->chars);

    if (evalaln != NULL) {
      if (mod == NULL) { /* do this the first time through; need a tree to initialize */        
        rmat = mm_new(strlen(DEFAULT_ALPHABET), DEFAULT_ALPHABET, CONTINUOUS);
        mod = tm_new(tree, rmat, NULL, kappa > 0 ? HKY85 : JC69, DEFAULT_ALPHABET,
                     1, 1, NULL, -1);
        if (is_crispr) {  /* tree model just a dummy in this case */
          crispr_mod->mod = mod;
          cpr_prep_model(crispr_mod);
        }
        else {
          tm_init_backgd(mod, evalaln, -1);
          if (kappa > 0) /* create HKY model */ {
            fprintf(stderr, "Using HKY85 with kappa = %f...\n", kappa);
            tm_set_HKY_matrix(mod, kappa, -1);
          }
          else {           /* create JC model */
            fprintf(stderr, "Using JC69...\n");
            tm_set_JC69_matrix(mod);
          }
        }
      }
      else
        nj_reset_tree_model(mod, tree);

      if (evalaln != NULL) {
        /* have to force index rebuild because node ids can change */
        sfree(mod->msa_seq_idx);
        tm_build_seq_idx(mod, evalaln);
        ll = nj_compute_log_likelihood(mod, data, NULL);
      }
      else { /* crispr case */
        sfree(crispr_mod->mod->msa_seq_idx);
        cpr_build_seq_idx(crispr_mod->mod, crispr_mod->mut);
        ll = cpr_compute_log_likelihood(data->crispr_mod, NULL);
      }

      /* occasionally get -inf from 0-length branches; let's just
         ignore those for now */
      if (isfinite(ll))  
        lst_push_dbl(lldists, ll);
    }

    else if (topol_ref != NULL) {
      double d = tr_robinson_foulds(tree, topol_ref);
      lst_push_dbl(rfdists, d);
    }

    else {  /* collect distance statistics */
      if (names == NULL) {  /* first time through get canonical list
                               of leaf names */
        List *l = tr_leaf_names(tree);
        lst_qsort_str(l, ASCENDING);
        nleaves = lst_size(l);
        npairs = nleaves * (nleaves-1) / 2;
        names = smalloc(nleaves * sizeof(char*)); 
        for (i = 0; i < lst_size(l); i++) {
          String *s = lst_get_ptr(l, i);
          names[i] = s->chars;
        }
        /* also set up lists of pairwise distances across trees */
        Dij_list = smalloc(npairs * sizeof(List*));
        for (i = 0; i < nleaves; i++) 
          for (j = i+1; j < nleaves; j++)
            Dij_list[nj_i_j_to_dist(i, j, nleaves)] = lst_new_dbl(1000);
        fprintf(stderr, "Extracting pairwise distances for all trees...\n");
      }
      
      /* get distance matrix for all pairs of leaves for this tree */
      D = nj_tree_to_distances(tree, names, nleaves);
      /* add distances to corresponding lists */
      for (i = 0; i < nleaves; i++) {
        for (j = i+1; j < nleaves; j++) {
          double d = mat_get(D, i, j);
          lst_push_dbl(Dij_list[nj_i_j_to_dist(i, j, nleaves)], d);
        }
      }
      
      mat_free(D);
    }
  }

  /* output results */
  fprintf(stderr, "Done processing %d trees.\n", lineno);
   
  if (evalaln != NULL) {
    printf("Successfully processed %d trees from %s.\n", lineno, treefname);
    printf("Log likelihood evaluated on %s:\n", msafname);
    lst_dbl_stats(lldists, &mean, &stdev, &median, &min, &max,
                  &min_95CI, &max_95CI, &q25, &q75);
    print_stats(stdout, mean, stdev, median, min, max, min_95CI,
                max_95CI, q25, q75);
    printf("Mean per site: %f\n", mean/evalaln->length);
  }
  else if (topol_ref != NULL) {
    printf("Successfully processed %d trees from %s.\n", lineno, treefname);
    lst_dbl_stats(rfdists, &mean, &stdev, &median, &min, &max,
                  &min_95CI, &max_95CI, &q25, &q75);
    printf("Robinson Foulds distances against %s:\n", topolfname);
    print_stats(stdout, mean, stdev, median, min, max, min_95CI,
                max_95CI, q25, q75);
  }
  else {
    printf("#leaf1\tleaf2\tmean\tstd\tmed\tmin\tmax\tlow95CI\thigh95CI\tlow50CI\thigh50CI\n");
    for (i = 0; i < nleaves; i++) {
      for (j = i+1; j < nleaves; j++) {
        lst_dbl_stats(Dij_list[nj_i_j_to_dist(i, j, nleaves)], &mean, &stdev,
                      &median, &min, &max, &min_95CI, &max_95CI, &q25, &q75);
        /* in this case, print a table */
        printf("%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
               names[i], names[j], mean, stdev, median, min, max, min_95CI,
               max_95CI, q25, q75);
      }
    }
  }
}


  
