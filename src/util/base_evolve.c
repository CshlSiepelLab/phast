#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <msa.h>
#include <category_map.h>
#include <tree_model.h>
#include <time.h>
#include "base_evolve.help"

int main(int argc, char *argv[]) {
  /* variables for options with default */
  int nsites = 1000, embed_len = -1;
  CategoryMap *cm = NULL;
  char *features_fname = NULL;
  msa_format_type msa_format = FASTA;
  TreeModel *embed_mod = NULL;
  
  /* other variables */
  FILE *F;
  TreeModel **mods;
  HMM *hmm;
  MSA *msa;
  int *labels = NULL, *path_to_cat, *reverse_compl;
  GFF_Set *feats;
  char c;
  int opt_idx, i, j;
  List *l;
  struct timeval now;

  struct option long_opts[] = {
    {"nsites", 1, 0, 'n'},
    {"msa-format", 1, 0, 'o'},
    {"features", 1, 0, 'f'},
    {"catmap", 1, 0, 'c'},
    {"embed", 1, 0, 'e'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "n:o:f:c:e:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'n':
      nsites = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'o':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == -1)
        die("ERROR: unrecognized alignment format.\n");
      break;
    case 'f':
      features_fname = optarg;
      break;
    case 'c':
      cm = cm_new_string_or_file(optarg);
      break;
    case 'e':
      l = get_arg_list(optarg);
      embed_mod = tm_new_from_file(fopen_fname(((String*)lst_get_ptr(l, 0))->chars, "r"));
      embed_len = get_arg_dbl_bounds(((String*)lst_get_ptr(l, 1))->chars, 1, INFTY);
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try base_evolve -h'.\n");
    }
  }

  if (optind < argc - 1)
    die("ERROR: base_evolve requires one or more arguments.  Try base_evolve -h.\n");

  if (optind == argc - 1)       /* single tree model */
    hmm = hmm_create_trivial();
  else {
    hmm = hmm_new_from_file(fopen_fname(argv[optind], "r"));
    optind++;
  }

  if (optind != argc - hmm->nstates)
    die("ERROR: number of tree models must equal number of states in HMM.\n");

  mods = smalloc(hmm->nstates * sizeof(void*));
  for (i = 0; i < hmm->nstates; i++) {
    F = fopen_fname(argv[optind+i], "r");
    mods[i] = tm_new_from_file(F);
    if (mods[i]->nratecats != 1)
      die("ERROR: rate variation currently not supported.\n");
    fclose(F);
  }  

  /* generate alignment and labels */
  if (features_fname != NULL)
    labels = smalloc(nsites * sizeof(int));
  gettimeofday(&now, NULL);
  srandom(now.tv_usec);         /* use microseconds to avoid using
                                   same seed in rapidly repeated
                                   calls */
  msa = tm_generate_msa(nsites, hmm, mods, labels);

  /* generate features, if necessary */
  if (features_fname != NULL) {
    if (cm == NULL) 
      cm = cm_create_trivial(hmm->nstates, "model_");

    path_to_cat = smalloc(hmm->nstates * sizeof(int));
    reverse_compl = smalloc(hmm->nstates * sizeof(int));
    for (i = 0; i < hmm->nstates; i++) {
      path_to_cat[i] = i; 
      reverse_compl[i] = FALSE;
    }
    feats = cm_labeling_as_gff(cm, labels, msa->length, path_to_cat, 
                             reverse_compl, "sim", "base_evolve", 
                             NULL, NULL, NULL);
    free(path_to_cat);
    free(reverse_compl);

    F = fopen_fname(features_fname, "w+");
    gff_print_set(F, feats);
    fclose(F);
  }

  /* add embedded element, if necessary */
  if (embed_mod != NULL) {
    MSA *embed_msa = tm_generate_msa(embed_len, NULL, &embed_mod, NULL);
    int startidx = (msa->length - embed_len)/2 + 1; 
    for (i = 0; i < embed_msa->length; i++)
      for (j = 0; j < msa->nseqs; j++)
        msa->seqs[j][startidx+i] = embed_msa->seqs[j][i];
  }

  /* print alignment */
  msa_print(stdout, msa, msa_format, 0);

  return 0;
}
