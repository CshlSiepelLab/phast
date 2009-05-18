#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <getopt.h>
#include <misc.h>
#include <maf.h>
#include <sufficient_stats.h>
#include <tree_likelihoods.h>
#include <phylo_hmm.h>
#include <indel_history.h>
#include <dmotif_indel_mod.h>
#include <subst_distrib.h>
#include <dmotif_phmm.h>
#include <pssm.h>
#include "dmsimulate.help"

#define DEFAULT_RHO 0.3
#define DEFAULT_PHI 0.5
#define DEFAULT_MU 0.01
#define DEFAULT_NU 0.01
#define DEFAULT_ZETA 0.001
#define DEFAULT_MSA_LEN 100000

int main(int argc, char *argv[]) {
  char c, *msa_fname;
  int i, j, k, opt_idx, *path, len;
  MSA *msa, *indel_msa;
  List *tmpl, *leaf_seqs;
  DMotifPhyloHmm *dm;
  GFF_Set *motifs, *tmp_gff;
  GFF_Feature *f;
  PSSM *motif;
  IndelHistory *ih;
  TreeNode *n;
  List *zeroed_states;

  struct option long_opts[] = {
    {"refseq", 1, 0, 'M'},
    {"msa-format", 1, 0, 'i'},
    {"refidx", 1, 0, 'r'},
    {"rho", 1, 0, 'R'},
    {"phi", 1, 0, 'p'},
    {"zeta", 1, 0, 'z'},
    {"transitions", 1, 0, 't'},    
    {"expected-length", 1, 0, 'E'},
    {"target-coverage", 1, 0, 'C'},
    {"seqname", 1, 0, 'N'},
    {"idpref", 1, 0, 'P'},
    {"indel-model", 1, 0, 'I'},
    {"msa-length", 1, 0, 'L'},
    {"random-lengths", 1, 0, 'l'},
    {"keep-ancestral", 0, 0, 'k'},
    {"cond-on-subs", 0, 0, 'X'},
    {"cond-on-species", 1, 0, 'x'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };
  
  /* arguments and defaults for options */
  FILE *refseq_f = NULL, *msa_f = NULL, *motif_f = NULL, *cond_spec_f = NULL;
  msa_format_type msa_format = FASTA;
  TreeModel *source_mod;
  double rho = DEFAULT_RHO, mu = DEFAULT_MU, nu = DEFAULT_NU, 
    phi = DEFAULT_PHI, zeta = DEFAULT_ZETA, gamma = -1, omega = -1, 
    alpha_c = -1, beta_c = -1, tau_c = -1, epsilon_c = -1,
    alpha_n = -1, beta_n = -1, tau_n = -1, epsilon_n = -1,
    lambda = -1;
  int set_transitions = FALSE, refidx = 1, msa_len = DEFAULT_MSA_LEN,
    max_len = 0, min_len = 0, do_ih = FALSE, keep_ancestral = FALSE, 
    do_zeroed = FALSE;
  char *seqname = NULL, *idpref = NULL, *seqname_root = NULL;
  
  while ((c = getopt_long(argc, argv, "R:t:p:z:E:C:r:M:i:N:P:I:L:l:k:h", 
			  long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'R':
      rho = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 't':
      set_transitions = TRUE;
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 2) 
        die("ERROR: bad argument to --transitions.\n");
      mu = lst_get_dbl(tmpl, 0);
      nu = lst_get_dbl(tmpl, 1);
      if (mu <= 0 || mu >= 1 || nu <= 0 || nu >= 1)
        die("ERROR: bad argument to --transitions.\n");
      lst_free(tmpl);
      break;
    case 'p':
      phi = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'z':
      zeta = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'E':
      omega = get_arg_dbl_bounds(optarg, 1, INFTY);
      mu = 1/omega;
      break;
    case 'C':
      gamma = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'r':
      refidx = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'M':
      refseq_f = fopen_fname(optarg, "r");
      break;
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == -1)
        die("ERROR: unrecognized alignment format.\n");
      break;
    case 'N':
      seqname_root = optarg;
      break;
    case 'P':
      idpref = optarg;
      break;
    case 'I':
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 4 && lst_size(tmpl) != 8)
	die("ERROR: bad argument to --indel-model.\n");
      alpha_n = lst_get_dbl(tmpl, 0);
      beta_n = lst_get_dbl(tmpl, 1);
      tau_n = lst_get_dbl(tmpl, 2);
      epsilon_n = lst_get_dbl(tmpl, 3);
      if (lst_size(tmpl) == 6) {
	alpha_c = lst_get_dbl(tmpl, 4);
	beta_c = lst_get_dbl(tmpl, 5);
	tau_c = lst_get_dbl(tmpl, 6);
	epsilon_c = lst_get_dbl(tmpl, 7);
      }
      else {
	alpha_c = alpha_n; beta_c = beta_n; tau_c = tau_n;
	epsilon_c = epsilon_n;
      }
      if (alpha_c <= 0 || alpha_c >= 1 || beta_c <= 0 || beta_c >= 1 ||
	  tau_c <= 0 || tau_c >= 1 || epsilon_c <= 0 || epsilon_c >= 1 ||
	  alpha_n <= 0 || alpha_n >= 1 || beta_n <= 0 || beta_n >= 1 ||
	  tau_n <= 0 || tau_n >= 1 || epsilon_n <= 0 || epsilon_n >= 1)
	die("ERROR: bad argument to --indel-model.\n");
      do_ih = TRUE;
      break;
    case 'L':
      msa_len = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'l':
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 3)
	die("ERROR: bad argument to --random-lengths.\n");
      min_len = (int)lst_get_dbl(tmpl, 0);
      max_len = (int)lst_get_dbl(tmpl, 1);
      lambda = lst_get_dbl(tmpl, 2);
      break;  
    case 'k':
      keep_ancestral = TRUE;
      break;
    case 'X':
      do_zeroed = TRUE;
      break;
    case 'x':
      cond_spec_f = fopen_fname(optarg, "r");
      do_zeroed = TRUE;
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'dmsimulate -h'.\n");
    }
  }

  if (optind != argc - 3) 
    die("Three arguments required.  Try 'dmsimulate -h'.\n");

  if (set_transitions && (gamma != -1 || omega != -1))
    die("ERROR: --transitions and --target-coverage/--expected-length cannot be used together.\n");

  if ((gamma != -1 && omega == -1) || (gamma == -1 && omega != -1))
    die("ERROR: --target-coverage and --expecteed-length must be used together.\n");

  if (gamma != -1)
    nu = gamma/(1-gamma) * mu;

  fprintf(stderr, "Reading tree model from %s...\n", argv[optind]);
  source_mod = tm_new_from_file(fopen_fname(argv[optind], "r"));

  fprintf(stderr, "Reading motif model from %s...\n", argv[optind+1]);
  motif_f = fopen_fname(argv[optind+1], "r");
  motif = mot_read(motif_f);

  if (source_mod->nratecats > 1) 
    die("ERROR: rate variation not currently supported.\n");

  if (source_mod->order > 0)
    die("ERROR: only single nucleotide models are currently supported.\n");

  if (!tm_is_reversible(source_mod->subst_mod))
    fprintf(stderr, "WARNING: p-value computation assumes reversibility and your model is non-reversible.\n");

  /* Name ancestral nodes */
  tr_name_ancestors(source_mod->tree);

  /* Instantiate the dmotif phmm */
  dm = dm_new(source_mod, motif, rho, mu, nu, phi, zeta, alpha_c, beta_c, 
              tau_c, epsilon_c, alpha_n, beta_n, tau_n, epsilon_n, FALSE,
	      FALSE, FALSE, FALSE);

  /* Read in the zeroed states if used and condition the model on site
     presence */
  if (do_zeroed) {
    zeroed_states = dms_read_zeroed_states(cond_spec_f);
    fclose(cond_spec_f);
    dms_condition_transitions(dm, zeroed_states);
  }
  
  /* set seqname and idpref, if necessary */
  if (seqname_root == NULL || idpref == NULL) {
    /* derive default from msa file name root */
    String *tmp = str_new_charstr(argv[optind+2]);
    if (!str_equals_charstr(tmp, "-")) {
      str_remove_path(tmp);
      if (idpref == NULL) idpref = strdup(tmp->chars);
      if (seqname_root == NULL) seqname_root = tmp->chars;    
    }
    else if (seqname_root == NULL) seqname_root = "refseq";
  }

  if (keep_ancestral == FALSE) {
    leaf_seqs = lst_new_int( ((source_mod->tree->nnodes - 1) / 2));
    for (i = 0; i < source_mod->tree->nnodes; i++) {
      n = lst_get_ptr(source_mod->tree->nodes, i);
      if (n->lchild == NULL) /* leaf node */
	lst_push_int(leaf_seqs, n->id);
    }			 
  }
  
  /* Simulate the alignment(s) */
  fprintf(stderr, "Simulating sequences...\n");
  /* Seed the random number generator */
  srandom(time(0));

  j = 1;
  msa_fname = smalloc(STR_MED_LEN * sizeof(char));
  seqname = smalloc(STR_MED_LEN * sizeof(char));
  motifs = gff_new_set();
  for (i = 0; i < msa_len; /* incrementing done inside loop */) {
    fprintf(stderr, "\tSimulating sequence %d...\n", j);

    if (lambda != -1) {
      len = ((int)gamma_draw(max_len, lambda) + min_len);
      sprintf(seqname, "%s.%d", seqname_root, j);
    } else {
      len = msa_len;
      sprintf(seqname, "%s", seqname_root);
    }
    i += len;

    path = smalloc(len * sizeof(int));
    msa = dm_generate_msa(len, dm, dm->phmm->mods, path);
    
    /*   for (i = 0; i < msa_len; i++) */
    /*     fprintf(stderr, "%d", path[i]); */
    /*   fprintf(stderr, "\n"); */

    if (do_ih) {
      ih = ih_new(source_mod->tree, len);
      fprintf(stderr, "\tSimulating indels for sequence %d...\n", j);
      indel_msa = dm_indel_mask(dm, msa, ih, path);
      msa_free(msa);
      if (!keep_ancestral) {
	fprintf(stderr, "\tPruning away ancestral sequences...\n");
	msa = msa_sub_alignment(indel_msa, leaf_seqs, TRUE, 0, len-1);
	msa_free(indel_msa);
	lst_free(leaf_seqs);
      } else {
	msa = indel_msa;
      }
    }
    
    /* Build the motif features GFF */
    tmp_gff = dm_labeling_as_gff(dm->phmm->cm, path, len, 
				 dm->m->width, dm->phmm->state_to_cat,
				 dm->state_to_motifpos, 
				 dm->phmm->reverse_compl,
				 seqname, "DMSIMULATE", NULL, NULL, idpref);
    
    for (k = 0; k < lst_size(tmp_gff->features); k++) {
      f = lst_get_ptr(tmp_gff->features, k);
      lst_push_ptr(motifs->features, gff_new_feature_copy(f));
    }
    
/*     fprintf(stderr, "i %d, j %d, len %d, tmp_gff %d, motifs %d\n", i, j, len,  */
/* 	    lst_size(tmp_gff->features), */
/* 	    lst_size(motifs->features)); */
    
    /* Print the alignment to the msa file */
    sprintf(msa_fname, "%s.%d.%s", argv[optind+2], j, 
	    msa_suffix_for_format(msa_format));
    fprintf(stderr, "\tWriting alignment to %s...\n", msa_fname);
    msa_f = fopen_fname(msa_fname, "w");
    msa_print(msa_f, msa, msa_format, 0);
    if (do_ih) {
      sprintf(msa_fname, "%s.%d.%s", argv[optind+2], j, "ih");
      fprintf(stderr, "\tWriting indel history to %s...\n", msa_fname);
      fclose(msa_f);
      msa_f = fopen_fname(msa_fname, "w");
      ih_print(ih, msa_f, seqname, "dmsimulate");
      ih_free(ih);
    }
    fclose(msa_f);
    free(path);
    msa_free(msa);
    gff_free_set(tmp_gff);
    j++;
  }

  /* Print motif features to stdout */
  fprintf(stderr, "Writing GFF to stdout...\n");
  gff_print_set(stdout, motifs);

  fprintf(stderr, "Done.\n");
  free(seqname);
  free(msa_fname);
  
  return 0;
}

