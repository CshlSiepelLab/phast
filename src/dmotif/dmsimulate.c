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
#define DEFAULT_XI 0.0001
#define DEFAULT_MSA_LEN 1000000
#define DEFAULT_MMOD_TYPE "F81"

int main(int argc, char *argv[]) {
  char c, *msa_fname;
  int i, j, k, opt_idx, *path, msa_len, cat, branch;
  unsigned long int seed;
  FILE *devrandom;
  MSA *msa, *indel_msa;
  List *tmpl, *leaf_seqs, *zeroed_states, *motif_states;
  DMotifPhyloHmm *dm;
  GFF_Set *motifs, *tmp_gff;
  GFF_Feature *f;
  PSSM *motif;
  IndelHistory *ih;
  TreeNode *n;
  String *cname;

  struct option long_opts[] = {
    {"refseq", 1, 0, 'M'},
    {"msa-format", 1, 0, 'o'},
    {"refidx", 1, 0, 'r'},
    {"rho", 1, 0, 'R'},
    {"phi", 1, 0, 'p'},
    {"zeta", 1, 0, 'z'},
    {"xi", 1, 0, 'Z'},
    {"xi-off", 0, 0, 'F'},
    {"transitions", 1, 0, 't'},    
    {"expected-length", 1, 0, 'E'},
    {"target-coverage", 1, 0, 'C'},
    {"seqname", 1, 0, 'N'},
    {"idpref", 1, 0, 'P'},
    {"indel-model", 1, 0, 'I'},
    {"nc-mot-indel-mode", 0, 0, 'g'},
    {"msa-length", 1, 0, 'L'},
    {"random-lengths", 1, 0, 'l'},
    {"keep-ancestral", 0, 0, 'k'},
    {"cond-on-subs", 0, 0, 'X'},
    {"cond-on-species", 1, 0, 'x'},
    {"mot-mod-type", 1, 0, 'S'},
    {"scale-by-branch", 0, 0, 'B'},
    {"single-branch", 1, 0, 'b'},
    {"event", 1, 0, 'e'},
    {"nseqs", 1, 0, 'n'},
    {"no-require-subs", 0, 0, 's'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };
  
  /* arguments and defaults for options */
  FILE *refseq_f = NULL, *msa_f = NULL, *motif_f = NULL, *cond_spec_f = NULL;
  msa_format_type msa_format = FASTA;
  TreeModel *source_mod;
  double rho = DEFAULT_RHO, mu = DEFAULT_MU, nu = DEFAULT_NU, 
    phi = DEFAULT_PHI, zeta = DEFAULT_ZETA, xi = DEFAULT_XI,
    gamma = -1, omega = -1, 
    alpha_c = -1, beta_c = -1, tau_c = -1, epsilon_c = -1,
    alpha_n = -1, beta_n = -1, tau_n = -1, epsilon_n = -1,
    lambda = -1;
  int set_transitions = FALSE, refidx = 1, len = DEFAULT_MSA_LEN,
    max_len = 0, min_len = 0, do_ih = FALSE, keep_ancestral = FALSE, 
    do_zeroed = FALSE, xi_mode = TRUE, scale_by_branch = FALSE,
    nseqs = -1, fixed_path = FALSE, ncm_idl_mode = 0, require_subs = TRUE;
  char *seqname = NULL, *idpref = NULL, *seqname_root = NULL, *ename = NULL,
    *bname = NULL;
  subst_mod_type mmod_type = tm_get_subst_mod_type(DEFAULT_MMOD_TYPE);
  dmevent_t event;
  
  while ((c = getopt_long(argc, argv,
			  "R:t:p:z:Z:F:E:C:r:M:o:N:P:I:L:l:k:S:B:b:e:n:j:s:h", 
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
    case 'o':
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
      len = get_arg_int_bounds(optarg, 0, INFTY);
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
    case 'Z':
      xi = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 'F':
      xi_mode = FALSE;
      break;
    case 'S':
      mmod_type = tm_get_subst_mod_type(optarg);
      break;
    case 'B':
      scale_by_branch = TRUE;
      break;
    case 'b':
      bname = optarg;
      break;
    case 'e':
      ename = optarg;
      event = dm_get_event_type(ename);
      fixed_path = TRUE;
      break;
    case 'n':
      nseqs = atoi(optarg);
      break;
    case 'j':
      ncm_idl_mode = 1;
      break;
    case 's':
      require_subs = FALSE;
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
  lst_free(tmpl);

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

  if (ename != NULL) {
    if (event == -1)
      die("ERROR: --event requires a valid event type as an argument (CONS|NEUT|BIRTH|DEATH). Try dmsimulate -h.\n");
    if ((event == BIRTH || event == DEATH) && bname == NULL)
      die("ERROR: A branch must be specified when simulating sequences with the --event option and type BIRTH or DEATH.\n");
  }

  if (bname != NULL) {
    if (ename == NULL)
      die("ERROR: An event must be specified with --event when using --single-branch. Try dmsimulate -h.\n");
    if (event == CONS || event == NEUT)
      die("ERROR: --single-branch can only be used with --event BIRTH|DEATH. Try dmsimulate -h.\n");
  }

  if (nseqs != -1) {
    if (max_len != 0)
      die("ERROR: --nseqs cannot be used with --random-legths! Try dmsimulate -h.\n");
  }

  /* Name ancestral nodes */
  tr_name_ancestors(source_mod->tree);

  /* Instantiate the dmotif phmm */
  dm = dm_new(source_mod, motif, rho, mu, nu, phi, zeta, xi, xi_mode,
	      alpha_c, beta_c, tau_c, epsilon_c, alpha_n, beta_n, tau_n, 
	      epsilon_n, FALSE, FALSE, FALSE, FALSE, FALSE, mmod_type,
	      scale_by_branch, ncm_idl_mode);

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
    leaf_seqs = lst_new_int( ((source_mod->tree->nnodes + 1) / 2));
    for (i = 0; i < source_mod->tree->nnodes; i++) {
      n = lst_get_ptr(source_mod->tree->nodes, i);
      if (n->lchild == NULL) /* leaf node */
	lst_push_int(leaf_seqs, n->id);
    }
  }

  path = smalloc(len * sizeof(int));
  if (fixed_path == TRUE) { /* Simualting a single site type and fixed sequence
			       length */
    cname = str_new(STR_SHORT_LEN);
    if (bname != NULL && (event == BIRTH || event == DEATH)) {
      str_append_charstr(cname, event == BIRTH ? "birth" : "death");
      str_append_charstr(cname, "-");
      str_append_charstr(cname, bname);
    } else {
      str_append_charstr(cname, event == CONS ? "conserved" : "nonconserved");
    }
    str_append_charstr(cname, "-motif");
    cat = cm_get_category(dm->phmm->cm, cname);
    if (cat == 0)
      die("ERROR: Invalid branch given as argument to --branch!\n");
    branch = dm->state_to_branch[cat];

    motif_states = lst_new_int(dm->m->width);
    for (i = 0; i < lst_size(dm->branch_to_states[branch]); i++) {
      cat = lst_get_int(dm->branch_to_states[branch], i);
      if (dm->state_to_event[cat] == event &&
	  dm->state_to_motifpos[cat] != -1)
	lst_push_int(motif_states, cat);
    }

    j = 0;
    for (i = 0; i < len; i++) {
      if (i < ((len - dm->m->width) / 2) ||
	  i > ((len - dm->m->width) / 2) + (dm->m->width - 1)) {
	path[i] = 0; /* Neutral background states */
      } else { /* Motif states */
	path[i] = lst_get_int(motif_states, j);
	j++;
      }
/*       fprintf(stderr, "%d ", path[i]); */
    }
/*     fprintf(stderr, "\n"); */

    str_free(cname);
    lst_free(motif_states);
  }

  /* Simulate the alignment(s) */
  fprintf(stderr, "Simulating sequences...\n");
  /* Seed the random number generator */
  devrandom = fopen("/dev/random", "r");
  fread(&seed, sizeof(seed), 1, devrandom);
  srandom(abs(seed));
  fclose(devrandom);

  j = 1;
  if (nseqs != -1)
    msa_len = len * nseqs;
  else
    msa_len = len;
  
  msa_fname = smalloc(STR_MED_LEN * sizeof(char));
  seqname = smalloc(STR_MED_LEN * sizeof(char));
  motifs = gff_new_set();
  for (i = 0; i < msa_len; /* incrementing done inside loop */) {
    fprintf(stderr, "\tSimulating sequence %d...\n", j);
    
    if (nseqs != -1) {
      sprintf(seqname, "%s.%d", seqname_root, j);
    } else if (lambda != -1) {
      len = ((int)gamma_draw(max_len, lambda) + min_len);
      sprintf(seqname, "%s.%d", seqname_root, j);
    } else {
      sprintf(seqname, "%s", seqname_root);
    }
    i += len;

    msa = dm_generate_msa(len, dm, dm->phmm->mods, path, keep_ancestral,
			  fixed_path, require_subs);

    /*   for (i = 0; i < msa_len; i++) */
    /*     fprintf(stderr, "%d", path[i]); */
    /*   fprintf(stderr, "\n"); */

    if (do_ih) {
      ih = ih_new(source_mod->tree, len);
      fprintf(stderr, "\tSimulating indels for sequence %d...\n", j);
      indel_msa = dm_indel_mask(dm, msa, ih, path, FALSE);
      msa_free(msa);
      if (keep_ancestral == FALSE) {
	fprintf(stderr, "\tPruning away ancestral sequences...\n");
	msa = msa_sub_alignment(indel_msa, leaf_seqs, TRUE, 0, len-1);
	msa_free(indel_msa);
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
    msa_free(msa);
    gff_free_set(tmp_gff);
    j++;
  }

  /* Print motif features to stdout */
  fprintf(stderr, "Writing GFF to stdout...\n");
  gff_print_set(stdout, motifs);

  fprintf(stderr, "Done.\n");
  if (keep_ancestral == FALSE)
    lst_free(leaf_seqs);
  free(seqname);
  free(msa_fname);
  free(path);
  
  return 0;
}

