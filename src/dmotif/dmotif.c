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
#include "dmotif.help"

#define DEFAULT_RHO 0.3
#define DEFAULT_PHI 0.5
#define DEFAULT_MU 0.01
#define DEFAULT_NU 0.01
#define DEFAULT_ZETA 0.001

int main(int argc, char *argv[]) {
  char c;
  int opt_idx, i, old_nnodes;
  MSA *msa;
  List *pruned_names = lst_new_ptr(5), *tmpl;
  DMotifPhyloHmm *dm;
  GFF_Set *predictions;
  int found = FALSE;
  List *ignore_types = lst_new_ptr(1);
  PSSM *motif;

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
    {"indel-history", 1, 0, 'H'},
    {"score-list", 0, 0, 'l'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* arguments and defaults for options */
  FILE *refseq_f = NULL, *msa_f = NULL, *motif_f = NULL;
  msa_format_type msa_format = FASTA;
  TreeModel *source_mod;
  double rho = DEFAULT_RHO, mu = DEFAULT_MU, nu = DEFAULT_NU, 
    phi = DEFAULT_PHI, zeta = DEFAULT_ZETA, gamma = -1, omega = -1, 
    alpha_c = -1, beta_c = -1, tau_c = -1, epsilon_c = -1,
    alpha_n = -1, beta_n = -1, tau_n = -1, epsilon_n = -1;
  int set_transitions = FALSE, refidx = 1, estim_phi = TRUE, 
    estim_gamma = TRUE, estim_omega = TRUE, estim_zeta = TRUE,
    score_list = FALSE;
  char *seqname = NULL, *idpref = NULL;
  IndelHistory *ih = NULL;
  
  while ((c = getopt_long(argc, argv, "R:t:p:z:E:C:r:M:i:N:P:I:H:l:h", 
			  long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'R':
      rho = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 't':
      if (optarg[0] != '~') estim_gamma = estim_omega = FALSE;
      else optarg = &optarg[1];
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
      if (optarg[0] != '~') estim_phi = FALSE;
      else optarg = &optarg[1];
      phi = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'z':
      if (optarg[0] != '~') estim_zeta = FALSE;
      else optarg = &optarg[1];
      zeta = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'E':
      if (optarg[0] != '~') estim_omega = FALSE;
      else optarg = &optarg[1];
      omega = get_arg_dbl_bounds(optarg, 1, INFTY);
      mu = 1/omega;
      break;
    case 'C':
      if (optarg[0] != '~') estim_gamma = FALSE;
      else optarg = &optarg[1];
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
      seqname = optarg;
      break;
    case 'P':
      idpref = optarg;
      break;
    case 'I':
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 3 && lst_size(tmpl) != 6)
        die("ERROR: bad argument to --indel-model.\n");
      alpha_n = lst_get_dbl(tmpl, 0);
      beta_n = lst_get_dbl(tmpl, 1);
      tau_n = lst_get_dbl(tmpl, 2);
      if (lst_size(tmpl) == 6) {
        alpha_c = lst_get_dbl(tmpl, 3);
        beta_c = lst_get_dbl(tmpl, 4);
        tau_c = lst_get_dbl(tmpl, 5);
      }
      else {
        alpha_c = alpha_n; beta_c = beta_n; tau_c = tau_n;
      }
      if (alpha_c <= 0 || alpha_c >= 1 || beta_c <= 0 || beta_c >= 1 || 
          tau_c <= 0 || tau_c >= 1 || alpha_n <= 0 || alpha_n >= 1 || 
          beta_n <= 0 || beta_n >= 1 || tau_n <= 0 || tau_n >= 1)
        die("ERROR: bad argument to --indel-model.\n");
      break;
    case 'H':
      fprintf(stderr, "Reading indel history from %s...\n", optarg);
      ih = ih_new_from_file(fopen_fname(optarg, "r"));
      break;
    case 'l':
      score_list = TRUE;
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'dmotif -h'.\n");
    }
  }

  if (optind != argc - 3) 
    die("Three arguments required.  Try 'dmotif -h'.\n");

  if (set_transitions && (gamma != -1 || omega != -1))
    die("ERROR: --transitions and --target-coverage/--expected-length cannot be used together.\n");

  if ((gamma != -1 && omega == -1) || (gamma == -1 && omega != -1))
    die("ERROR: --target-coverage and --expecteed-length must be used together.\n");

  if (gamma != -1)
    nu = gamma/(1-gamma) * mu;

  fprintf(stderr, "Reading tree model from %s...\n", argv[optind+1]);
  source_mod = tm_new_from_file(fopen_fname(argv[optind+1], "r"));

  fprintf(stderr, "Reading motif model from %s...\n", argv[optind+2]);
  motif_f = fopen_fname(argv[optind+2], "r");
  motif = mot_read(motif_f);

  if (source_mod->nratecats > 1) 
    die("ERROR: rate variation not currently supported.\n");

  if (source_mod->order > 0)
    die("ERROR: only single nucleotide models are currently supported.\n");

  if (!tm_is_reversible(source_mod->subst_mod))
    fprintf(stderr, "WARNING: p-value computation assumes reversibility and your model is non-reversible.\n");

  /* read alignment */
  fprintf(stderr, "Reading alignment from %s...\n", argv[optind]);
  msa_f = fopen_fname(argv[optind], "r");
  if (msa_format == MAF) {
    msa = maf_read(msa_f, refseq_f, 1, NULL, NULL, NULL, -1, TRUE, NULL, 
                   NO_STRIP, FALSE); 
  }
  else 
    msa = msa_new_from_file(msa_f, msa_format, NULL);

  if (msa_alph_has_lowercase(msa)) msa_toupper(msa); 
  msa_remove_N_from_alph(msa);

  if (msa->ss == NULL) {
    fprintf(stderr, "Extracting sufficient statistics...\n");
    ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1);
  }
  else if (msa->ss->tuple_idx == NULL)
    die("ERROR: ordered representation of alignment required unless --suff-stats.\n");

  /* prune tree, if necessary */
  old_nnodes = source_mod->tree->nnodes;
  tm_prune(source_mod, msa, pruned_names);

  if (lst_size(pruned_names) == (old_nnodes + 1) / 2)
    die("ERROR: no match for leaves of tree in alignment (leaf names must match alignment names).\n");
  if (lst_size(pruned_names) > 0) {
    fprintf(stderr, "WARNING: pruned away leaves of tree with no match in alignment (");
    for (i = 0; i < lst_size(pruned_names); i++)
      fprintf(stderr, "%s%s", ((String*)lst_get_ptr(pruned_names, i))->chars, 
              i < lst_size(pruned_names) - 1 ? ", " : ").\n");
  }

  /* this has to be done after pruning tree */
  tr_name_ancestors(source_mod->tree);

  /* also make sure match for reference sequence in tree */
  if (refidx > 0) {
    for (i = 0, found = FALSE; !found && i < source_mod->tree->nnodes; i++) {
      TreeNode *n = lst_get_ptr(source_mod->tree->nodes, i);
      if (!strcmp(n->name, msa->names[refidx-1]))
        found = TRUE;
    }
    if (!found) die("ERROR: no match for reference sequence in tree.\n");
  }

  /* checks for indel model */
  if (alpha_c > 0) {
    if (ih == NULL) {
      fprintf(stderr, "Reconstructing indel history by parsimony...\n");
      ih = ih_reconstruct(msa, source_mod->tree);
    }
    else {
      if (ih->ncols != msa->length)
        die("ERROR: indel history doesn't seem to match alignment.\n");
      if (ih->tree->nnodes != source_mod->tree->nnodes)
        die("ERROR: indel history doesn't seem to match tree model.\n");
    }
  }

  dm = dm_new(source_mod, motif, rho, mu, nu, phi, zeta, alpha_c, beta_c, 
              tau_c, epsilon_c, alpha_n, beta_n, tau_n, epsilon_n, estim_gamma,
	      estim_omega, estim_phi, estim_zeta);

  /* compute emissions */
  phmm_compute_emissions(dm->phmm, msa, TRUE);

  /* If in score-list mode, just compute/print window scores and stop */
  if (score_list) {
    dm_print_motif_scores(dm);
    return 0;
  }

  /* add emissions for indel model, if necessary */
  if (alpha_c > 0) {
    fprintf(stderr, "Adjusting emissions for indels...\n");
    dm_add_indel_emissions(dm, ih);
  }

  /* postprocess for missing data (requires special handling) */
  fprintf(stderr, "Adjusting emissions for missing data...\n");
  dm_handle_missing_data(dm, msa); 

  if (estim_gamma || estim_omega || estim_phi) {
    fprintf(stderr, "Estimating free parameters...\n");
    dm_estimate_transitions(dm, msa);
  }

  /* set seqname and idpref, if necessary */
  if (seqname == NULL || idpref == NULL) {
    /* derive default from file name root */
    String *tmp = str_new_charstr(argv[optind]);
    if (!str_equals_charstr(tmp, "-")) {
      str_remove_path(tmp);
      str_root(tmp, '.');
      if (idpref == NULL) idpref = strdup(tmp->chars);
      if (seqname == NULL) seqname = tmp->chars;    
    }
    else if (seqname == NULL) seqname = "refseq";
  }

  /* obtain predictions */
  fprintf(stderr, "Running Viterbi algorithm...\n");
  predictions = dm_phmm_predict_viterbi(dm, seqname, NULL, idpref, NULL);
  lst_push_ptr(ignore_types, str_new_charstr("nonconserved"));
  gff_filter_by_type(predictions, ignore_types, TRUE, NULL);

  /* score predictions */
  fprintf(stderr, "Scoring predictions...\n");
  dm_score_predictions(dm, predictions);
  
  /* can free emissions now */
  for (i = 0; i < dm->phmm->hmm->nstates; i++)
    free(dm->phmm->emissions[i]);
  free(dm->phmm->emissions);
  dm->phmm->emissions = NULL;

  /* convert GFF to coord frame of reference sequence and adjust
     coords by idx_offset, if necessary  */
  if (refidx != 0 || msa->idx_offset != 0)
    msa_map_gff_coords(msa, predictions, 0, refidx, msa->idx_offset, NULL);

  if (refidx != 0) 
    gff_flatten(predictions);	
  /* necessary because coord conversion might create overlapping
     features (can happen in deletions in reference sequence) */

  /* now output predictions */
  fprintf(stderr, "Writing GFF to stdout...\n");
  gff_print_set(stdout, predictions);

  fprintf(stderr, "Done.\n");
  
  return 0;
}

