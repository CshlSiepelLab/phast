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
#include <indel_mod.h>
#include <subst_distrib.h>
#include <bd_phylo_hmm.h>
#include "dless.help"

#define DEFAULT_RHO 0.3
#define DEFAULT_PHI 0.5
#define DEFAULT_MU 0.01
#define DEFAULT_NU 0.01

int main(int argc, char *argv[]) {
  char c;
  char *msa_fname = NULL;
  int opt_idx, i, old_nnodes;
  MSA *msa;
  List *pruned_names = lst_new_ptr(5), *tmpl;
  BDPhyloHmm *bdphmm;
  GFF_Set *predictions;
  int found = FALSE;
  List *ignore_types = lst_new_ptr(1);

  struct option long_opts[] = {
    {"refseq", 1, 0, 'M'},
    {"msa-format", 1, 0, 'i'},
    {"refidx", 1, 0, 'r'},
    {"rho", 1, 0, 'R'},
    {"phi", 1, 0, 'p'},
    {"transitions", 1, 0, 't'},    
    {"expected-length", 1, 0, 'E'},
    {"target-coverage", 1, 0, 'C'},
    {"seqname", 1, 0, 'N'},
    {"idpref", 1, 0, 'P'},
    {"indel-model", 1, 0, 'I'},
    {"indel-history", 1, 0, 'H'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* arguments and defaults for options */
  FILE *refseq_f = NULL, *msa_f = NULL;
  msa_format_type msa_format = UNKNOWN_FORMAT;
  TreeModel *source_mod;
  double rho = DEFAULT_RHO, mu = DEFAULT_MU, nu = DEFAULT_NU, 
    phi = DEFAULT_PHI, gamma = -1, omega = -1, 
    alpha_c = -1, beta_c = -1, tau_c = -1,
    alpha_n = -1, beta_n = -1, tau_n = -1;
  int set_transitions = FALSE, refidx = 1, estim_phi = TRUE, 
    estim_gamma = TRUE, estim_omega = TRUE;
  char *seqname = NULL, *idpref = NULL;
  IndelHistory *ih = NULL;

  while ((c = getopt_long(argc, argv, "R:t:p:E:C:r:M:i:N:P:I:H:h", long_opts, &opt_idx)) != -1) {
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
      refseq_f = phast_fopen(optarg, "r");
      break;
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == UNKNOWN_FORMAT)
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
      ih = ih_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'dless -h'.\n");
    }
  }

  if (optind != argc - 1)
    die("Missing alignment file or model file.  Try 'dless -h'.\n");

  if (set_transitions && (gamma != -1 || omega != -1))
    die("ERROR: --transitions and --target-coverage/--expected-length cannot be used together.\n");

  if ((gamma != -1 && omega == -1) || (gamma == -1 && omega != -1))
    die("ERROR: --target-coverage and --expecteed-length must be used together.\n");

  set_seed(-1);

  if (gamma != -1)
    nu = gamma/(1-gamma) * mu;

  fprintf(stderr, "Reading tree model from %s...\n", argv[optind]);
  source_mod = tm_new_from_file(phast_fopen(argv[optind], "r"), 1);

  if (source_mod->nratecats > 1) 
    die("ERROR: rate variation not currently supported.\n");

  if (source_mod->order > 0)
    die("ERROR: only single nucleotide models are currently supported.\n");

  if (!tm_is_reversible(source_mod))
    phast_warning("WARNING: p-value computation assumes reversibility and your model is non-reversible.\n");

  /* read alignment */
  msa_f = phast_fopen(argv[optind], "r");

  fprintf(stderr, "Reading alignment from %s...\n", argv[optind]);
  if (msa_format == UNKNOWN_FORMAT) 
    msa_format = msa_format_for_content(msa_f, 1);

  if (msa_format == MAF) {
    msa = maf_read(msa_f, refseq_f, 1, NULL, NULL, NULL, -1, TRUE, NULL, 
                   NO_STRIP, FALSE); 
  }
  else 
    msa = msa_new_from_file_define_format(msa_f, msa_format, NULL);

  if (msa_alph_has_lowercase(msa)) msa_toupper(msa); 
  msa_remove_N_from_alph(msa);

  if (msa->ss == NULL) {
    fprintf(stderr, "Extracting sufficient statistics...\n");
    ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1, 0);
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

  bdphmm = bd_new(source_mod, rho, mu, nu, phi, alpha_c, beta_c, tau_c, 
                  alpha_n, beta_n, tau_n, estim_gamma, estim_omega, 
                  estim_phi);

  /* compute emissions */
  phmm_compute_emissions(bdphmm->phmm, msa, FALSE);

  /* add emissions for indel model, if necessary */
  if (alpha_c > 0) {
    fprintf(stderr, "Adjusting emissions for indels...\n");
    bd_add_indel_emissions(bdphmm, ih);
  }

  /* postprocess for missing data (requires special handling) */
  fprintf(stderr, "Adjusting emissions for missing data...\n");
  bd_handle_missing_data(bdphmm, msa);

  if (estim_gamma || estim_omega || estim_phi) {
    fprintf(stderr, "Estimating free parameters...\n");
    bd_estimate_transitions(bdphmm, msa);
  }

  /* set seqname and idpref, if necessary */
  if (seqname == NULL || idpref == NULL) {
    /* derive default from file name root */
    String *tmp = str_new_charstr(msa_fname);
    if (!str_equals_charstr(tmp, "-")) {
      str_remove_path(tmp);
      str_root(tmp, '.');
      if (idpref == NULL) idpref = copy_charstr(tmp->chars);
      str_root(tmp, '.');         /* apply one more time for double suffix */
      if (seqname == NULL) seqname = tmp->chars;    
    }
    else if (seqname == NULL) seqname = "refseq";
  }

  /* obtain predictions */
  fprintf(stderr, "Running Viterbi algorithm...\n");
  predictions = phmm_predict_viterbi(bdphmm->phmm, seqname, NULL, idpref, NULL);
  lst_push_ptr(ignore_types, str_new_charstr("nonconserved"));
  gff_filter_by_type(predictions, ignore_types, TRUE, NULL);

  /* score predictions */
  fprintf(stderr, "Scoring predictions...\n");
  bd_score_predictions(bdphmm, predictions);
  
  /* can free emissions now */
  for (i = 0; i < bdphmm->phmm->hmm->nstates; i++)
    sfree(bdphmm->phmm->emissions[i]);
  sfree(bdphmm->phmm->emissions);
  bdphmm->phmm->emissions = NULL;

  /* convert GFF to coord frame of reference sequence and adjust
     coords by idx_offset, if necessary  */
  if (refidx != 0 || msa->idx_offset != 0)
    msa_map_gff_coords(msa, predictions, 0, refidx, msa->idx_offset);

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

