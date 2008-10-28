/* Dmotif implementation using a sampling strategy for both parameter
   estimation and path prediction. Written by Adam Diehl, Copyright 2008. */

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
#include <hashtable.h>
#include "dmsample.help"

#define DEFAULT_RHO 0.3
#define DEFAULT_PHI 0.5
#define DEFAULT_MU 0.01
#define DEFAULT_NU 0.01
#define DEFAULT_ZETA 0.001
#define DEFAULT_BSAMPLES 5000
#define DEFAULT_NSAMPLES 100000
#define DEFAULT_SAMPLE_INTERVAL 1

int main(int argc, char *argv[]) {
  char c, *key;
  int opt_idx, i, j, old_nnodes, **priors, *counts, cbstate;
  Multi_MSA *blocks;
  MSA *msa;
  List *pruned_names = lst_new_ptr(5), *tmpl;
  DMotifPhyloHmm *dm;
  GFF_Set *predictions, *reference = NULL;
  GFF_Feature *f;
  int found = FALSE;
  List *keys;
  PSSM *motif;
  double ***emissions;
  Hashtable *path_counts;
  String *cbname;

  struct option long_opts[] = {
    {"refseq", 1, 0, 'M'},
    {"refidx", 1, 0, 'r'},
    {"rho", 1, 0, 'R'},
    {"burn-in-samples", 1, 0, 'b'},
    {"samples", 1, 0, 's'},
    {"sample-interval", 1, 0, 'v'},
    {"seqname", 1, 0, 'N'},
    {"idpref", 1, 0, 'P'},
    {"indel-model", 1, 0, 'I'},
    {"log", 1, 0, 'l'},
    {"reference-gff", 1, 0, 'g'},
    {"ref-as-prior", 0, 0, 'u'},
    {"force_priors", 0, 0, 'p'},
    {"dump-hash", 1, 0, 'D'},
    {"precomputed-hash", 1, 0, 'd'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* arguments and defaults for options */
  FILE *refseq_f = NULL, *msa_f = NULL, *motif_f = NULL, *prior_f = NULL,
    *log = NULL, *ref_gff_f = NULL, *hash_f = NULL;
  TreeModel *source_mod;
  double rho = DEFAULT_RHO, mu = DEFAULT_MU, nu = DEFAULT_NU, 
    phi = DEFAULT_PHI, zeta = DEFAULT_ZETA,
    alpha_c = -1, beta_c = -1, tau_c = -1, epsilon_c = -1,
    alpha_n = -1, beta_n = -1, tau_n = -1, epsilon_n = -1;
  int refidx = 1, bsamples = DEFAULT_BSAMPLES, nsamples = DEFAULT_NSAMPLES, 
    sample_interval = DEFAULT_SAMPLE_INTERVAL, do_ih = 0, 
    ref_as_prior = FALSE, force_priors = FALSE, precomputed_hash = FALSE;
  char *seqname = NULL, *idpref = NULL;
  IndelHistory *ih = NULL;
  
  while ((c = getopt_long(argc, argv, "R:b:s:r:M:N:P:I:l:v:g:u:p:D:d:h",
			  long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'R':
      rho = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'b':
      bsamples = atoi(optarg);
      break;
    case 's':
      nsamples = atoi(optarg);
      break;
    case 'r':
      refidx = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'M':
      refseq_f = fopen_fname(optarg, "r");
      break;
    case 'N':
      seqname = optarg;
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
      do_ih = 1;
      break;
    case 'l':
      log = fopen_fname(optarg, "w");
      break;
    case 'v':
      sample_interval = get_arg_int_bounds(optarg, 0, INFTY);
    case 'g':
      ref_gff_f = fopen_fname(optarg, "r");
      fprintf(stderr, "Reading reference features from %s...\n", optarg);
      reference = gff_read_set(ref_gff_f);
      break;
    case 'u':
      ref_as_prior = TRUE;
      break;
    case 'p':
      force_priors = TRUE;
      break;
    case 'D':
      hash_f = fopen_fname(optarg, "w");
      break;
    case 'd':
      precomputed_hash = TRUE;
      hash_f = fopen_fname(optarg, "r");
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'dmsample -h'.\n");
    }
  }

  if (optind != argc - 4)
    die("Four arguments required.  Try 'dmsample -h'.\n");

  if (ref_gff_f != NULL && log == NULL) {
    if (!ref_as_prior)
      fprintf(stderr, "WARNING: Useless application of --reference-gff. Try 'dmsample -h'\n");
  }
  
  if (ref_as_prior && ref_gff_f == NULL)
    die("ERROR: --ref-as-prior requires --reference-gff. Try 'dmsample -h'\n");

  if (force_priors) {
    if (!ref_gff_f)
      die("ERROR: --force-priors requires --reference-gff. Try 'dmsample -h'\n");
    if (!ref_as_prior)
      ref_as_prior = TRUE;
  }

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

  /* Instead of reading in a sequence at a time, all sequences must be read in
     and stored in memory -- a path must be sampled for all alignments during
     each sampling iteration. Emissions only need to be computed once for all
     alignments, however, since we are not sampling rho. */

  /* read alignments */
  fprintf(stderr, "Reading alignments from %s...\n", argv[optind]);
  msa_f = fopen_fname(argv[optind], "r");
  blocks = msa_multimsa_new(msa_f, do_ih);

  fprintf(stderr, "Processing data in alignments...\n");
  for (i = 0; i < blocks->nblocks; i++) {
    if (msa_alph_has_lowercase(blocks->blocks[i]))
      msa_toupper(blocks->blocks[i]);
    msa_remove_N_from_alph(blocks->blocks[i]);

    if (blocks->blocks[i]->ss == NULL) {
      fprintf(stderr, "\tExtracting sufficient statistics for %s (%d of %d)...\n",
	      (((String*)lst_get(blocks->seqnames, i))->chars),
	      (i+1), blocks->nblocks);
      ss_from_msas(blocks->blocks[i], 1, TRUE, NULL, NULL, NULL, -1);
    }
    else if (msa->ss->tuple_idx == NULL)
      die("ERROR: ordered representation of alignment required unless --suff-stats.\n");
  }

  /* Read in priors for parameter estimation */
  fprintf(stderr, "Reading transition parameter priors from %s...\n", argv[optind + 3]);
  priors = smalloc(4 * sizeof(int*));
  for (i = 0; i < 4; i++) {
    priors[i] = smalloc(2 * sizeof(int));
    priors[i][0] = priors[i][1] = 0;
  }
  prior_f = fopen_fname(argv[optind + 3], "r");

  dms_read_priors(priors, prior_f);
  
  /* prune tree, if necessary */
  old_nnodes = source_mod->tree->nnodes;
  tm_prune(source_mod, blocks->blocks[0], pruned_names);

  if (lst_size(pruned_names) == (old_nnodes + 1) / 2)
    die("ERROR: no match for leaves of tree in alignment (leaf names must match alignment names).\n");
  if (lst_size(pruned_names) > 0) {
    fprintf(stderr, "WARNING: pruned away leaves of tree with no match in alignment (");
    for (i = 0; i < lst_size(pruned_names); i++)
      fprintf(stderr, "%s%s", ((String*)lst_get_ptr(pruned_names, i))->chars,
              i < lst_size(pruned_names) - 1 ? ", " : ").\n");
  }
  lst_free(pruned_names);

  /* this has to be done after pruning tree */
  tr_name_ancestors(source_mod->tree);

  /* also make sure match for reference sequence in tree */
  if (refidx > 0) {
    for (i = 0, found = FALSE; !found && i < source_mod->tree->nnodes; i++) {
      TreeNode *n = lst_get_ptr(source_mod->tree->nodes, i);
      if (!strcmp(n->name, blocks->blocks[0]->names[refidx-1]))
        found = TRUE;
    }
    if (!found) die("ERROR: no match for reference sequence in tree.\n");
  }
  
  dm = dm_new(source_mod, motif, rho, mu, nu, phi, zeta, alpha_c, beta_c,
              tau_c, epsilon_c, alpha_n, beta_n, tau_n, epsilon_n,
	      FALSE, FALSE, FALSE, FALSE);
  
  /* Cycle through msa's to compute emissions, adjust for missing data and
     adjust for indels */
  fprintf(stderr, "Computing emission probabilities...\n");

  emissions = smalloc(blocks->nblocks * sizeof(double*));
  /*  Needed for the emissions computation */
  dm->phmm->state_pos = smalloc(dm->phmm->nmods * sizeof(int));
  dm->phmm->state_neg = smalloc(dm->phmm->nmods * sizeof(int));

  for (i = 0; i < blocks->nblocks; i++) {

    fprintf(stderr, "\t%s (%d of %d)...\n",
	    (((String*)lst_get(blocks->seqnames, i))->chars),
	    (i+1), blocks->nblocks);

    msa = blocks->blocks[i];

    /* Allocate the emissions array slice */
    emissions[i] = smalloc(dm->phmm->hmm->nstates * sizeof(double*));
    for (j = 0; j < dm->phmm->hmm->nstates; j++)
      emissions[i][j] = smalloc(msa->length * sizeof(double));
    
    /* Have to hack the phmm so the internal emissions pointer references the
       current slice of the 3-dimensional emissions structure we need. */
    dm->phmm->emissions = emissions[i];

    fprintf(stderr, "\t\tComputing emissions.\n");
    phmm_compute_emissions(dm->phmm, msa, TRUE);
    
    /* add emissions for indel model, if necessary */
    if (do_ih) {
      fprintf(stderr, "\t\tAdjusting for indels.\n");
      ih = blocks->ih[i];
      dm_add_indel_emissions(dm, ih);
    }

    /* postprocess for missing data (requires special handling) */
    fprintf(stderr, "\t\tAdjusting for missing data.\n");
    dm_handle_missing_data(dm, msa);

    fprintf(stderr, "\t\tDone.\n");
  }
  
  /** Call the sampler **/
  if (!precomputed_hash) {
    fprintf(stderr, "Sampling state paths...\n");
    path_counts = hsh_new(max(blocks->nblocks, 10000));
    dms_sample_paths(dm, blocks, emissions, bsamples, nsamples, 
		     sample_interval, path_counts, priors, log, reference, 
		     ref_as_prior, force_priors);
  } else {
    fprintf(stderr, "Reading sampling data from disk...\n");
/*     fprintf(stderr, "nsamples_init %d\t", nsamples); */
    path_counts = dms_read_hash(hash_f, dm->phmm->hmm->nstates, &nsamples);
  }

/*   fprintf(stderr, "nsamples_new %d\n", nsamples); */
  
  /* Can free priors */
  for (i = 0; i < 4; i++)
    free(priors[i]);
  free(priors);

  /* Dump hash, for debugging purposes. */
  if (hash_f != NULL && !precomputed_hash) {
    dms_write_hash(path_counts, hash_f, dm->phmm->hmm->nstates, nsamples);
    return 0;
  }

  /* Generate a GFF from the features hash */
  fprintf(stderr, "Formatting output as GFF...\n");
  predictions = gff_new_set();
  cbname = str_new(STR_SHORT_LEN);
  str_append_charstr(cbname, "conserved-background");
  cbstate = cm_get_category(dm->phmm->cm, cbname);
  str_free(cbname);
  keys = hsh_keys(path_counts);
  for (i = 0; i < lst_size(keys); i++) {
    /* go through entries, build data for gff feture from each */
    key = lst_get_ptr(keys, i);
/*     fprintf(stderr, "i %d lst_size %d key %s\n", i, lst_size(keys), key); */
    counts = hsh_get(path_counts, key);
    f = dms_motif_as_gff_feat(dm, emissions, blocks, key, counts, nsamples,
			      sample_interval, cbstate);
    lst_push_ptr(predictions->features, f);
  }

  /* Free up some memory */
  lst_free(keys);  
  for (i = 0; i < blocks->nblocks; i++) {
/*     for (j = 0; j < dm->phmm->hmm->nstates; j++) */
/*       free(emissions[j]); */
    free(emissions[i]);
  }
  free(emissions);
  dm->phmm->emissions = NULL;
  
  /* convert GFF to coord frame of reference sequence and adjust
     coords by idx_offset, if necessary  */
  if (refidx != 0 || msa->idx_offset != 0)
    msa_map_gff_coords(msa, predictions, 0, refidx, msa->idx_offset, NULL);
  
  /* now output predictions */
  fprintf(stderr, "Writing GFF to stdout...\n");
  gff_print_set(stdout, predictions);

  if (reference != NULL)
    gff_free_set(reference);
  gff_free_set(predictions);
  
  fprintf(stderr, "Done.\n");
  
  return 0;
}

