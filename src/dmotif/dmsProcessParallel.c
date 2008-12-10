/* Post-processor for dmsample hashes computed by parallel markov chains. 
   Written by Adam Diehl, Copyright 2008. */

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
#define DEFAULT_SAMPLE_INTERVAL 1
#define MEDIUM 25

int main(int argc, char *argv[]) {
  char c, *key;
  int opt_idx, i, old_nnodes, *counts, nsamples, nsamples_this,
    max_seqlen;
/*   int cbstate, j; */
  DMotifPmsaStruct *dmpmsa;
  PooledMSA *blocks;
/*   MSA *msa; */
  List *pruned_names = lst_new_ptr(5), *tmpl;
  DMotifPhyloHmm *dm;
  GFF_Set *predictions;
  GFF_Feature *f;
  int found = FALSE;
  List *keys, *seqnames;
  PSSM *motif;
  Hashtable *path_counts, *tmp;
/*   String *cbname; */
  FILE *hash_f;

  struct option long_opts[] = {
    {"refseq", 1, 0, 'M'},
    {"refidx", 1, 0, 'r'},
    {"rho", 1, 0, 'R'},
    {"seqname", 1, 0, 'N'},
    {"idpref", 1, 0, 'P'},
    {"indel-model", 1, 0, 'I'},
    {"dump-hash", 1, 0, 'D'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* arguments and defaults for options */
  FILE *refseq_f = NULL, *msa_f = NULL, *motif_f = NULL, *dump_f = NULL;
  TreeModel *source_mod;
  double rho = DEFAULT_RHO, mu = DEFAULT_MU, nu = DEFAULT_NU, 
    phi = DEFAULT_PHI, zeta = DEFAULT_ZETA,
    alpha_c = -1, beta_c = -1, tau_c = -1, epsilon_c = -1,
    alpha_n = -1, beta_n = -1, tau_n = -1, epsilon_n = -1;
  int refidx = 1, sample_interval = DEFAULT_SAMPLE_INTERVAL, do_ih = 0;
  char *seqname = NULL, *idpref = NULL, delim[2];
/*   IndelHistory *ih = NULL; */
  String *hash_f_string;
  List *hash_files;  

  while ((c = getopt_long(argc, argv, "R:r:M:N:P:I:D:h",
			  long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'R':
      rho = get_arg_dbl_bounds(optarg, 0, 1);
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
    case 'D':
      dump_f = fopen_fname(optarg, "w");
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'dmsample -h'.\n");
    }
  }

  if (optind != argc - 4)
    die("Four arguments required.  Try 'dmsProcessParallel -h'.\n");

  /* Load up the list of hash files */
  hash_f_string = str_new_charstr(argv[optind]);
  hash_files = lst_new_ptr(MEDIUM);
  sprintf(delim, ",");
  str_split(hash_f_string, delim, hash_files);
/*   fprintf(stderr, "%s, %d\n", hash_f_string->chars, lst_size(hash_files)); */
/*   for (i = 0; i < lst_size(hash_files); i++) */
  /*     fprintf(stderr, "%s\n", ((String*)lst_get_ptr(hash_files, i))->chars); */
  str_free(hash_f_string);

  /* read alignments -- we really only need sequence lengths! */
  fprintf(stderr, "Reading alignments from %s...\n", argv[optind+1]);
  msa_f = fopen_fname(argv[optind+1], "r");
  dmpmsa = dms_read_alignments(msa_f, do_ih);
  blocks = dmpmsa->pmsa;
  seqnames = dmpmsa->seqnames;
  max_seqlen = dmpmsa->max_seqlen;

  /* Read phylogenetic model file */
  fprintf(stderr, "Reading tree model from %s...\n", argv[optind+2]);
  source_mod = tm_new_from_file(fopen_fname(argv[optind+2], "r"));

  fprintf(stderr, "Reading motif model from %s...\n", argv[optind+3]);
  motif_f = fopen_fname(argv[optind+3], "r");
  motif = mot_read(motif_f);

  if (source_mod->nratecats > 1)
    die("ERROR: rate variation not currently supported.\n");

  if (source_mod->order > 0)
    die("ERROR: only single nucleotide models are currently supported.\n");

  if (!tm_is_reversible(source_mod->subst_mod))
    fprintf(stderr, "WARNING: p-value computation assumes reversibility and your model is non-reversible.\n");

  /* prune tree, if necessary */
  old_nnodes = source_mod->tree->nnodes;
  tm_prune(source_mod, lst_get_ptr(blocks->source_msas, 0), pruned_names);

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
      if (!strcmp(n->name, 
		  ((MSA*)lst_get_ptr(blocks->source_msas, 0))->names[refidx-1]))
        found = TRUE;
    }
    if (!found) die("ERROR: no match for reference sequence in tree.\n");
  }
  
  dm = dm_new(source_mod, motif, rho, mu, nu, phi, zeta, alpha_c, beta_c,
              tau_c, epsilon_c, alpha_n, beta_n, tau_n, epsilon_n,
	      FALSE, FALSE, FALSE, FALSE);

  /* Not using emissions, so no need to compute */

  /* Read hashes from each file, flattening into the merged hash as we go */
  path_counts = hsh_new(lst_size(hash_files));
  for (i = 0; i < lst_size(hash_files); i++) {
    hash_f = fopen_fname(((String*)lst_get_ptr(hash_files, i))->chars, "r");
    fprintf(stderr, "Reading sampling data from file %s...\n", 
	    ((String*)lst_get_ptr(hash_files, i))->chars);
/*     fprintf(stderr, "nsamples_init %d\t", nsamples); */
    tmp = dms_read_hash(hash_f, dm->phmm->hmm->nstates, &nsamples_this);
    nsamples += nsamples_this;
/*     fprintf(stderr, "nsamples_this = %d, nsamples = %d\n", nsamples_this, nsamples); */
    dms_combine_hashes(path_counts, tmp, dm->phmm->hmm->nstates);
    fclose(hash_f);
    hsh_free(tmp);
  }

  /* Dump hash, for debugging purposes. */
  if (dump_f != NULL) {
    dms_write_hash(path_counts, dump_f, dm->phmm->hmm->nstates, nsamples);
    return 0;
  }
    
  /* Generate a GFF from the features hash */
  fprintf(stderr, "Formatting output as GFF...\n");
  predictions = gff_new_set();
  keys = hsh_keys(path_counts);
  for (i = 0; i < lst_size(keys); i++) {
    /* go through entries, build data for gff feture from each */
    key = lst_get_ptr(keys, i);
/*     fprintf(stderr, "i %d lst_size %d key %s\n", i, lst_size(keys), key); */
    counts = hsh_get(path_counts, key);
    f = dms_motif_as_gff_feat(dm, blocks, seqnames, key, counts, nsamples,
			      sample_interval);
    lst_push_ptr(predictions->features, f);
  }

  /* Free up some memory */
  lst_free(keys);  
    
  /* now output predictions */
  fprintf(stderr, "Writing GFF to stdout...\n");
  gff_print_set(stdout, predictions);

  gff_free_set(predictions);
  
  fprintf(stderr, "Done.\n");
  
  return 0;
}

