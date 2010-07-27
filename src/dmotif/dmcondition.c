/* Software to produce lists of states/positions in states to zero out when
   conditioning dmsample on presence of sites in a given species and/or
   requiring at least one substitution to make a gain or loss call. Works on
   a single MSA file and outputs a plain-text representation of a List of
   DMzeroedState objects.
   Written by Adam Diehl, Copyright 2009. */

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <getopt.h>
#include <misc.h>
#include <maf.h>
#include <sufficient_stats.h>
#include <phylo_hmm.h>
#include <dmotif_phmm.h>
#include <pssm.h>
#include "dmcondition.help"

#define DEFAULT_RHO 0.3
#define DEFAULT_PHI 0.5
#define DEFAULT_MU 0.01
#define DEFAULT_NU 0.01
#define DEFAULT_ZETA 0.001
#define DEFAULT_XI 0.0001

int main(int argc, char *argv[]) {
  int opt_idx, i, old_nnodes, found = FALSE;
  char c;
  MSA *msa;
  List *pruned_names = lst_new_ptr(5), *zeroed_states = NULL;
  DMzeroedState *z, **zeroed_states_ar = NULL;
  DMotifPhyloHmm *dm;
  TreeModel *source_mod;
  TreeNode *tree, *n;
  PSSM *motif;
 
  struct option long_opts[] = {
    {"refseq", 1, 0, 'M'},
    {"refidx", 1, 0, 'r'},
    {"cond-on-subs", 0, 0, 'X'},
    {"cond-on-species", 1, 0, 'x'},
    {"msa-format", 1, 0, 'i'},
    {"nosubs", 0, 0, 'n'},
    {"as-bed", 1, 0, 'b'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };
  
  /* arguments and defaults for options */
  FILE *refseq_f = NULL, *msa_f = NULL, *source_mod_f, *motif_f = NULL,
    *bedfile = NULL;
  int refidx = 1, cond_on_subs = FALSE, cond_spec = -1, nosubs = FALSE;
  String *cond_spec_str = NULL;
  msa_format_type msa_format = FASTA;  

#ifdef RPHAST
  GetRNGstate(); //seed R's random number generator
#endif

  while ((c = getopt_long(argc, argv, "M:r:X:x:i:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'M':
      refseq_f = fopen_fname(optarg, "r");
      break;
    case 'r':
      refidx = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'X':
      cond_on_subs = TRUE;
      break;
    case 'x':
      cond_spec_str = str_new_charstr(optarg);
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == -1)
        die("ERROR: unrecognized alignment format.\n");
      break;
    case 'n':
      nosubs = TRUE;
      break;
    case 'b':
      bedfile = fopen_fname(optarg, "w");
      break;
    case '?':
      die("Bad argument.  Try 'dmsample -h'.\n");
    }
  }

  /* Some sanity checks for proper usage */
  if (optind != argc - 3)
    die("Three arguments required.  Try 'dmcondition -h'.\n");

  if (cond_on_subs == FALSE && cond_spec_str == NULL)
    die("dmcondition requires either or both of the options --cond-on-subs or\n--cond-on-species! Try 'dmcondition -h'.\n");

  /* Handle arguments */
  msa_f = fopen_fname(argv[optind], "r");
  source_mod_f = fopen_fname(argv[optind+1], "r");
  motif_f = fopen_fname(argv[optind+2], "r");

  /* Read in the tree model and do some sanity checks */
  fprintf(stderr, "Reading tree model from %s...\n", argv[optind+1]);
  source_mod = tm_new_from_file(source_mod_f);
  
  /* Read in the motif model */
  fprintf(stderr, "Reading motif model from %s...\n", argv[optind+2]);
  motif = mot_read(motif_f);
  
  /* read alignment */
  fprintf(stderr, "Reading alignment from %s...\n", argv[optind]);
  if (msa_format == MAF) {
    if (refseq_f == NULL)
      die("ERROR: MAF input requires a reference sequence supplied with the\n--refseq option! Try 'dmcondition -h'.\n");
    msa = maf_read(msa_f, refseq_f, 1, NULL, NULL, NULL, -1, TRUE, NULL, 
                   NO_STRIP, FALSE);
    fclose(refseq_f);
  }
  else 
    msa = msa_new_from_file(msa_f, msa_format, NULL);
  
  if (msa_alph_has_lowercase(msa)) msa_toupper(msa); 
  msa_remove_N_from_alph(msa);
  
  if (msa->ss == NULL) {
    fprintf(stderr, "\tExtracting sufficient statistics...\n");
    ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1);
  }
  else if (msa->ss->tuple_idx == NULL)
    die("ERROR: ordered representation of alignment required!\n");
  
  fclose(msa_f);
  fclose(source_mod_f);
  fclose(motif_f);
  
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
  lst_free(pruned_names);
  
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
  
  dm = dm_new(source_mod, motif, 0.3, 0.01, 0.01, 0.01, 0.001, 0.0001, FALSE, 
	      0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
	      FALSE, FALSE, FALSE, FALSE, FALSE, F81, FALSE, 0);

  tree = source_mod->tree;

  /* Build the array of states and/or positions within states to zero */
/*   fprintf(stderr, "nstates %d, nnodes %d, cond_spec %d\n",  */
/* 	  dm->phmm->hmm->nstates, source_mod->tree->nnodes, cond_spec); */
  
  if (cond_spec_str != NULL) {
    fprintf(stderr, "Conditioning on site presence...\n");
    n = tr_get_node(tree, cond_spec_str->chars);
    cond_spec = n->id;
    zeroed_states_ar = dms_condition_on_species(dm, tree, cond_spec);    
  }

  if (cond_on_subs == TRUE) {
    fprintf(stderr, "Conditioning on substitutions...\n");
    if (zeroed_states_ar == NULL) {
      zeroed_states_ar = smalloc(dm->phmm->hmm->nstates *
				 sizeof(DMzeroedState*));
      for (i = 0; i < dm->phmm->hmm->nstates; i++)
	zeroed_states_ar[i] = NULL;
    }
    dms_condition_on_subs(dm, source_mod, msa, zeroed_states_ar, nosubs);
  }

  /* Convert the array to a list */
  zeroed_states = dms_zeroed_states_array_to_list(zeroed_states_ar, 
						  dm->phmm->hmm->nstates,
						  10000);

  /* Write the structure to stdout */
  fprintf(stderr, "Writing data to standard out...\n");
  dms_write_zeroed_states(stdout, zeroed_states);

  /* Write bedfile, if called for */
  if (bedfile != NULL) {
    fprintf(stderr, "Writing bed file...\n");
    dms_zeroed_as_bed(bedfile, dm, zeroed_states, argv[optind]);
  }
    
  /* Clean up after ourselves */
  dm_free(dm);
  for (i = 0; i < lst_size(zeroed_states); i++) {
    z = lst_get_ptr(zeroed_states, i);
    dms_free_zeroed_state(z);
  }
  if (zeroed_states != NULL)
    lst_free(zeroed_states);
  free(zeroed_states_ar);
  
  fprintf(stderr, "Done.\n");
  return 0;
}
  
  
