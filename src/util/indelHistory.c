/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: indelHistory.c,v 1.7 2008-11-12 02:07:58 acs Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <msa.h>
#include <tree_model.h>
#include <hashtable.h>
#include <sufficient_stats.h>
#include <indel_history.h>
#include "indelHistory.help"

int main(int argc, char *argv[]) {
  List *pruned_names = lst_new_ptr(5);
  TreeModel *source_mod;
  MSA *msa = NULL, *out_msa;
  IndelHistory *ih;
  char *read_hist_fname = NULL;
  char c;
  int opt_idx, old_nnodes, i;

  msa_format_type msa_format = UNKNOWN_FORMAT;
  int output_alignment = FALSE, ia_names = FALSE;
  
  struct option long_opts[] = {
    {"msa-format", 1, 0, 'i'},
    {"output-alignment", 0, 0, 'A'},
    {"read-history", 1, 0, 'H'},
    {"ia-names", 0, 0, 'I'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "i:H:AIh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == -1)
        die("ERROR: unrecognized alignment format.\n");
      break;
    case 'A':
      output_alignment = TRUE;
      break;
    case 'H':
      read_hist_fname = optarg;
      break;
    case 'I':
      ia_names = TRUE;
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'indelHistory -h'.\n");
    }
  }

  set_seed(-1);

  if (read_hist_fname != NULL) {
    fprintf(stderr, "Reading indel history from %s...\n", read_hist_fname);
    ih = ih_new_from_file(phast_fopen(read_hist_fname, "r"));
  }

  else {
    FILE *mfile;
    if (optind != argc - 2) 
      die("Two arguments required.  Try 'indelHistory -h'.\n");

    fprintf(stderr, "Reading alignment from %s...\n", argv[optind]);
    mfile = phast_fopen(argv[optind], "r");
    if (msa_format == UNKNOWN_FORMAT) 
      msa_format = msa_format_for_content(mfile, 1);
    msa = msa_new_from_file_define_format(mfile, msa_format, "ACGTNB^.-");
    phast_fclose(mfile);

    if (msa->seqs == NULL && (msa->ss == NULL || msa->ss->tuple_idx == NULL))
      die("ERROR: ordered representation of alignment required.\n");

    fprintf(stderr, "Reading tree from %s...\n", argv[optind+1]);
    source_mod = tm_new_from_file(phast_fopen(argv[optind+1], "r"), 1);
    
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
    
    tr_name_ancestors(source_mod->tree);

    if (msa->nseqs > (source_mod->tree->nnodes + 1) / 2) { /* assume ancestral 
							      seqs specified 
							      in this case */
      if (ia_names) {
        fprintf(stderr, "Converting sequence names...\n");
        ih_convert_ia_names(msa, source_mod->tree);
      }

      fprintf(stderr, "Extracting indel history from alignment...\n");
      ih = ih_extract_from_alignment(msa, source_mod->tree);
    }
   
    else {                        /* infer by parsimony */
      if (msa->ss == NULL) {
        fprintf(stderr, "Extracting sufficient statistics...\n");
        ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1, 0);
      }
      
      fprintf(stderr, "Inferring indel history by parsimony...\n");
      ih = ih_reconstruct(msa, source_mod->tree);
    }
  }

  if (output_alignment) {
    out_msa = ih_as_alignment(ih, msa);
    msa_print(stdout, out_msa, FASTA, FALSE);
  }
  else
    ih_print(ih, stdout, 
             read_hist_fname != NULL ? read_hist_fname : argv[optind], 
             "indelHistory");

  fprintf(stderr, "Done.\n");
  return 0;
}

