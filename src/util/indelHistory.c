/* $Id: indelHistory.c,v 1.5 2005-09-04 21:28:22 acs Exp $
   Written by Adam Siepel, 2005
   Copyright 2005, Adam Siepel, University of California */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <msa.h>
#include <trees.h>
#include <hashtable.h>
#include <sufficient_stats.h>
#include <indel_history.h>
#include "indelHistory.help"

int main(int argc, char *argv[]) {
  TreeNode *tree;
  MSA *msa = NULL, *out_msa;
  IndelHistory *ih;
  char *read_hist_fname = NULL;
  char c;
  int opt_idx;

  msa_format_type msa_format = FASTA;
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
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'indelHistory -h'.\n");
    }
  }

  if (read_hist_fname != NULL) {
    fprintf(stderr, "Reading indel history from %s...\n", read_hist_fname);
    ih = ih_new_from_file(fopen_fname(read_hist_fname, "r"));
  }

  else {
    if (optind != argc - 2) 
      die("Two arguments required.  Try 'indelHistory -h'.\n");

    fprintf(stderr, "Reading alignment from %s...\n", argv[optind]);
    msa = msa_new_from_file(fopen_fname(argv[optind], "r"), msa_format, "ACGTNB^.-");

    if (msa->seqs == NULL && (msa->ss == NULL || msa->ss->tuple_idx == NULL))
      die("ERROR: ordered representation of alignment required.\n");

    fprintf(stderr, "Reading tree from %s...\n", argv[optind+1]);
    tree = tr_new_from_file(fopen_fname(argv[optind+1], "r"));

    tr_name_ancestors(tree);

    if (msa->nseqs > (tree->nnodes + 1) / 2) { /* assume ancestral seqs
                                                  specified in this case */
      if (ia_names) {
        fprintf(stderr, "Converting sequence names...\n");
        ih_convert_ia_names(msa, tree);
      }

      fprintf(stderr, "Extracting indel history from alignment...\n");
      ih = ih_extract_from_alignment(msa, tree);
    }
   
    else {                        /* infer by parsimony */
      if (msa->ss == NULL) {
        fprintf(stderr, "Extracting sufficient statistics...\n");
        ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1);
      }
      
      fprintf(stderr, "Inferring indel history by parsimony...\n");
      ih = ih_reconstruct(msa, tree);
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

