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
#include <getopt.h>
#include <phast/misc.h>
#include <phast/msa.h>
#include <phast/maf.h>
#include <phast/nj.h>
#include "nj_var.help"

int main(int argc, char *argv[]) {
  signed char c;
  int opt_idx, i, ntips = 0;
  MSA *msa = NULL;

  char *alphabet = "ACGT";
  char **names = NULL;
  msa_format_type format = UNKNOWN_FORMAT;
  FILE *infile = NULL, *indistfile = NULL, *outdistfile = NULL;
  Matrix *D;
  TreeNode *tree;
  List *namestr;

  struct option long_opts[] = {
    {"format", 1, 0, 'i'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "i:d:o:n:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'i':
      format = msa_str_to_format(optarg);
      if (format == UNKNOWN_FORMAT) 
        die("ERROR: bad input format.\n");
      break;
    case 'd':
      indistfile = phast_fopen(optarg, "r");
      break;
    case 'o':
      outdistfile = phast_fopen(optarg, "w");
      break;
    case 'n':
      namestr = get_arg_list(optarg);
      ntips = lst_size(namestr);
      names = smalloc(ntips * sizeof(char*));
      for (i = 0; i < ntips; i++)
	names[i] = ((String*)lst_get_ptr(namestr, i))->chars;
      break;
    case 'h':
      printf("%s", HELP); 
      exit(0);
    case '?':
      die("Bad argument.  Try 'nj_var -h'.\n");
    }
  }

  if ((indistfile == NULL && optind != argc - 1) ||
      (indistfile != NULL && optind != argc))
    die("Bad arguments.  Try 'nj_var -h'.\n");

  if (indistfile == NULL) {
    infile = phast_fopen(argv[optind], "r");
    if (format == UNKNOWN_FORMAT)
      format = msa_format_for_content(infile, 1);
    if (format == MAF) 
      msa = maf_read(phast_fopen(argv[optind], "r"), NULL, 1, alphabet,
		     NULL, NULL, -1, TRUE, NULL, NO_STRIP, FALSE);
    else
      msa = msa_new_from_file_define_format(phast_fopen(argv[optind], "r"), 
					    format, alphabet);

    D = nj_compute_JC_matr(msa);
    names = msa->names;
  }
  
  else {
    if (names == NULL || ntips <= 0)
      die("Bad arguments.  Try 'nj_var -h'.\n");
    D = mat_new_from_file(indistfile, ntips, ntips);
  }

  if (outdistfile != NULL)
    mat_print(D, outdistfile);
  
  tree = nj_infer_tree(D, names);
  tr_print(stdout, tree, TRUE);
  
  return (0);
}
