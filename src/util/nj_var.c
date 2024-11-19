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
  int opt_idx;
  MSA *msa;

  char *alphabet = "ACGT";
  msa_format_type format = UNKNOWN_FORMAT;
  FILE *infile;
  Matrix *D;
  TreeNode *tree;

  struct option long_opts[] = {
    {"format", 1, 0, 'i'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "i:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'i':
      format = msa_str_to_format(optarg);
      if (format == UNKNOWN_FORMAT) 
        die("ERROR: bad input format.\n");
      break;      
    case 'h':
      printf("%s", HELP); 
      exit(0);
    case '?':
      die("Bad argument.  Try 'program -h'.\n");
    }
  }

  if (optind != argc - 1) 
    die("Alignment filename required.  Try 'nj_var -h'.\n");


  infile = phast_fopen(argv[optind], "r");
  if (format == UNKNOWN_FORMAT)
    format = msa_format_for_content(infile, 1);
  if (format == MAF) 
    msa = maf_read(phast_fopen(argv[optind], "r"), NULL, 1, alphabet,
		   NULL, NULL, -1, TRUE, NULL, NO_STRIP, FALSE);
  else
    msa = msa_new_from_file_define_format(phast_fopen(argv[optind], "r"), 
                             format, alphabet);


  /* FIXME
     allow a distance matrix to be specified directly
     option to output mod file?  maybe not really necessary
   */


  D = nj_compute_JC_matr(msa);

  /*  mat_print(D, stdout);*/

  
  tree = nj_infer_tree(D, msa->names);
  tr_print(stdout, tree, TRUE);
  
  return (0);
}
