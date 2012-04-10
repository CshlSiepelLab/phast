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
#include <misc.h>
#include <pbs_code.h>
#include "pbsDecode.help"

int main(int argc, char *argv[]) {
  FILE *prob_f;
  char c;
  int opt_idx, i, nvals = 0, ngaps = 0, pos = 0;
  unsigned int idx;
  PbsCode *code;

  struct option long_opts[] = {
    {"start", 1, 0, 's'},
    {"end", 1, 0, 'e'},
    {"discard-gaps", 0, 0, 'G'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* options and defaults */
  int start = -1, end = -1, discard_gaps = FALSE;

  while ((c = getopt_long(argc, argv, "s:e:Gh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 's':
      start = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'e':
      end = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'G':
      discard_gaps = TRUE;
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'pbsDecode -h'.\n", argv[0]);
    }
  }

  if (optind != argc - 2) 
    die("Two arguments required.  Try 'pbsDecode -h'.\n");

  set_seed(-1);
    
  prob_f = phast_fopen(argv[optind], "rb");
  code = pbs_new_from_file(phast_fopen(argv[optind+1], "r"));

  if (start > 1) {
    if (fseek(prob_f, (start - 1) * code->nbytes, SEEK_SET) != 0)
      die("ERROR: fseek failed.\n");
    pos = start - 1;
  }

  while (pbs_read_binary(code, &idx, prob_f) != EOF) {
    
    pos++;

    if (end != -1 && pos > end) break;

    if (idx == code->gap_code) {
      if (!discard_gaps) 
	printf("-\n");
      ngaps++;
      continue;
    }

    if (idx >= code->code_size)
      die("ERROR: bad code ('%d')\n", idx); 

    for (i = 0; i < code->sg->d; i++)  
      printf("%f%s", code->rp[idx]->data[i],   
  	     i == code->sg->d - 1 ? "\n" : "\t");

    nvals++;
  }

  fprintf(stderr, "Dimensions: %d\n\
Rows per dimension: %d\n\
Code size: %d\n\
Bytes per vector: %d\n\
Vectors processed: %d\n\
Gaps: %d%s\n", 
	  code->sg->d, code->sg->nrows, code->code_size, 
	  code->nbytes, nvals, ngaps, discard_gaps ? " (discarded)" : "");
  
  return 0;
}
