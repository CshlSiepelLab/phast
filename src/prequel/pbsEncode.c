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
#include "pbsEncode.help"

int main(int argc, char *argv[]) {
  FILE *prob_f;
  char c;
  int opt_idx, i, nlines = 0, ngaps = 0;
  unsigned idx;
  PbsCode *code;
  List *fields = lst_new_ptr(10);
  double error, tot_error = 0, prob;
  Vector *v;
  String *line = str_new(STR_MED_LEN);

  struct option long_opts[] = {
    {"discard-gaps", 0, 0, 'G'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* variables for options, with defaults */
  int discard_gaps = FALSE;

  set_seed(-1);

  while ((c = getopt_long(argc, argv, "Gh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'G':
      discard_gaps = TRUE;
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'pbsEncode -h'.\n");
    }
  }

  if (optind != argc - 2) 
    die("Two arguments required.  Try 'pbsEncode -h'.\n");
    
  prob_f = phast_fopen(argv[optind], "r");
  code = pbs_new_from_file(phast_fopen(argv[optind+1], "r"));
  v = vec_new(code->sg->d);
  
  while (str_readline(line, prob_f) != EOF) {

    if (line->length == 0 || line->chars[0] == '#') continue;

    str_split(line, NULL, fields);

    if (lst_size(fields) == 1 && str_equals_charstr(lst_get_ptr(fields, 0), "-")) {
      ngaps++;
      if (!discard_gaps)
	pbs_write_binary(code, code->gap_code, stdout);
    }
    else {			/* ordinary prob vector */
      if (lst_size(fields) != code->sg->d)
	die("ERROR: number of columns must equal dimension of code (%d).\n",
	    code->sg->d);
    
      for (i = 0; i < code->sg->d; i++) {
	if (str_as_dbl(lst_get_ptr(fields, i), &prob) != 0 || prob < 0 || prob > 1)
	  die("ERROR: bad value ('%s')\n", lst_get_ptr(fields, i));

	vec_set(v, i, prob);
      }

      idx = pbs_get_index(code, v, &error);
      tot_error += error;
      pbs_write_binary(code, idx, stdout);
      nlines++;
    }

    lst_free_strings(fields);
  }

  fprintf(stderr, "Dimensions: %d\n\
Rows per dimension: %d\n\
Code size: %d\n\
Bytes per vector: %d\n\
Vectors processed: %d\n\
Gaps: %d%s\n\
Average approximation error: %f bits\n", 
	  code->sg->d, code->sg->nrows, code->code_size, code->nbytes,  
	  nlines, ngaps, discard_gaps ? " (discarded)" : "", tot_error/nlines);

  return 0;
}
