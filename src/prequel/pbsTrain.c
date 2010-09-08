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
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include <misc.h>
#include <stringsplus.h>
#include <pbs_code.h>
#include "pbsTrain.help"

int main(int argc, char *argv[]) {
  FILE *STATSF;
  char c;
  int opt_idx, i, max_nrows;
  String *line = str_new(STR_MED_LEN), *args = str_new(STR_MED_LEN);
  List *fields = lst_new_ptr(5), *vectors = lst_new_ptr(1000), 
    *counts = lst_new_int(1000);
  double dim = -1, error = -1;
  PbsCode *code;
  char comment[1000];
  time_t t;
  int have_data = TRUE;

  /* argument variables and defaults */
  int nrows = -1, nbytes = 1;
  training_mode mode = FULL;
  FILE *logf = NULL;
  struct timeval now;

  struct option long_opts[] = {
    {"nrows", 1, 0, 'n'},
    {"nbytes", 1, 0, 'b'},
    {"no-greedy", 0, 0, 'G'},
    {"no-train", 1, 0, 'x'},
    {"log", 1, 0, 'l'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  set_seed(-1);

  /* first capture arg list for comment in output */
  for (i = 1; i < argc; i++) {
    str_append_charstr(args, argv[i]);
    if (i < argc - 1) str_append_char(args, ' ');
  }

  while ((c = getopt_long(argc, argv, "n:b:l:Gxh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'n':
      nrows = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'b':
      nbytes = get_arg_int_bounds(optarg, 1, MAX_NBYTES);
      break;
    case 'G':
      mode = NO_GREEDY;
      break;
    case 'x':
      mode = NO_TRAIN;
      dim = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'l':
      logf = fopen_fname(optarg, "w+");
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'pbsTrain -h'.\n");
    }
  }

  if (mode == NO_TRAIN && optind == argc) 
    have_data = FALSE;		/* data optional */

  if (have_data) {
    if (optind != argc - 1) 
      die("ERROR: Bad arguments.  Try 'pbsTrain -h'.\n");

    STATSF = fopen_fname(argv[optind], "r");

    /* read stats */
    while (str_readline(line, STATSF) != EOF) {
      int count;
      double prob, norm_const;
      Vector *v;

      str_trim(line);
      if (line->length == 0 || line->chars[0] == '#') continue;

      str_split(line, NULL, fields);

      if (str_as_int(lst_get_ptr(fields, 0), &count) != 0)
	die("ERROR: Bad count in stats file ('%s')\n", lst_get_ptr(fields, 0));
      lst_push_int(counts, count);

      if (dim == -1) dim = lst_size(fields) - 1;
      else if (dim != lst_size(fields) - 1)
	die("ERROR: Each probability vector must have the same dimension\n");

      v = vec_new(dim);
      for (i = 0; i < dim; i++) {
	if (str_as_dbl(lst_get_ptr(fields, i+1), &prob) != 0 || 
	    prob < 0 || prob > 1)
	  die("ERROR: Bad probability in stats file ('%s')\n", 
	      lst_get_ptr(fields, i+1));

	vec_set(v, i, prob);
      }
    
      /* normalize to avoid problems from rounding errors */
      norm_const = normalize_probs(v->data, dim);
      if (fabs(1-norm_const) > 1e-2)
	die("ERROR: Probabilities in stats file don't sum to one.\nOffending line: '%s'\n",
	    line->chars);

      lst_push_ptr(vectors, v);
      lst_free_strings(fields);
    }
  }

  max_nrows = sxg_max_nrows(dim, ~(~0 << (8*nbytes)));
  if (nrows == -1) nrows = max_nrows;
  else if (nrows > max_nrows)
    die("ERROR: nrows exceeds maximum of %d for nbytes = %d and dimension = %d\n", 
	max_nrows, nbytes, dim);

  code = pbs_new(dim, nrows, nbytes);
  
  if (mode != NO_TRAIN)
    error = pbs_estimate_from_data(code, vectors, counts, logf, mode);

  else if (have_data) {		/* not training but need error */
    int tot_count = 0;
    double this_error;
    error = 0;
    for (i = 0; i < lst_size(vectors); i++) {
      pbs_get_index(code, lst_get_ptr(vectors, i), &this_error);
      error += this_error * lst_get_int(counts, i);
      tot_count += lst_get_int(counts, i);
    }    
    error /= tot_count;
  }

  /* generate comment */
  t = time(NULL);
  sprintf(comment, "# Code generated by pbsTrain, with argument(s) \"%s\"\n\
#  %s\n\
# Average training error = %f bits\n", 
	  args->chars, ctime(&t), error);

  pbs_write(code, stdout, comment);

  return 0;
}
