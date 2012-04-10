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
#include <tree_model.h>
#include <prob_vector.h>
#include "modFreqs.help"

int main(int argc, char *argv[]) {
  char c;
  int opt_idx;
  Vector *newfreqs = vec_new(4);
  TreeModel *mod;

  struct option long_opts[] = {
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'modFreqs -h'.\n");
    }
  }

  if (optind != argc - 2 && optind != argc - 5) 
    die("ERROR: Wrong number of arguments.  Try 'modFreqs -h'.\n");

  set_seed(-1);

  mod = tm_new_from_file(phast_fopen(argv[optind], "r"), 1);

  if (!tm_is_reversible(mod)) 
    die("ERROR: reversible input model required.\n");
  if (mod->order != 0)
    die("ERROR: single nucleotide model required.\n");
  if (strcmp(mod->rate_matrix->states, DEFAULT_ALPHABET) != 0)
    die("ERROR: default nucleotide alphabet required.\n");

  if (optind == argc - 2) {
    double gc = get_arg_dbl_bounds(argv[optind+1], 0, 1);
    vec_set(newfreqs, 0, (1-gc)/2);
    vec_set(newfreqs, 1, gc/2);
    vec_set(newfreqs, 2, gc/2);
    vec_set(newfreqs, 3, (1-gc)/2);
  }
  else {
    vec_set(newfreqs, 0, get_arg_dbl_bounds(argv[optind+1], 0, 1));
    vec_set(newfreqs, 1, get_arg_dbl_bounds(argv[optind+2], 0, 1));
    vec_set(newfreqs, 2, get_arg_dbl_bounds(argv[optind+3], 0, 1));
    vec_set(newfreqs, 3, get_arg_dbl_bounds(argv[optind+4], 0, 1));
    pv_normalize(newfreqs);
  }

  tm_mod_freqs(mod, newfreqs);

  tm_print(stdout, mod);

  return 0;
}
